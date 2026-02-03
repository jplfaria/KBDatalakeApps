"""
BERDL API implementation of QueryPangenomeABC.

This module provides pangenome queries against the BERDL Data Lake API.
It implements the same interface as QueryPangenomeLocal but retrieves data
from the remote API instead of local Parquet/CSV files.

Author: Jose P. Faria (jplfaria@gmail.com)
Date: February 2026

=== KNOWN ISSUES & WORKAROUNDS ===

Issue 1: gtdb_species_clade_id with '--' chars causes 500 errors with '=' operator
    Workaround: Use LIKE with species name pattern instead
    Example: 's__Staphylococcus_lugdunensis--RS_GCF_002901705.1'
          -> LIKE '%Staphylococcus_lugdunensis%'

Issue 2: gene_genecluster_junction table (1B rows) may timeout
    Workaround: Use LIMIT 100 instead of 1000 for pagination
    Note: This table may still be unreliable - use with caution
"""

import time
import requests
import pandas as pd
from typing import List, Optional, Any

from .query_pangenome import QueryPangenomeABC


class QueryPangenomeBERDL(QueryPangenomeABC):
    """
    Query KBase pangenome data from BERDL Data Lake API.

    Schema: kbase_ke_pangenome
    Tables:
        - genome (293K rows)
        - gene_cluster (132M rows)
        - genome_ani (421M rows)
        - gene_genecluster_junction (1B rows) - may timeout
        - gtdb_species_clade (27K rows)
        - ncbi_env (4M rows)
    """

    BERDL_API_URL = "https://hub.berdl.kbase.us/apis/mcp/delta/tables/query"
    SCHEMA = "kbase_ke_pangenome"

    # API settings
    PAGE_LIMIT = 1000
    JUNCTION_PAGE_LIMIT = 100  # Smaller limit for 1B row junction table
    MAX_RETRIES = 3
    RETRY_DELAY = 2
    TIMEOUT = 120

    def __init__(self, token: str):
        """
        Initialize the BERDL pangenome query client.

        Args:
            token: BERDL API token (from KBase environment or MFA login).
        """
        if not token:
            raise ValueError("Token is required for BERDL API access.")

        self.token = token
        self.headers = {
            "accept": "application/json",
            "Authorization": f"Bearer {token}",
            "Content-Type": "application/json",
            "User-Agent": "Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36"
        }

    @staticmethod
    def clade_id_to_like_pattern(clade_id: str) -> str:
        """
        Convert clade_id to a LIKE pattern that works with BERDL API.

        Problem: Full clade IDs with '--' cause 500 errors with '=' operator
        Solution: Extract species name and use LIKE pattern

        Example:
            Input:  's__Staphylococcus_lugdunensis--RS_GCF_002901705.1'
            Output: '%Staphylococcus_lugdunensis%'
        """
        if '--' in clade_id:
            species_part = clade_id.split('--')[0]
            if species_part.startswith('s__'):
                species_part = species_part[3:]
            return f'%{species_part}%'
        return f'%{clade_id}%'

    def _execute_query(self, query: str, limit: int = 10000) -> pd.DataFrame:
        """
        Execute a SQL query against the BERDL API with pagination and retry logic.
        """
        all_results = []
        offset = 0

        while len(all_results) < limit:
            page_limit = min(self.PAGE_LIMIT, limit - len(all_results))
            payload = {
                "query": query,
                "limit": page_limit,
                "offset": offset
            }

            retries = 0
            while retries <= self.MAX_RETRIES:
                try:
                    response = requests.post(
                        self.BERDL_API_URL,
                        headers=self.headers,
                        json=payload,
                        timeout=self.TIMEOUT
                    )

                    if response.status_code == 200:
                        data = response.json()
                        results = data.get('result', [])
                        all_results.extend(results)

                        if not data.get('pagination', {}).get('has_more', False):
                            return pd.DataFrame(all_results)

                        offset += page_limit
                        break

                    elif response.status_code == 408:
                        retries += 1
                        if retries <= self.MAX_RETRIES:
                            time.sleep(self.RETRY_DELAY * (2 ** (retries - 1)))
                        else:
                            raise Exception(f"Query timeout after {self.MAX_RETRIES} retries")

                    elif response.status_code == 403:
                        error_msg = response.json().get('message', '')
                        raise Exception(f"Access denied: {error_msg}")

                    elif response.status_code == 500:
                        error_msg = response.json().get('message', '')
                        raise Exception(f"Server error (500): {error_msg[:100]}")

                    else:
                        raise Exception(f"API error {response.status_code}: {response.text[:200]}")

                except requests.exceptions.Timeout:
                    retries += 1
                    if retries <= self.MAX_RETRIES:
                        time.sleep(self.RETRY_DELAY * (2 ** (retries - 1)))
                    else:
                        raise Exception(f"Request timeout after {self.MAX_RETRIES} retries")

        return pd.DataFrame(all_results)

    def _execute_query_junction(self, query: str, max_results: int = 10000) -> pd.DataFrame:
        """
        Execute a SQL query with smaller page size for junction table.

        Uses JUNCTION_PAGE_LIMIT (100) instead of PAGE_LIMIT (1000) for the
        gene_genecluster_junction table (1B rows) to avoid timeouts.
        """
        all_results = []
        offset = 0

        while len(all_results) < max_results:
            page_limit = min(self.JUNCTION_PAGE_LIMIT, max_results - len(all_results))
            payload = {
                "query": query,
                "limit": page_limit,
                "offset": offset
            }

            retries = 0
            while retries <= self.MAX_RETRIES:
                try:
                    response = requests.post(
                        self.BERDL_API_URL,
                        headers=self.headers,
                        json=payload,
                        timeout=self.TIMEOUT
                    )

                    if response.status_code == 200:
                        data = response.json()
                        results = data.get('result', [])
                        all_results.extend(results)

                        if not data.get('pagination', {}).get('has_more', False):
                            return pd.DataFrame(all_results)

                        offset += page_limit
                        break

                    elif response.status_code == 408:
                        retries += 1
                        if retries <= self.MAX_RETRIES:
                            time.sleep(self.RETRY_DELAY * (2 ** (retries - 1)))
                        else:
                            raise Exception(f"Query timeout after {self.MAX_RETRIES} retries")

                    else:
                        raise Exception(f"API error {response.status_code}: {response.text[:200]}")

                except requests.exceptions.Timeout:
                    retries += 1
                    if retries <= self.MAX_RETRIES:
                        time.sleep(self.RETRY_DELAY * (2 ** (retries - 1)))
                    else:
                        raise Exception(f"Request timeout after {self.MAX_RETRIES} retries")

        return pd.DataFrame(all_results)

    # === QueryPangenomeABC interface implementation ===

    def get_clade_gene_clusters(self, clade_id: str) -> pd.DataFrame:
        """
        Get all gene clusters for a given GTDB species clade.

        Note: Uses LIKE pattern workaround due to '--' character issue.
        """
        like_pattern = self.clade_id_to_like_pattern(clade_id)
        query = f"""
            SELECT * FROM {self.SCHEMA}.gene_cluster
            WHERE gtdb_species_clade_id LIKE '{like_pattern}'
        """
        return self._execute_query(query)

    def get_clade_by_member(self, member_id: str) -> pd.DataFrame:
        """
        Get the clade information for a given genome/member.
        """
        clade_id = self.get_member_representative(member_id)
        return self.get_clade_metadata(clade_id)

    def get_cluster_members(self, cluster_id: str) -> pd.DataFrame:
        """
        Get all genes belonging to a specific gene cluster.

        Note: Uses smaller page size (100) for this 1B row junction table.
              This method may timeout on very large clusters.
        """
        query = f"""
            SELECT * FROM {self.SCHEMA}.gene_genecluster_junction
            WHERE gene_cluster_id = '{cluster_id}'
        """
        return self._execute_query_junction(query)

    def get_clusters_members(self, cluster_ids) -> pd.DataFrame:
        """
        Get all genes belonging to multiple gene clusters.

        Note: Queries one cluster at a time with smaller page size.
              This avoids timeouts on the 1B row junction table.
        """
        if not cluster_ids:
            return pd.DataFrame(columns=['gene_id', 'gene_cluster_id'])

        all_results = []
        for cluster_id in cluster_ids:
            try:
                df = self.get_cluster_members(cluster_id)
                all_results.append(df)
            except Exception as e:
                print(f"Warning: Failed to get members for cluster {cluster_id}: {e}")

        if all_results:
            return pd.concat(all_results, ignore_index=True)
        return pd.DataFrame(columns=['gene_id', 'gene_cluster_id'])

    def get_clade_members(self, clade_id: str) -> pd.DataFrame:
        """
        Get all genomes belonging to a GTDB species clade.

        Note: Uses LIKE pattern workaround due to '--' character issue.
        """
        like_pattern = self.clade_id_to_like_pattern(clade_id)
        query = f"""
            SELECT * FROM {self.SCHEMA}.genome
            WHERE gtdb_species_clade_id LIKE '{like_pattern}'
        """
        return self._execute_query(query)

    def get_member_representative(self, member_id: str) -> str:
        """
        Get the GTDB species clade ID for a given genome.
        """
        query = f"""
            SELECT gtdb_species_clade_id FROM {self.SCHEMA}.genome
            WHERE genome_id = '{member_id}'
        """
        df = self._execute_query(query, limit=10)

        if len(df) == 1:
            return df['gtdb_species_clade_id'].iloc[0]
        elif len(df) == 0:
            raise ValueError(f"Genome not found: {member_id}")
        else:
            raise ValueError(f"Expected exactly 1 result, got {len(df)}")

    def get_member_ani_matrix(self, member_id: str) -> pd.DataFrame:
        """
        Get ANI comparison data for a given genome.
        """
        query = f"""
            SELECT * FROM {self.SCHEMA}.genome_ani
            WHERE genome1_id = '{member_id}' OR genome2_id = '{member_id}'
        """
        return self._execute_query(query)

    def get_clade_metadata(self, clade_id: str) -> pd.DataFrame:
        """
        Get metadata for a GTDB species clade.

        Note: Uses LIKE pattern workaround due to '--' character issue.
        """
        like_pattern = self.clade_id_to_like_pattern(clade_id)
        query = f"""
            SELECT * FROM {self.SCHEMA}.gtdb_species_clade
            WHERE gtdb_species_clade_id LIKE '{like_pattern}'
        """
        return self._execute_query(query)

    # === Additional utility methods ===

    def get_genome_info(self, genome_id: str) -> pd.DataFrame:
        """
        Get information for a specific genome.
        """
        query = f"""
            SELECT * FROM {self.SCHEMA}.genome
            WHERE genome_id = '{genome_id}'
        """
        return self._execute_query(query)

    def get_genomes_by_ids(self, genome_ids: List[str]) -> pd.DataFrame:
        """
        Get information for multiple genomes by their IDs.
        """
        if not genome_ids:
            return pd.DataFrame()

        ids_str = ", ".join([f"'{gid}'" for gid in genome_ids])
        query = f"""
            SELECT * FROM {self.SCHEMA}.genome
            WHERE genome_id IN ({ids_str})
        """
        return self._execute_query(query)
