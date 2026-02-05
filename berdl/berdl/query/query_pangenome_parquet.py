"""
Local Parquet implementation of QueryPangenomeABC.

This module provides pangenome queries against local Parquet files.
It implements the same interface as QueryPangenomeBERDL but retrieves data
from local files instead of the remote API.

Author: Jose P. Faria (jplfaria@gmail.com)
Date: February 2026

=== DATA FILES ===

Expected files at base_path (default: /data/reference_data/):
- genome.parquet: Genome metadata (293K rows)
- gene_cluster.parquet: Gene clusters (132M rows)
- genome_ani.parquet: ANI comparisons (421M rows)
- gene_genecluster_junction.parquet: Gene-cluster mappings (1B rows)
- gtdb_species_clade.parquet: GTDB species clade info (27K rows)

=== USAGE ===

    from berdl.query.query_pangenome_parquet import QueryPangenomeParquet

    query = QueryPangenomeParquet(base_path="/data/reference_data")

    # Get gene clusters for a clade
    clusters = query.get_clade_gene_clusters("s__Escherichia_coli")
"""

import pandas as pd
import pyarrow.parquet as pq
from pathlib import Path
from typing import List, Optional, Iterable

from .query_pangenome import QueryPangenomeABC


class QueryPangenomeParquet(QueryPangenomeABC):
    """
    Query KBase pangenome data from local Parquet files.

    Tables:
        - genome.parquet (293K rows)
        - gene_cluster.parquet (132M rows)
        - genome_ani.parquet (421M rows)
        - gene_genecluster_junction.parquet (1B rows)
        - gtdb_species_clade.parquet (27K rows)
    """

    DEFAULT_BASE_PATH = "/data/reference_data"

    # Cache for loaded DataFrames
    _genome_df: pd.DataFrame = None
    _gene_cluster_df: pd.DataFrame = None
    _genome_ani_df: pd.DataFrame = None
    _gene_genecluster_df: pd.DataFrame = None
    _gtdb_species_clade_df: pd.DataFrame = None

    def __init__(self, base_path: str = None):
        """
        Initialize the local parquet pangenome query client.

        Args:
            base_path: Path to directory containing parquet files.
                       Default: /data/reference_data
        """
        raise NotImplementedError("AI hallucinated. Files do not exists")
        self.base_path = Path(base_path or self.DEFAULT_BASE_PATH)

        # File paths
        self.genome_path = self.base_path / "genome.parquet"
        self.gene_cluster_path = self.base_path / "gene_cluster.parquet"
        self.genome_ani_path = self.base_path / "genome_ani.parquet"
        self.gene_genecluster_path = self.base_path / "gene_genecluster_junction.parquet"
        self.gtdb_species_clade_path = self.base_path / "gtdb_species_clade.parquet"

    def _load_genome(self) -> pd.DataFrame:
        """Load genome parquet file (cached)."""
        if QueryPangenomeParquet._genome_df is None:
            if not self.genome_path.exists():
                raise FileNotFoundError(f"Genome file not found: {self.genome_path}")
            print(f"    Loading genome data from {self.genome_path}...")
            QueryPangenomeParquet._genome_df = pq.read_table(self.genome_path).to_pandas()
            print(f"    Loaded {len(QueryPangenomeParquet._genome_df)} genomes")
        return QueryPangenomeParquet._genome_df

    def _load_gene_cluster(self) -> pd.DataFrame:
        """Load gene_cluster parquet file (cached)."""
        if QueryPangenomeParquet._gene_cluster_df is None:
            if not self.gene_cluster_path.exists():
                raise FileNotFoundError(f"Gene cluster file not found: {self.gene_cluster_path}")
            print(f"    Loading gene cluster data from {self.gene_cluster_path}...")
            QueryPangenomeParquet._gene_cluster_df = pq.read_table(self.gene_cluster_path).to_pandas()
            print(f"    Loaded {len(QueryPangenomeParquet._gene_cluster_df)} gene clusters")
        return QueryPangenomeParquet._gene_cluster_df

    def _load_genome_ani(self) -> pd.DataFrame:
        """Load genome_ani parquet file (cached)."""
        if QueryPangenomeParquet._genome_ani_df is None:
            if not self.genome_ani_path.exists():
                raise FileNotFoundError(f"Genome ANI file not found: {self.genome_ani_path}")
            print(f"    Loading genome ANI data from {self.genome_ani_path}...")
            QueryPangenomeParquet._genome_ani_df = pq.read_table(self.genome_ani_path).to_pandas()
            print(f"    Loaded {len(QueryPangenomeParquet._genome_ani_df)} ANI records")
        return QueryPangenomeParquet._genome_ani_df

    def _load_gene_genecluster(self) -> pd.DataFrame:
        """Load gene_genecluster_junction parquet file (cached)."""
        if QueryPangenomeParquet._gene_genecluster_df is None:
            if not self.gene_genecluster_path.exists():
                raise FileNotFoundError(f"Gene-cluster junction file not found: {self.gene_genecluster_path}")
            print(f"    Loading gene-cluster junction data from {self.gene_genecluster_path}...")
            QueryPangenomeParquet._gene_genecluster_df = pq.read_table(self.gene_genecluster_path).to_pandas()
            print(f"    Loaded {len(QueryPangenomeParquet._gene_genecluster_df)} gene-cluster mappings")
        return QueryPangenomeParquet._gene_genecluster_df

    def _load_gtdb_species_clade(self) -> pd.DataFrame:
        """Load gtdb_species_clade parquet file (cached)."""
        if QueryPangenomeParquet._gtdb_species_clade_df is None:
            if not self.gtdb_species_clade_path.exists():
                raise FileNotFoundError(f"GTDB species clade file not found: {self.gtdb_species_clade_path}")
            print(f"    Loading GTDB species clade data from {self.gtdb_species_clade_path}...")
            QueryPangenomeParquet._gtdb_species_clade_df = pq.read_table(self.gtdb_species_clade_path).to_pandas()
            print(f"    Loaded {len(QueryPangenomeParquet._gtdb_species_clade_df)} GTDB species clades")
        return QueryPangenomeParquet._gtdb_species_clade_df

    # === QueryPangenomeABC interface implementation ===

    def get_clade_gene_clusters(self, clade_id: str) -> pd.DataFrame:
        """
        Get all gene clusters for a given GTDB species clade.
        """
        df = self._load_gene_cluster()
        return df[df['gtdb_species_clade_id'] == clade_id].copy()

    def get_clade_by_member(self, member_id: str) -> pd.DataFrame:
        """
        Get the clade information for a given genome/member.
        """
        clade_id = self.get_member_representative(member_id)
        return self.get_clade_metadata(clade_id)

    def get_cluster_members(self, cluster_id: str) -> pd.DataFrame:
        """
        Get all genes belonging to a specific gene cluster.
        """
        df = self._load_gene_genecluster()
        return df[df['gene_cluster_id'] == cluster_id].copy()

    def get_clusters_members(self, cluster_ids: Iterable[str]) -> pd.DataFrame:
        """
        Get all genes belonging to multiple gene clusters.
        """
        if not cluster_ids:
            return pd.DataFrame(columns=['gene_id', 'gene_cluster_id'])

        cluster_ids_set = set(cluster_ids)
        df = self._load_gene_genecluster()
        return df[df['gene_cluster_id'].isin(cluster_ids_set)].copy()

    def get_clade_members(self, clade_id: str) -> pd.DataFrame:
        """
        Get all genomes belonging to a GTDB species clade.
        """
        df = self._load_genome()
        return df[df['gtdb_species_clade_id'] == clade_id].copy()

    def get_member_representative(self, member_id: str) -> str:
        """
        Get the GTDB species clade ID for a given genome.
        """
        df = self._load_genome()
        result = df[df['genome_id'] == member_id]

        if len(result) == 1:
            return result['gtdb_species_clade_id'].iloc[0]
        elif len(result) == 0:
            raise ValueError(f"Genome not found: {member_id}")
        else:
            raise ValueError(f"Expected exactly 1 result, got {len(result)}")

    def get_member_ani_matrix(self, member_id: str) -> pd.DataFrame:
        """
        Get ANI comparison data for a given genome.
        """
        df = self._load_genome_ani()
        mask = (df['genome1_id'] == member_id) | (df['genome2_id'] == member_id)
        return df[mask].copy()

    def get_clade_metadata(self, clade_id: str) -> pd.DataFrame:
        """
        Get metadata for a GTDB species clade.
        """
        df = self._load_gtdb_species_clade()
        return df[df['gtdb_species_clade_id'] == clade_id].copy()

    # === Additional utility methods ===

    def get_genome_info(self, genome_id: str) -> pd.DataFrame:
        """
        Get information for a specific genome.
        """
        df = self._load_genome()
        return df[df['genome_id'] == genome_id].copy()

    def get_genomes_by_ids(self, genome_ids: List[str]) -> pd.DataFrame:
        """
        Get information for multiple genomes by their IDs.
        """
        if not genome_ids:
            return pd.DataFrame()

        df = self._load_genome()
        return df[df['genome_id'].isin(genome_ids)].copy()

    @classmethod
    def clear_cache(cls):
        """Clear all cached data to free memory."""
        cls._genome_df = None
        cls._gene_cluster_df = None
        cls._genome_ani_df = None
        cls._gene_genecluster_df = None
        cls._gtdb_species_clade_df = None
