"""
Local Parquet implementation of QueryPangenomeABC.

=============================================================================
DISABLED: Pangenome parquet files are too large to load into memory
=============================================================================

The pangenome parquet files are 100+ GB in size:
- gene_cluster.parquet: 132M rows
- genome_ani.parquet: 421M rows
- gene_genecluster_junction.parquet: 1B rows

These cannot be loaded into pandas DataFrames. Philippe has a separate
implementation that uses streaming/lazy evaluation for these queries.

This module is kept for reference but should NOT be used.
=============================================================================

Author: Jose P. Faria (jplfaria@gmail.com)
Date: February 2026
"""

# =============================================================================
# DISABLED - DO NOT USE - Files too large (100+ GB) to load into memory
# =============================================================================

# Stub class to prevent import errors - raises error if actually used
class QueryPangenomeParquet:
    """
    DISABLED: Pangenome parquet files are too large (100+ GB) to load into memory.

    Use Philippe's streaming implementation instead.
    This stub class exists only to prevent import errors.
    """

    def __init__(self, *args, **kwargs):
        raise NotImplementedError(
            "QueryPangenomeParquet is DISABLED. "
            "Pangenome parquet files are too large (100+ GB) to load into memory. "
            "Use Philippe's streaming implementation instead."
        )
