import polars as pl
from pathlib import Path
from query_genome import QueryGenomeABC


class QueryGenomeLocal(QueryGenomeABC):

    def __init__(self):
        self.root_reference = Path('/data/reference_data/berdl_db/cdm_genomes/ke-pangenomes')
        self.ldf_name = pl.scan_parquet(f'{self.root_reference}/block_*/name.parquet')
        pass

    def get_genome_features(self, genome_id: str):
        """Return all gene clusters for a given clade."""
        raise NotImplementedError

    def get_genome_contigs(self, genome_id: str):
        """Return all gene clusters for a given clade."""
        raise NotImplementedError
