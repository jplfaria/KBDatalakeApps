from query_pangenome import QueryPangenomeABC
import os
from pathlib import Path
import polars as pl
from modelseedpy import MSGenome


_STRIP_FNA_ = len('/global/cfs/cdirs/kbase/jungbluth/Projects/Project_Pangenome_GTDB/GTDB_v214_download/ftp.ncbi.nlm.nih.gov/')
_STRIP_FAA_ = len('/global/cfs/cdirs/kbase/jungbluth/Projects/Project_Pangenome_GTDB/')


class QueryPangenomeLocal(QueryPangenomeABC):

    def __init__(self, base_dir='/home/fliu/scratch/data/biodb/ke-pangenomes'):
        self.base_dir = base_dir
        parquet_dir = Path(f"{self.base_dir}/parquet/table_gene_cluster_V1.0")
        self.lf_gene_cluster = pl.scan_parquet(str(parquet_dir / "*.parquet"))
        parquet_dir = Path(f"{self.base_dir}/parquet/table_genome_ANI_V1.1")
        self.lf_genome_ani = pl.scan_parquet(str(parquet_dir / "*.parquet"))
        self.lf_gtdb_species_clade = pl.scan_csv(f'{self.base_dir}/csv/table_gtdb_species_clade_V1.1.tsv',
                                                 separator='\t')
        self.lf_genome = pl.scan_csv(f'{self.base_dir}/csv/table_genome_V1.1.tsv', separator='\t')
        self.lf_gene_genecluster = pl.scan_csv(f'{self.base_dir}/csv/table_gene_genecluster_junction_V1.0.tsv',
                                               separator='\t')
        self.lf_gtdb_metadata = pl.scan_csv(f'{self.base_dir}/csv/table_gtdb_species_clade_V1.1.tsv', separator='\t')
        self.lf_ncbi_env = pl.scan_csv(f'{self.base_dir}/csv/table_NCBI_env_V1.1.tsv', separator='\t')

    def get_clade_gene_clusters(self, clade_id):
        return self.lf_gene_cluster.filter((pl.col("gtdb_species_clade_id") == clade_id)).collect()

    def get_clade_by_member(self, member_id):
        pass

    def get_cluster_members(self, cluster_id: str):
        return self.lf_gene_genecluster.filter((pl.col("gene_cluster_id") == cluster_id)).collect()

    def get_clusters_members(self, cluster_ids):
        return self.lf_gene_genecluster.filter((pl.col("gene_cluster_id").is_in(cluster_ids))).collect()

    def get_clade_members(self, clade_id):
        return self.lf_genome.filter(
            pl.col("gtdb_species_clade_id") == clade_id
        ).collect()

    def get_member_representative(self, member_id):
        res = self.lf_genome.filter(pl.col("genome_id") == member_id).select("gtdb_species_clade_id").collect()
        if res.shape == (1, 1):
            return res.to_dict()['gtdb_species_clade_id'][0]
        else:
            raise ValueError(f'expected exactly 1 result. Got - {res.shape}')

    def get_member_ani_matrix(self, member_id):
        return self.lf_genome_ani.filter(
            (pl.col("genome1_id") == member_id) |
            (pl.col("genome2_id") == member_id)
        ).collect()

    def get_clade_metadata(self, clade_id):

        pass

    def get_member_genome(self, member_id):
        faa_file_path_nersc = self.lf_genome.filter(
            pl.col("genome_id") == member_id
        ).select("faa_file_path_nersc").collect().item()
        _path = faa_file_path_nersc[_STRIP_FAA_:]
        f = f'/home/fliu/scratch/data/biodb/ke-pangenomes/{_path}'
        if os.path.exists(f):
            return MSGenome.from_fasta(f)
        return None

    def get_member_assembly(self, member_id):
        fna_file_path_nersc = self.lf_genome.filter(
            pl.col("genome_id") == member_id
        ).select("fna_file_path_nersc").collect().item()
        _path = fna_file_path_nersc[_STRIP_FNA_:]
        f = f'/home/fliu/scratch/data/biodb/ke-pangenomes/{_path}'
        if os.path.exists(f):
            return MSGenome.from_fasta(f)
        return None
