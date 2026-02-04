import polars as pl
from pathlib import Path
from berdl.query.query_genome import QueryGenomeABC


class QueryGenomeLocal(QueryGenomeABC):

    def __init__(self):
        self.root_reference = Path('/data/reference_data/berdl_db/cdm_genomes/ke-pangenomes')
        self.ldf_name = pl.scan_parquet(f'{self.root_reference}/block_*/*name*.parquet')
        #self.ldf_name_feature = pl.scan_parquet(f'{self.root_reference}/block_*/name_feature.parquet')
        self.ldf_ccol_x_protein = pl.scan_parquet(f'{self.root_reference}/block_*/contig_collection_x_protein.parquet')
        self.ldf_ccol_x_feature = pl.scan_parquet(f'{self.root_reference}/block_*/contig_collection_x_feature.parquet')
        self.ldf_contig_x_feature = pl.scan_parquet(f'{self.root_reference}/block_*/contig_x_feature.parquet')
        self.ldf_protein = pl.scan_parquet(f'{self.root_reference}/block_*/protein.parquet')
        self.ldf_feature = pl.scan_parquet(f'{self.root_reference}/block_*/feature.parquet')
        self.ldf_feature_x_protein = pl.scan_parquet(f'{self.root_reference}/block_*/feature_x_protein.parquet')

    def get_cdm_genome_id_by_name(self, name: str):
        df = self.ldf_name.filter(pl.col("name") == name).collect()
        if df.height == 1:
            return df.row(0, named=True)
        return None

    def get_genome_features(self, genome_id: str):
        res = self.get_cdm_genome_id_by_name(genome_id)
        if res:
            ccol_id = res['entity_id']
            df = self.ldf_ccol_x_feature.filter(pl.col("contig_collection_id") == ccol_id).select(
                pl.col("feature_id")).collect()
            cdm_feat_ids = {o[0] for o in df.rows()}

            df_features = self.ldf_feature.filter(pl.col("feature_id").is_in(cdm_feat_ids)).select(
                ["feature_id", "start", "end", "type"]).collect()

            df_feature_to_contig = self.ldf_contig_x_feature.filter(
                pl.col("feature_id").is_in(cdm_feat_ids)).collect()
            d_feature_to_contig = {row['feature_id']: row['contig_id'] for row in df_feature_to_contig.rows(named=True)}

            df_feature_to_protein = self.ldf_feature_x_protein.filter(
                pl.col("feature_id").is_in(cdm_feat_ids)).collect()
            d_feature_to_protein = {row['feature_id']: row['protein_id'] for row in
                                    df_feature_to_protein.rows(named=True)}

            cdm_prot_ids = set(d_feature_to_protein.values())
            df_proteins = self.ldf_protein.filter(pl.col("protein_id").is_in(cdm_prot_ids)).collect()
            cdm_proteins = {}
            for row in df_proteins.rows(named=True):
                cdm_proteins[row['protein_id']] = row

            cntg_ids = set(d_feature_to_contig.values())
            df_contig_names = self.ldf_name.filter(pl.col("entity_id").is_in(cntg_ids)).select(
                ["entity_id", "name"]).collect()
            d_cntg_id_to_name = {row['entity_id']: row['name'] for row in df_contig_names.rows(named=True)}

            df_names = self.ldf_name.filter(pl.col("entity_id").is_in(cdm_feat_ids)).select(
                ["entity_id", "name"]).collect()

            d_names = {row['entity_id']: row['name'] for row in df_names.rows(named=True)}

            from modelseedpy.core.msgenome import MSFeature, MSGenome
            features = {}
            for row in df_features.rows(named=True):
                start = row['start']
                end = row['end']
                feature_type = row['type']
                cdm_feature_id = row['feature_id']
                cdm_protein_id = d_feature_to_protein[cdm_feature_id]
                cdm_protein = cdm_proteins[cdm_protein_id]
                protein_hash = cdm_protein['hash']
                sequence = cdm_protein['sequence']
                cdm_contig_id = d_feature_to_contig[cdm_feature_id]
                contig_name = d_cntg_id_to_name[cdm_contig_id]
                feature_to_kpg_id = d_names.get(cdm_feature_id, cdm_feature_id)
                desc = f'{contig_name} {start} {end} {feature_type} {protein_hash} {cdm_contig_id} {cdm_feature_id} {cdm_protein_id}'
                genome_feature = MSFeature(feature_to_kpg_id, sequence, description=desc)
                features[genome_feature.id] = genome_feature

            genome = MSGenome()
            genome.add_features(list(features.values()))
            return genome

        else:
            return None

    def get_genome_contigs(self, genome_id: str):
        """Return all gene clusters for a given clade."""
        raise NotImplementedError
