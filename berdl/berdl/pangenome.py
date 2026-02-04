from modelseedpy.core.msgenome import MSFeature, MSGenome
from berdl.pangenome.paths_pangenome import PathsPangenome


class BERDLPangenome:

    def __init__(self, query_pg, query_g, paths: PathsPangenome):
        self.pg = query_pg  # pan-genome query api
        self.query_g = query_g
        self.paths = paths.ensure()

    def run(self, genome_id):
        clade_id = self.pg.get_member_representative(genome_id)
        clade_members = self.pg.get_clade_members(clade_id)
        clade_gene_clusters = self.pg.get_clade_gene_clusters(clade_id)
        clade_cluster_ids = set(clade_gene_clusters['gene_cluster_id'])
        df_gene_genecluster = self.pg.get_clusters_members(clade_cluster_ids)
        d_gene_to_cluster = {o[0]: o[1] for o in df_gene_genecluster.iter_rows()}
        d_cluster_to_genes = {}
        for k, v in d_gene_to_cluster.items():
            if v not in d_cluster_to_genes:
                d_cluster_to_genes[v] = set()
            d_cluster_to_genes[v].add(k)

        # get clade member faa and fna
        members = {}
        for row in clade_members.rows(named=True):
            member_id = row['genome_id']
            members[member_id] = self.query_g.get_genome_features(member_id)

        # build master protein user_genome + pangenome

        u_proteins = {}
        for k, g in members.items():
            print(k, len(g.features))
            for feature in g.features:
                _parts = feature.description.split(' ')
                h = _parts[4]
                u_proteins[h] = MSFeature(h, feature.seq)

        genome_master_faa = MSGenome()
        genome_master_faa.add_features(list(u_proteins.values()))
        genome_master_faa.to_fasta(str(self.paths.out_master_faa_pangenome_members))

        #  collect user genome and add to u_proteins
        

        #  collect phenotype and fitness
        genome_master_faa_fitness = MSGenome.from_fasta(str(self.paths.ref_master_faa_protein_fitness))
        for feature in genome_master_faa_fitness:
            if feature.id not in u_proteins:
                u_proteins[feature.id] = MSFeature(feature.id, feature.seq)
        genome_master_faa_phenotype = MSGenome.from_fasta(str(self.paths.ref_master_faa_protein_phenotype))
        for feature in genome_master_faa_phenotype:
            if feature.id not in u_proteins:
                u_proteins[feature.id] = MSFeature(feature.id, feature.seq)

        # rebuild master faa genome with proteins from
        #  user / pangenome / fitness / phenotypes
        genome_master_faa = MSGenome()
        genome_master_faa.add_features(list(u_proteins.values()))
        genome_master_faa.to_fasta(str(self.paths.out_master_faa))
