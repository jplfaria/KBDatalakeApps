

class BERDLPangenome:

    def __init__(self):
        self.pg = None  # pan-genome query api
        pass

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

        # build master protein user_genome + pangenome

        pass
