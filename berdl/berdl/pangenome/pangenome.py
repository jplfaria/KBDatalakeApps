import subprocess
import polars as pl
from pathlib import Path
from modelseedpy.core.msgenome import MSFeature, MSGenome
from berdl.pangenome.paths_pangenome import PathsPangenome
from berdl.hash_seq import ProteinSequence
from berdl.fitness import create_genome_fitness_table, map_protein_hash_to_fitness_records


_STRIP_FNA_ = len('/global/cfs/cdirs/kbase/jungbluth/Projects/Project_Pangenome_GTDB/GTDB_v214_download/ftp.ncbi.nlm.nih.gov/')


class BERDLPangenome:

    def __init__(self, query_pg, query_g, paths: PathsPangenome):
        self.pg = query_pg  # pan-genome query api
        self.query_g = query_g
        self.paths = paths.ensure()

    @staticmethod
    def read_cluster_tsv(file_mmseqs_clusters):
        r_to_m = {}
        m_to_r = {}
        with open(str(file_mmseqs_clusters), 'r') as fh:
            line = fh.readline()
            while line:
                r, m = line.strip().split('\t')
                if r not in r_to_m:
                    r_to_m[r] = set()
                r_to_m[r].add(m)
                m_to_r[m] = r
                line = fh.readline()
        return r_to_m, m_to_r

    def extend_gene_to_cluster(self,  m_to_r, d_gene_to_cluster,):
        clade_members = pl.read_csv(self.paths.out_members_tsv, separator='\t')
        data = {
            'genome_id': [],
            'feature_id': [],
            'cluster_id': [],
            'protein_hash': [],
            'mmseqs_rep_hash': [],
        }
        for member_id in {r[0] for r in clade_members.select("genome_id").rows()}:
            genome = MSGenome.from_fasta(str(self.paths.genome_dir / f'{member_id}.faa'))
            for feature in genome.features:
                if feature.seq:
                    protein = ProteinSequence(feature.seq)
                    protein_h = protein.hash_value
                    data['genome_id'].append(member_id)
                    data['feature_id'].append(feature.id)
                    data['cluster_id'].append(d_gene_to_cluster[feature.id])
                    data['protein_hash'].append(protein_h)
                    data['mmseqs_rep_hash'].append(m_to_r[protein_h])
        df = pl.DataFrame(data)
        df.write_parquet(self.paths.root / 'pangenome_cluster_with_mmseqs.parquet')

    def mmseqs2(self, filename_faa: Path, t=30):
        work_dir = self.paths.out_mmseqs_dir
        log_stdout = work_dir / 'log.out'
        log_stderr = work_dir / 'log.err'
        cmd = [
            'mmseqs', 'easy-cluster',
            '--min-seq-id', '0.0',
            '--cov-mode', '0',
            '-c', '0.80',
            '--threads', str(t),
            str(filename_faa.resolve()),
            filename_faa.name[:-4],
            'mmseqs2_tmp',
        ]
        print(f'mmseqs2 run {" ".join(cmd)}')
        with open(log_stdout, 'w') as fh_out, open(log_stderr, 'w') as fh_err:
            process = subprocess.Popen(
                cmd,
                cwd=str(work_dir),
                stdout=fh_out,  # inherit parent stdout
                stderr=fh_err,  # inherit parent stderr
            )
            ret = process.wait()

        if ret != 0:
            print('mmseqs failed exec')
            print('out')
            with open(log_stdout) as fh:
                print(fh.read())
            print('err')
            with open(log_stderr) as fh:
                print(fh.read())
            raise RuntimeError(
                f"mmseqs2 failed with exit code {ret}"
            )

        print("mmseqs2 completed")

    def run(self, selected_clade_member_id):
        clade_id = self.pg.get_member_representative(selected_clade_member_id)
        print('clade_id', clade_id)
        clade_members = self.pg.get_clade_members(clade_id)

        clade_members.write_csv(self.paths.out_members_tsv, separator='\t')

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
            genome_features = self.query_g.get_genome_features(member_id)
            genome_features.to_fasta(self.paths.genome_dir / f'{member_id}.faa')
            members[member_id] = genome_features

            # write contigs
            # FIXME: cant fetch from parquet use temp ref data for contigs
            path_ref_contig = Path('/data/reference_data/contigs')
            file_contigs = path_ref_contig / row['fna_file_path_nersc'][_STRIP_FNA_:]
            if file_contigs.exists():
                genome_contigs = MSGenome.from_fasta(str(file_contigs))
                genome_contigs.to_fasta(self.paths.assembly_dir / f'{member_id}.fna')

        # build master protein user_genome + pangenome

        u_proteins = {}
        for k, g in members.items():
            print(k, len(g.features))
            for feature in g.features:
                _parts = feature.description.split(' ')
                h = _parts[4]
                u_proteins[h] = MSFeature(h, feature.seq)

        print(f"write {clade_id} master faa")
        genome_master_faa = MSGenome()
        genome_master_faa.add_features(list(u_proteins.values()))
        genome_master_faa.to_fasta(str(self.paths.out_master_faa_pangenome_members))

        #  collect user genome and add to u_proteins
        import json
        with open(self.paths.genome_prep_clade_data, 'r') as fh:
            clade_assignment = json.load(fh)

        user_genomes = {k for k, v in clade_assignment.items() if v == selected_clade_member_id}
        for g in user_genomes:
            genome_path = self.paths.genome_prep_genome_folder / f'user_{g}.faa'
            genome_user = MSGenome.from_fasta(str(genome_path.resolve()))
            for feature in genome_user.features:
                protein = ProteinSequence(feature.seq)
                h = protein.hash_value
                u_proteins[h] = MSFeature(h, feature.seq)

        #  collect phenotype and fitness
        genome_master_faa_fitness = MSGenome.from_fasta(str(self.paths.ref_master_faa_protein_fitness))
        for feature in genome_master_faa_fitness.features:
            if feature.id not in u_proteins:
                u_proteins[feature.id] = MSFeature(feature.id, feature.seq)
        genome_master_faa_phenotype = MSGenome.from_fasta(str(self.paths.ref_master_faa_protein_phenotype))
        for feature in genome_master_faa_phenotype.features:
            if feature.id not in u_proteins:
                u_proteins[feature.id] = MSFeature(feature.id, feature.seq)

        # rebuild master faa genome with proteins from
        #  user / pangenome / fitness / phenotypes
        print(f"write {clade_id} user / pangenome / fitness / phenotypes master faa")
        genome_master_faa = MSGenome()
        genome_master_faa.add_features(list(u_proteins.values()))
        genome_master_faa.to_fasta(str(self.paths.out_master_faa))

        # run mmseqs
        self.mmseqs2(self.paths.out_master_faa)

        # read clusters
        print('write extended clusters')
        filename_clusters = self.paths.out_mmseqs_dir / f'{self.paths.out_master_faa.name[:-4]}_cluster.tsv'
        r_to_m, m_to_r = self.read_cluster_tsv(filename_clusters)
        if filename_clusters.exists():
            self.extend_gene_to_cluster(m_to_r, d_gene_to_cluster)

        # map user_genome_to_pangenomes
        print('map user genome to clusters')
        with open(self.paths.genome_prep_clade_data, 'r') as fh:
            user_to_clade = {f'user_{k}.faa' for k, v in json.load(fh).items() if v == selected_clade_member_id}

        # map cluster to fitness
        m_to_fitness_feature = map_protein_hash_to_fitness_records()
        for filename in user_to_clade:
            path_input_genomes = (self.paths.root / '../..' / 'genome').resolve()
            path_genome_faa = path_input_genomes / filename
            input_genome = MSGenome.from_fasta(str(path_genome_faa))
            input_genome_id = path_genome_faa.name[:-4]
            df_input_genome_fitness = create_genome_fitness_table(input_genome, input_genome_id, m_to_r, r_to_m,
                                                                  m_to_fitness_feature)
            df_input_genome_fitness.write_parquet(path_input_genomes / f'user_{input_genome_id}_fitness.parquet')
