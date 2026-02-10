import json
import os
from modelseedpy import MSGenome
from berdl.genome_paths import GenomePaths
import pyarrow as pa
import pandas as pd


def run_ani(query, library, output_file, threads=20):
    import subprocess
    cmd = [
        'skani', 'search', "-t", str(threads),
        "--ql", str(query),
        "-d", str(library),
        "-o", str(output_file),
    ]
    print(' '.join(cmd))
    output = subprocess.run(
        cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True
    )
    return output


def _read_search_output_as_parquet(filename, sep='\t'):
    float_types = {'ANI', 'Align_fraction_ref', 'Align_fraction_query'}

    with open(str(filename), 'r') as fh:
        header = fh.readline()

        names = header.strip().split(sep)
        col_number = len(names)
        data = {n: [] for n in names}

        line = fh.readline()
        while line:
            row = line.strip().split(sep)
            for i in range(col_number):
                _d = row[i]
                if names[i] in float_types:
                    _d = float(_d)
                data[names[i]].append(_d)
            line = fh.readline()

        if len(data[names[0]]) == 0:
            return None

        pa_table = pa.Table.from_arrays([pa.array(data[n]) for n in names], names=names)
        return pa_table


class BERDLPreGenome:

    def __init__(self, kbase, paths: GenomePaths):
        self.kbase = kbase
        self.paths = paths.ensure()

    def get_genome_assembly(self, g):
        assembly = self.kbase.get_from_ws(g.assembly_ref)
        handle_local_path = f'/tmp/{assembly.fasta_handle_ref}'
        if not os.path.exists(handle_local_path):
            print('fetch from kbase')
            self.kbase.download_file_from_kbase2(assembly.fasta_handle_ref, handle_local_path)
        return MSGenome.from_fasta(handle_local_path)

    def pre_user_genomes(self, genomes):
        user_genome_files = {}
        with open(str(self.paths.library_user_genomes), 'w') as fh:
            for g in genomes:
                genome_id = g.info.id
                file_genome = self.paths.genome_dir / f'user_{genome_id}.faa'
                file_assembly = self.paths.assembly_dir / f'user_{genome_id}.fna'
                if not file_genome.exists():
                    g.to_fasta(str(file_genome))
                if not file_assembly.exists():
                    g_assembly = self.get_genome_assembly(g)
                    g_assembly.to_fasta(str(file_assembly))

                fh.write(f'{file_assembly.resolve()}\n')
                user_genome_files[genome_id] = (file_assembly.name, file_genome.name)
        return user_genome_files

    def run_ani_databases(self):
        if not self.paths.ani_kepangenomes_out.exists():
            program_out = run_ani(self.paths.library_user_genomes,
                                  self.paths.skani_fast_kepangenome_db,
                                  self.paths.ani_kepangenomes_out)
        if not self.paths.ani_fitness_out.exists():
            program_out = run_ani(self.paths.library_user_genomes,
                                  self.paths.skani_fast_fitness_db,
                                  self.paths.ani_fitness_out)
        if not self.paths.ani_phenotypes_out.exists():
            program_out = run_ani(self.paths.library_user_genomes,
                                  self.paths.skani_fast_phenotypes_db,
                                  self.paths.ani_phenotypes_out)
        df_ani_clade = _read_search_output_as_parquet(self.paths.ani_kepangenomes_out).to_pandas()
        df_ani_fitness = _read_search_output_as_parquet(self.paths.ani_fitness_out).to_pandas()
        df_ani_phenotype = _read_search_output_as_parquet(self.paths.ani_phenotypes_out).to_pandas()
        return df_ani_clade, df_ani_fitness, df_ani_phenotype

    @staticmethod
    def t_ncbi_to_gtdb_id(s):
        a = s.split('_')
        if s.startswith('GCA'):
            return f'GB_{a[0]}_{a[1]}'
        elif s.startswith('GCF'):
            return f'RS_{a[0]}_{a[1]}'

    @staticmethod
    def ani_transform(df, fn_q_transform, fn_r_transform):
        query_to_ref = {}
        for row_id, d in df.iterrows():
            q = fn_q_transform(d['Query_file'])
            r = fn_r_transform(d['Ref_file'])
            if q not in query_to_ref:
                query_to_ref[q] = {}
            query_to_ref[q][r] = [d['ANI'], d['Align_fraction_ref'], d['Align_fraction_query']]
        return query_to_ref

    def ani_translate_clade(self, df, assembly_to_user_id):
        def q_transform(s):
            return assembly_to_user_id[s.split('/')[-1]]

        def r_transform(s):
            return self.t_ncbi_to_gtdb_id(s.split('/')[-1].rsplit('_', 1)[0])

        t = self.ani_transform(df, q_transform, r_transform)
        return t

    def ani_translate_fitness(self, df, assembly_to_user_id):
        def q_transform(s):
            return assembly_to_user_id[s.split('/')[-1]]
        t = self.ani_transform(df, q_transform, lambda x: x.split('/')[-1])
        return t

    def ani_translate_phenotype(self, df, assembly_to_user_id):
        def q_transform(s):
            return assembly_to_user_id[s.split('/')[-1]]

        df_anl_ecoli = pd.read_csv('/data/reference_data/meta_phenotype_ecoli.tsv', sep='\t')
        df_pmi = pd.read_csv('/data/reference_data/meta_phenotype_pmi.tsv', sep='\t')
        df_leaf = pd.read_csv('/data/reference_data/meta_phenotype_leaf.tsv', sep='\t')
        contig_h_to_genome_id = {}
        for row_id, row in df_anl_ecoli.iterrows():
            if row['h'] not in contig_h_to_genome_id:
                contig_h_to_genome_id[row['h']] = row['genome_id']
            else:
                raise ValueError('dup')
        for row_id, row in df_pmi.iterrows():
            h = row['hash_contigset']
            if h not in contig_h_to_genome_id:
                contig_h_to_genome_id[h] = row['genome_id']
            else:
                raise ValueError('dup')
        for row_id, row in df_leaf.iterrows():
            h = row['hash_contigset']
            if h not in contig_h_to_genome_id:
                contig_h_to_genome_id[h] = row['genome_id']
            else:
                raise ValueError('dup')

        def r_transform_phenotypes(s):
            return contig_h_to_genome_id[s.split('/')[-1][:-4]]

        t = self.ani_transform(df, q_transform, r_transform_phenotypes)
        return t

    def run(self, genomes: list):
        user_genome_files = self.pre_user_genomes(genomes)
        df_ani_clade, df_ani_fitness, df_ani_phenotype = self.run_ani_databases()

        assembly_to_user_id = {v[0]: k for k, v in user_genome_files.items()}

        def match_top_clade(ani_clades):
            top_matches = {}

            for clade, genomes in ani_clades.items():
                # Select the genome with the highest ANI (first value in the list)
                best_genome = max(genomes.items(), key=lambda item: item[1][0])[0]
                top_matches[clade] = best_genome

            return top_matches

        ani_clades = self.ani_translate_clade(df_ani_clade, assembly_to_user_id)
        user_to_clade = match_top_clade(ani_clades)

        with open(self.paths.json_user_to_clade, 'w') as fh:
            fh.write(json.dumps(user_to_clade))

        ani_fitness = self.ani_translate_fitness(df_ani_fitness, assembly_to_user_id)
        ani_phenotype = self.ani_translate_phenotype(df_ani_phenotype, assembly_to_user_id)

        with open(self.paths.ani_kepangenomes_json, 'w') as fh:
            fh.write(json.dumps(ani_clades))
        with open(self.paths.ani_fitness_json, 'w') as fh:
            fh.write(json.dumps(ani_fitness))
        with open(self.paths.ani_phenotypes_json, 'w') as fh:
            fh.write(json.dumps(ani_phenotype))

        return user_genome_files, user_to_clade, ani_clades, ani_fitness, ani_phenotype
