import json
import os
import sqlite3

import pandas as pd

from berdl.genome_paths import GenomePaths
from berdl.pangenome.paths_pangenome import PathsPangenome
from modelseedpy import MSGenome
from berdl.hash_seq import ProteinSequence
from pathlib import Path
import polars as pl


def convert_kbase_location(feature_location):
    contig, p0, strand, sz = feature_location
    start = p0
    end = start + sz - 1
    if strand == '-':
        end = p0
        start = end - sz + 1
    return contig, start, end, strand


class DatalakeTableBuilder:

    def __init__(self, root_genome: GenomePaths, root_pangenome: PathsPangenome):
        self.root_genome = root_genome
        self.root_pangenome = root_pangenome
        self.include_dna_sequence = False
        self.include_protein_sequence = False

    def build(self):
        #self.build_genome_table()
        df_user_features = self.build_user_genome_feature_parquet()
        self.build_user_genome_features_table(df_user_features)
        self.build_ani_table()
        #self.build_user_genome_features_table()
        #self.build_pangenome_genome_features_table()
         #self.build_phenotype_kegg_table()
        # model
        # ontology
        #

    def build_genome_table(self):
        conn = sqlite3.connect(str(self.root_pangenome.out_sqlite3_file))
        cur = conn.cursor()
        cur.execute("""
                CREATE TABLE IF NOT EXISTS genome (

                    genome TEXT NOT NULL,
                    gtdb_taxonomy TEXT NOT NULL,
                    ncbi_taxonomy TEXT NOT NULL,
                    n_contigs INT NOT NULL,
                    n_features INT NOT NULL,
                    PRIMARY KEY genome
                );
                """)
        pass

    @staticmethod
    def _update_ontology(feature_ontology_terms, df):
        """
        add df ontologies to dict
        :param feature_ontology_terms: dict to update
        :param df: df with ontologies
        :return:
        """
        for row in df.rows(named=True):
            feature_id = row['feature_id']
            if feature_id not in feature_ontology_terms:
                feature_ontology_terms[feature_id] = {}
            for ontology_term in row.keys():
                values = row[ontology_term]
                if values and ontology_term != 'feature_id':
                    if ontology_term not in feature_ontology_terms[feature_id]:
                        feature_ontology_terms[feature_id][ontology_term] = set()
                    for value in values.split('; '):
                        feature_ontology_terms[feature_id][ontology_term].add(value)

    def collect_annotation(self, genome_id: str, genome_root: Path):
        """
        :param genome_id: genome id
        :param genome_root: root path to scan annotation files related to genome_id
        :return:
        """
        feature_ontology_terms = {}
        for f in os.listdir(str(genome_root)):
            if genome_id in f and '_annotation_' in f and f.endswith('.tsv'):
                print(f'collecting {genome_root / f}')
                df = pl.read_csv(genome_root / f, separator='\t')
                self._update_ontology(feature_ontology_terms, df)

        path_rast = genome_root / f'user_{genome_id}_rast.tsv'
        if path_rast.exists():
            print(f'collecting {path_rast}')
            df = pl.read_csv(path_rast, separator='\t')
            self._update_ontology(feature_ontology_terms, df)

        return feature_ontology_terms

    def build_ani_table(self):
        conn = sqlite3.connect(str(self.root_pangenome.out_sqlite3_file))
        cur = conn.cursor()
        cur.execute("DROP TABLE IF EXISTS ani;")
        cur.execute("""
        CREATE TABLE IF NOT EXISTS ani (
            genome1 TEXT NOT NULL,
            genome2 TEXT NOT NULL,
            ani REAL NOT NULL,
            af1 REAL NOT NULL,
            af2 REAL NOT NULL,
            PRIMARY KEY (genome1, genome2)
        );
        """)
        if self.root_genome.ani_kepangenomes_json.exists():
            with open(self.root_genome.ani_kepangenomes_json) as fh:
                data_ani_clade = json.load(fh)
            for genome_user, ani_matches in data_ani_clade.items():
                for genome_ani, (ani, af1, af2) in ani_matches.items():
                    cur.execute(
                        "INSERT INTO ani (genome1, genome2, ani, af1, af2) VALUES (?, ?, ?, ?, ?)",
                        (f'user_{genome_user}', genome_ani, ani, af1, af2)
                    )

        if self.root_genome.ani_fitness_json.exists():
            with open(self.root_genome.ani_fitness_json) as fh:
                data_ani_fitness = json.load(fh)
            for genome_user, ani_matches in data_ani_fitness.items():
                for genome_ani, (ani, af1, af2) in ani_matches.items():
                    cur.execute(
                        "INSERT INTO ani (genome1, genome2, ani, af1, af2) VALUES (?, ?, ?, ?, ?)",
                        (f'user_{genome_user}', genome_ani, ani, af1, af2)
                    )

        if self.root_genome.ani_phenotypes_json:
            with open(self.root_genome.ani_phenotypes_json) as fh:
                data_ani_phenotypes = json.load(fh)
            for genome_user, ani_matches in data_ani_phenotypes.items():
                for genome_ani, (ani, af1, af2) in ani_matches.items():
                    cur.execute(
                        "INSERT INTO ani (genome1, genome2, ani, af1, af2) VALUES (?, ?, ?, ?, ?)",
                        (f'user_{genome_user}', genome_ani, ani, af1, af2)
                    )

        conn.commit()
        conn.close()

    def build_pangenome_member_feature_parquet(self):
        data = {
            'genome': [],
            'contig': [],
            'feature_id': [],
            #'aliases': [],
            #'length': [],
            'start': [],
            'end': [],
            'strand': [],
            'type': [],
            #'dna_sequence': [],
            'protein_sequence': [],
            'protein_sequence_hash': [],
            'cluster': [],
            'is_core': [],
        }

        genome_id_to_ontologies = {}
        all_ontology_terms = set()
        path_genome_dir = self.root_pangenome.genome_dir
        for f in os.listdir(str(path_genome_dir)):
            if f.endswith('.faa'):
                path_genome_faa = path_genome_dir / f
                genome_id = path_genome_faa.name[:-4]  # strip .faa
                feature_ontology_terms = self.collect_annotation(genome_id, path_genome_dir)
                for o in feature_ontology_terms.values():
                    all_ontology_terms |= set(o)
                genome_id_to_ontologies[genome_id] = feature_ontology_terms

        print(f'found these ontology terms: {all_ontology_terms}')
        for term in all_ontology_terms:
            data[f'ontology_{term}'] = []

        df_clusters = pl.read_parquet(self.root_pangenome / 'pangenome_cluster_with_mmseqs.parquet')
        for f in os.listdir(str(path_genome_dir)):
            if f.endswith('.faa'):
                path_genome_faa = path_genome_dir / f
                genome_id = path_genome_faa.name[:-4]  # strip .faa
                genome = MSGenome.from_fasta(str(path_genome_faa))
                feature_ontology_terms = genome_id_to_ontologies.get(genome_id, {})

                d_feature_cluster = {row['feature_id']: row for row in df_clusters.filter(
                    pl.col("genome_id") == genome_id).rows(named=True)}

                for feature in genome.features:
                    feature_id = feature.id
                    _parts = feature.description.split(' ')
                    contig_id = _parts[0]
                    start = _parts[1]
                    end = _parts[2]
                    strand = _parts[3]
                    feature_type = _parts[4]
                    protein_sequence = feature.seq if self.include_protein_sequence else None
                    protein_hash = ProteinSequence(feature.seq) if feature.seq else None

                    data['genome'].append(genome_id)
                    data['contig'].append(contig_id)
                    data['feature_id'].append(feature_id)
                    #data['aliases'].append(aliases)
                    #data['length'].append(feature_len)
                    data['start'].append(start)
                    data['end'].append(end)
                    data['strand'].append(strand)
                    data['type'].append(feature_type)
                    #data['dna_sequence'].append(dna_sequence)
                    data['protein_sequence'].append(protein_sequence)
                    data['protein_sequence_hash'].append(protein_hash)
                    data['cluster'].append(d_feature_cluster.get(feature_id, {}).get('cluster_id'))
                    data['is_core'].append(d_feature_cluster.get(feature_id, {}).get('is_core'))

                    for term in all_ontology_terms:
                        values = feature_ontology_terms.get(feature_id, {}).get(term)
                        if values is None:
                            data[f'ontology_{term}'].append(None)
                        else:
                            data[f'ontology_{term}'].append('; '.join(values))

        return pl.DataFrame(data)

    def build_user_genome_feature_parquet(self):
        input_genome_dir = Path(self.root_genome.genome_dir)
        data = {
            'genome': [],
            'contig': [],
            'feature_id': [],
            'aliases': [],
            'length': [],
            'start': [],
            'end': [],
            'strand': [],
            'type': [],
            'dna_sequence': [],
            'protein_sequence': [],
            'protein_sequence_hash': [],
            'pangenome_cluster': [],
            'pangenome_is_core': [],
        }

        genome_id_to_ontologies = {}
        all_ontology_terms = set()
        for f in os.listdir(str(input_genome_dir)):
            if f.endswith('.faa'):
                path_genome_faa = input_genome_dir / f
                genome_id = path_genome_faa.name[5:-4]  # strip user_<genome_id>.faa
                path_genome_tsv = input_genome_dir / f'{genome_id}_genome.tsv'
                if path_genome_tsv.exists():
                    print('build_user_genome_features_table - found', genome_id)
                    feature_ontology_terms = self.collect_annotation(genome_id, input_genome_dir)
                    for o in feature_ontology_terms.values():
                        all_ontology_terms |= set(o)
                    genome_id_to_ontologies[genome_id] = feature_ontology_terms

        print(f'found these ontology terms: {all_ontology_terms}')
        for term in all_ontology_terms:
            data[f'ontology_{term}'] = []

        for f in os.listdir(str(input_genome_dir)):
            if f.endswith('.faa'):
                path_genome_faa = input_genome_dir / f
                genome_id = path_genome_faa.name[5:-4]  # strip user_<genome_id>.faa
                path_genome_tsv = input_genome_dir / f'{genome_id}_genome.tsv'
                path_genome_pangenome_profile = input_genome_dir / f'{genome_id}_pangenome_profile.tsv'
                if path_genome_tsv.exists():
                    d_pangenome_profile = {}
                    if path_genome_pangenome_profile.exists():
                        d_pangenome_profile = {row['feature_id']: row for row in
                                               pl.read_csv(path_genome_pangenome_profile, separator='\t').rows(
                                                   named=True)}

                    df_genome_tsv = pl.read_csv(path_genome_tsv, separator='\t')
                    feature_ontology_terms = genome_id_to_ontologies[genome_id]
                    for row in df_genome_tsv.rows(named=True):
                        # contig, p0, strand, sz
                        feature_id = row['gene_id']
                        contig = row['contig']
                        df_start = int(row['start'])
                        df_end = row['end']
                        strand = row['strand']
                        start = df_start
                        end = df_end
                        if strand == '-':
                            start = df_end
                            end = df_start
                        aliases = row['aliases']
                        feature_len = end - start
                        feature_type = row['type']
                        dna_sequence = row['dna_sequence'] if self.include_dna_sequence else None
                        protein_sequence = row['protein_translation'] if self.include_protein_sequence else None
                        protein_hash = ProteinSequence(row['protein_translation']).hash_value \
                            if row.get("protein_translation", "") else None

                        data['genome'].append(genome_id)
                        data['contig'].append(contig)
                        data['feature_id'].append(feature_id)
                        data['aliases'].append(aliases)
                        data['length'].append(feature_len)
                        data['start'].append(start)
                        data['end'].append(end)
                        data['strand'].append(strand)
                        data['type'].append(feature_type)
                        data['dna_sequence'].append(dna_sequence)
                        data['protein_sequence'].append(protein_sequence)
                        data['protein_sequence_hash'].append(protein_hash)
                        data['pangenome_cluster'].append(
                            d_pangenome_profile.get('feature_id', {}).get('pangenome_cluster'))
                        data['pangenome_is_core'].append(
                            d_pangenome_profile.get('feature_id', {}).get('is_core')
                        )
                        for term in all_ontology_terms:
                            values = feature_ontology_terms.get(feature_id, {}).get(term)
                            if values is None:
                                data[f'ontology_{term}'].append(None)
                            else:
                                data[f'ontology_{term}'].append('; '.join(values))
        return pl.DataFrame(data)

    def build_user_genome_features_table(self, df_user_features):
        conn = sqlite3.connect(str(self.root_pangenome.out_sqlite3_file))
        cur = conn.cursor()
        cur.execute("DROP TABLE IF EXISTS user_feature;")
        sql_create_table = """
        CREATE TABLE IF NOT EXISTS user_feature (
                    genome TEXT NOT NULL,
                    contig TEXT NOT NULL,
                    feature_id TEXT NOT NULL,
                    aliases TEXT,
                    length INTEGER NOT NULL,
                    start INTEGER NOT NULL,
                    end INTEGER NOT NULL,
                    strand TEXT NOT NULL,
                    type TEXT NOT NULL,
                    dna_sequence TEXT,
                    protein_sequence TEXT,
                    protein_sequence_hash TEXT,
                    pangenome_cluster TEXT,
                    pangenome_is_core INTEGER,
        """
        ontology_cols = {col for col in df_user_features.columns if col.startswith('ontology_')}
        for col in ontology_cols:
            sql_create_table += f"{col} TEXT,\n"
        sql_create_table += "PRIMARY KEY (genome, contig, feature_id));"
        cur.execute(sql_create_table)

        df_user_features.to_pandas().to_sql("user_feature", conn, if_exists="append", index=False)

        conn.commit()
        conn.close()

    def build_pangenome_genome_features_table(self):
        conn = sqlite3.connect(str(self.root_pangenome.out_sqlite3_file))
        cur = conn.cursor()
        cur.execute("""
                CREATE TABLE IF NOT EXISTS pangenome_feature (
                    genome TEXT NOT NULL,
                    contig TEXT NOT NULL,
                    feature_id TEXT NOT NULL,
                    length INTEGER NOT NULL,
                    start INTEGER NOT NULL,
                    end INTEGER NOT NULL,
                    strand TEXT NOT NULL,
                    type TEXT NOT NULL,
                    sequence TEXT NOT NULL,
                    sequence_hash TEXT NOT NULL,
                    pangenome_cluster TEXT NOT NULL,
                    pangenome_is_core INTEGER NOT NULL,
                    PRIMARY KEY (genome, contig, feature_id)
                );
                """)

        # read genome member list
        # genome = MSGenome.from_fasta(fffff)
        pass

    def build_phenotype_kegg_table(self):
        pass
