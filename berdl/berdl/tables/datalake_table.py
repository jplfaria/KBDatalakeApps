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
from berdl.query.query_ontology_local import QueryOntologyLocal


def convert_kbase_location(feature_location):
    contig, p0, strand, sz = feature_location
    start = p0
    end = start + sz - 1
    if strand == '-':
        end = p0
        start = end - sz + 1
    return contig, start, end, strand


class DatalakeTableBuilder:

    def __init__(self, root_genome: GenomePaths, root_pangenome: PathsPangenome,
                 input_genomes: list,
                 include_dna_sequence=True, include_protein_sequence=True,
                 export_tables=False):
        self.input_genomes = set(input_genomes)
        self.root_genome = root_genome
        self.root_pangenome = root_pangenome
        self.include_dna_sequence = include_dna_sequence
        self.include_protein_sequence = include_protein_sequence

        path_pangenome_members = self.root_pangenome.root / 'members.tsv'
        self.df_members = pl.read_csv(path_pangenome_members, separator='\t')
        self.filter_genome_ids = {o[0] for o in self.df_members.select("genome_id").rows()}
        self.filter_genome_ids |= self.input_genomes
        self.export_tables = export_tables

    def build(self):
        self.build_genome_table()
        self.build_ani_table()
        df_user_features = self.build_user_genome_feature_parquet()
        self.build_user_genome_features_table(df_user_features)

        path_genome_dir = self.root_pangenome.genome_dir
        df_pangenome_features = None
        if path_genome_dir.exists():
            df_pangenome_features = self.build_pangenome_member_feature_parquet()
            self.build_pangenome_genome_features_table(df_pangenome_features)

        self.build_input_genome_reactions()
        self.build_gene_reaction_data()
        self.build_genome_phenotype()
        self.build_gene_phenotype()

        if self.export_tables:
            path_table_root = self.root_pangenome.root / 'table'
            path_table_root.mkdir(exist_ok=True)
            df_user_features.write_parquet(path_table_root / 'user_feature.parquet')
            if df_pangenome_features is not None:
                df_pangenome_features.write_parquet(path_table_root / 'pangenome_feature.parquet')

        #self.build_user_genome_features_table()
        #self.build_pangenome_genome_features_table()
         #self.build_phenotype_kegg_table()
        # model
        # ontology
        #

    def build_genome_table(self):
        table_name = 'genome'
        conn = sqlite3.connect(str(self.root_pangenome.out_sqlite3_file))
        cur = conn.cursor()
        cur.execute(f"DROP TABLE IF EXISTS {table_name};")
        cur.execute(f"""
                CREATE TABLE IF NOT EXISTS {table_name} (
                    genome TEXT NOT NULL,
                    gtdb_taxonomy TEXT,
                    ncbi_taxonomy TEXT,
                    ncbi_taxid INT,
                    checkm_completeness REAL,
                    checkm_contamination REAL,
                    size INT,
                    kind TEXT NOT NULL,
                    PRIMARY KEY (genome)
                );
                """)

        clade_ids = set()
        member_ids = {o[0] for o in self.df_members.select("genome_id").rows()}
        with open(self.root_genome.ani_kepangenomes_json, 'r') as fh:
            clades = json.load(fh)
        for _, v in clades.items():
            clade_ids |= set(v)

        gtdb_ids = clade_ids | member_ids
        print(gtdb_ids)

        path_gtdb_ar = Path('/data/reference_data/gtdb/ar53_metadata_r214.tsv')
        path_gtdb_bac = Path('/data/reference_data/gtdb/bac120_metadata_r214.tsv')
        if path_gtdb_bac.exists() and path_gtdb_ar.exists():
            ldf_gtdb_metadata = pl.scan_csv([path_gtdb_ar, path_gtdb_bac], separator='\t')
            res = ldf_gtdb_metadata.filter(pl.col("accession").is_in(gtdb_ids)).select(
                pl.col("accession").alias("genome"),
                pl.col("gtdb_taxonomy"),
                pl.col("ncbi_taxonomy"),
                pl.col("ncbi_taxid"),
                pl.col("checkm_completeness"),
                pl.col("checkm_contamination"),
                pl.when(pl.col("accession").is_in(member_ids)).then(
                    pl.lit("clade_member")).otherwise(pl.lit("clade")).alias("kind"),
                #  pl.col("gc_count"),
                #  pl.col("gc_percentage"),
                #  pl.col("l50_contigs"),
                #  pl.col("l50_scaffolds"),
                #  pl.col("longest_contig"),
                #  pl.col("longest_scaffold"),
                pl.col("genome_size").alias("size"),

            ).collect()

            #print(res)

            res.to_pandas().to_sql(table_name, conn, if_exists="append", index=False)
        else:
            print("error: gtdb metadata not found unable to add pangenome members data")

        _data_input_genomes = {
            'genome': [],
            'gtdb_taxonomy': [],
            'ncbi_taxonomy': [],
            'ncbi_taxid': [],
            'checkm_completeness': [],
            'checkm_contamination': [],
            'size': [],
            'kind': [],
        }

        for input_genome in self.input_genomes:
            _data_input_genomes['genome'].append(input_genome)
            _data_input_genomes['gtdb_taxonomy'].append(None)
            _data_input_genomes['ncbi_taxonomy'].append(None)
            _data_input_genomes['ncbi_taxid'].append(None)
            _data_input_genomes['checkm_completeness'].append(None)
            _data_input_genomes['checkm_contamination'].append(None)
            _data_input_genomes['size'].append(None)
            _data_input_genomes['kind'].append('user')

        pl.DataFrame(_data_input_genomes).to_pandas().to_sql(
            table_name, conn, if_exists="append", index=False)

        conn.commit()
        conn.close()

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
                        feature_ontology_terms[feature_id][ontology_term].add(
                            QueryOntologyLocal.clean_bakta_value(
                                ontology_term, value))  # nasty hack to clean KEGG:, COG:, uniref prefix

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

        if self.root_genome.ani_phenotypes_json.exists():
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

        df_clusters = None
        path_pangenome_cluster = self.root_pangenome.root / 'pangenome_cluster_with_mmseqs.parquet'
        if path_pangenome_cluster.exists():
            df_clusters = pl.read_parquet(path_pangenome_cluster)
        for f in os.listdir(str(path_genome_dir)):
            if f.endswith('.faa'):
                path_genome_faa = path_genome_dir / f
                genome_id = path_genome_faa.name[:-4]  # strip .faa
                genome = MSGenome.from_fasta(str(path_genome_faa))
                feature_ontology_terms = genome_id_to_ontologies.get(genome_id, {})

                d_feature_cluster = {}
                if df_clusters is not None:
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
                    protein_hash = ProteinSequence(feature.seq).hash_value if feature.seq else None

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
        faa_to_process = set()
        for f in os.listdir(str(input_genome_dir)):
            if f.endswith('.faa'):
                path_genome_faa = input_genome_dir / f
                genome_id = path_genome_faa.name[:-4]  # strip .faa
                if genome_id in self.input_genomes:
                    path_genome_tsv = input_genome_dir / f'{genome_id}_genome.tsv'
                    if path_genome_tsv.exists():
                        print('build_user_genome_features_table - found', genome_id)
                        feature_ontology_terms = self.collect_annotation(genome_id, input_genome_dir)
                        for o in feature_ontology_terms.values():
                            all_ontology_terms |= set(o)
                        genome_id_to_ontologies[genome_id] = feature_ontology_terms

                    faa_to_process.add(path_genome_faa)
                else:
                    print(f'skip {genome_id}. Does not belong to this pangenome')

        print(f'found these ontology terms: {all_ontology_terms}')
        for term in all_ontology_terms:
            data[f'ontology_{term}'] = []

        for path_genome_faa in faa_to_process:
            genome_id = path_genome_faa.name[:-4]  # strip user_<genome_id>.faa
            path_genome_tsv = input_genome_dir / f'{genome_id}_genome.tsv'
            path_genome_pangenome_profile = input_genome_dir / f'{genome_id}_pangenome_profile.tsv'
            if path_genome_tsv.exists():
                d_pangenome_profile = {}
                if path_genome_pangenome_profile.exists():
                    print(f'found pangenome profile: {path_genome_pangenome_profile}')
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
                        d_pangenome_profile.get(feature_id, {}).get('pangenome_cluster'))
                    data['pangenome_is_core'].append(
                        d_pangenome_profile.get(feature_id, {}).get('is_core')
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

    def build_pangenome_genome_features_table(self, df):
        conn = sqlite3.connect(str(self.root_pangenome.out_sqlite3_file))
        cur = conn.cursor()
        cur.execute("DROP TABLE IF EXISTS pangenome_feature;")
        sql_create_table = """
        CREATE TABLE IF NOT EXISTS pangenome_feature (
                    genome TEXT NOT NULL,
                    contig TEXT NOT NULL,
                    feature_id TEXT NOT NULL,
                    start INTEGER NOT NULL,
                    end INTEGER NOT NULL,
                    strand TEXT NOT NULL,
                    type TEXT NOT NULL,
                    protein_sequence TEXT,
                    protein_sequence_hash TEXT,
                    cluster TEXT NOT NULL,
                    is_core INTEGER NOT NULL,
        """
        ontology_cols = {col for col in df.columns if col.startswith('ontology_')}
        for col in ontology_cols:
            sql_create_table += f"{col} TEXT,\n"
        sql_create_table += "PRIMARY KEY (genome, contig, feature_id));"
        cur.execute(sql_create_table)

        df.to_pandas().to_sql("pangenome_feature", conn, if_exists="append", index=False)

        conn.commit()
        conn.close()

    def build_genome_phenotype(self):
        conn = sqlite3.connect(str(self.root_pangenome.out_sqlite3_file))
        cur = conn.cursor()
        cur.execute("DROP TABLE IF EXISTS genome_phenotype;")
        sql_create_table = """
        CREATE TABLE genome_phenotype (
          genome_id                  TEXT NOT NULL,
          phenotype_id               TEXT NOT NULL,
          phenotype_name             TEXT,
          class                      TEXT,     -- e.g. P / N / A / C
          simulated_objective        REAL,
          observed_objective         REAL,
          gap_count                  INTEGER,
          gapfilled_reactions        TEXT,     -- semicolon-separated reaction IDs
          reaction_count             INTEGER,
          transports_added           TEXT,     -- e.g. cpd00971_c0 (may be empty)
          closest_experimental_data  TEXT,
          source                     TEXT,     -- e.g. pangenome
        
          PRIMARY KEY (genome_id, phenotype_id)
        );
        """
        cur.execute(sql_create_table)

        path_to_data = self.root_genome.root / 'phenotypes' / 'genome_phenotypes.tsv'
        if path_to_data.exists():
            ldf = pl.scan_csv(path_to_data, separator='\t')
            genome_ids = set(
                ldf.select("genome_id").unique().collect().get_column("genome_id").to_list())
            excluded = genome_ids - self.filter_genome_ids
            print('genome_phenotype excluded', excluded)

            df_filter = ldf.filter(pl.col("genome_id").is_in(self.filter_genome_ids)).collect()
            df_filter.to_pandas().to_sql("genome_phenotype", conn, if_exists="append", index=False)

        conn.commit()
        conn.close()

    def build_gene_phenotype(self):
        conn = sqlite3.connect(str(self.root_pangenome.out_sqlite3_file))
        cur = conn.cursor()
        cur.execute("DROP TABLE IF EXISTS gene_phenotype;")
        sql_create_table = """
        CREATE TABLE gene_phenotype (
          genome_id              TEXT NOT NULL,
          gene_id                TEXT NOT NULL,
          phenotype_id           TEXT NOT NULL,
          phenotype_name         TEXT,
          association_sources    TEXT,   -- e.g. model_prediction
          model_pred_reactions   TEXT,   -- semicolon-separated reaction IDs
          model_pred_max_flux    REAL,
          fitness_match          TEXT,   -- e.g. no_fitness_ortholog
          fitness_max            REAL,
          fitness_min            REAL,
          fitness_avg            REAL,
          fitness_count          INTEGER,
          essentiality_fraction  REAL,
        
          PRIMARY KEY (genome_id, gene_id, phenotype_id)
        );
        """
        cur.execute(sql_create_table)

        path_to_data = self.root_genome.root / 'phenotypes' / 'gene_phenotypes.tsv'
        if path_to_data.exists():
            ldf = pl.scan_csv(path_to_data, separator='\t')
            genome_ids = set(
                ldf.select("genome_id").unique().collect().get_column("genome_id").to_list())
            excluded = genome_ids - self.filter_genome_ids
            print('gene_phenotype excluded', excluded)

            df_filter = ldf.filter(pl.col("genome_id").is_in(self.filter_genome_ids)).collect()
            df_filter.to_pandas().to_sql("gene_phenotype", conn, if_exists="append", index=False)

        conn.commit()
        conn.close()

    def build_gene_reaction_data(self):
        table_name = 'genome_gene_reaction_essentially_test'
        conn = sqlite3.connect(str(self.root_pangenome.out_sqlite3_file))
        cur = conn.cursor()
        cur.execute(f"DROP TABLE IF EXISTS {table_name};")
        sql_create_table = f"""
        CREATE TABLE IF NOT EXISTS {table_name} (
           genome_id           TEXT NOT NULL,
           gene_id             TEXT NOT NULL,
           reaction            TEXT,     -- semicolon-separated reaction IDs
           rich_media_flux     REAL,
           rich_media_class    TEXT,     -- e.g. blocked / variable / essential
           minimal_media_flux  REAL,
           minimal_media_class TEXT,
           PRIMARY KEY (genome_id, gene_id)
        );
        """
        cur.execute(sql_create_table)

        path_to_data = self.root_genome.root / 'models' / 'gene_reaction_data.tsv'
        if path_to_data.exists():
            ldf = pl.scan_csv(path_to_data, separator='\t')
            genome_ids = set(
                ldf.select("genome_id").unique().collect().get_column("genome_id").to_list())
            excluded = genome_ids - self.filter_genome_ids
            print(f'{table_name} excluded: {excluded}')

            df_filter = ldf.filter(pl.col("genome_id").is_in(self.filter_genome_ids)).collect()
            df_filter.to_pandas().to_sql(table_name, conn, if_exists="append", index=False)

        conn.commit()
        conn.close()

    def build_input_genome_reactions(self):
        conn = sqlite3.connect(str(self.root_pangenome.out_sqlite3_file))
        cur = conn.cursor()
        cur.execute("DROP TABLE IF EXISTS genome_reaction;")
        sql_create_table = """
        CREATE TABLE IF NOT EXISTS genome_reaction (
          genome_id            TEXT,
          reaction_id          TEXT,
          genes                TEXT,
          equation_names       TEXT,
          equation_ids         TEXT,
          directionality       TEXT,
          upper_bound          REAL,
          lower_bound          REAL,
          gapfilling_status    TEXT,
          rich_media_flux      REAL,
          rich_media_class     TEXT,
          minimal_media_flux   REAL,
          minimal_media_class  TEXT,
          PRIMARY KEY (genome_id, reaction_id)
        );
        """
        cur.execute(sql_create_table)

        path_data = self.root_genome.root / 'models' / 'genome_reactions.tsv'
        if path_data.exists():
            ldf = pl.scan_csv(path_data, separator='\t')
            genome_ids = set(
                ldf.select("genome_id").unique().collect().get_column("genome_id").to_list())
            excluded = genome_ids - self.filter_genome_ids
            print('genome_reaction excluded', excluded)

            df_filter = ldf.filter(pl.col("genome_id").is_in(self.filter_genome_ids)).collect()
            df_filter.to_pandas().to_sql("genome_reaction", conn, if_exists="append", index=False)

        conn.commit()
        conn.close()

    def build_model_performance(self):
        conn = sqlite3.connect(str(self.root_pangenome.out_sqlite3_file))
        cur = conn.cursor()
        cur.execute("DROP TABLE IF EXISTS model_performance;")
        sql_create_table = """
        CREATE TABLE model_performance (
          genome_id                  TEXT NOT NULL,
          taxonomy                   TEXT,
          false_positives            INTEGER,
          false_negatives            INTEGER,
          true_positives             INTEGER,
          true_negatives             INTEGER,
          accuracy                   REAL,
          positive_growth            INTEGER,
          negative_growth            INTEGER,
          avg_positive_growth_gaps   REAL,
          avg_negative_growth_gaps   REAL,
          closest_user_genomes       TEXT,
          source                     TEXT,   -- e.g. pangenome / user / experiment
        
          PRIMARY KEY (genome_id)
        );
        """
        cur.execute(sql_create_table)

        path_data = self.root_genome.root / 'phenotypes' / 'model_performance.tsv'
        if path_data.exists():
            ldf = pl.scan_csv(path_data, separator='\t')
            genome_ids = set(
                ldf.select("genome_id").unique().collect().get_column("genome_id").to_list())
            excluded = genome_ids - self.filter_genome_ids
            print('model_performance excluded', excluded)

            df_filter = ldf.filter(pl.col("genome_id").is_in(self.filter_genome_ids)).collect()
            df_filter.to_pandas().to_sql("model_performance", conn, if_exists="append", index=False)

        conn.commit()
        conn.close()

    def build_table_media_composition(self):
        conn = sqlite3.connect(str(self.root_pangenome.out_sqlite3_file))
        cur = conn.cursor()
        cur.execute("DROP TABLE IF EXISTS media_composition;")
        sql_create_table = """
                CREATE TABLE media_composition (
                  media_id              TEXT NOT NULL,
                  compound_id           TEXT,
                  max_uptake            REAL,
                  compound_name         TEXT,

                  PRIMARY KEY (media_id, compound_id)
                );
                """
        cur.execute(sql_create_table)

        path_data = self.root_genome.root / 'models' / 'media_compositions.tsv'
        if path_data.exists():
            pl.read_csv(path_data, separator='\t').to_pandas().to_sql(
                "media_composition", conn, if_exists="append", index=False)

        conn.commit()
        conn.close()

    def build_phenotype_kegg_table(self):
        pass
