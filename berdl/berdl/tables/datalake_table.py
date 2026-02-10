import json
import sqlite3
from berdl.genome_paths import GenomePaths
from berdl.pangenome.paths_pangenome import PathsPangenome
from modelseedpy import MSGenome
from pathlib import Path


class DatalakeTableBuilder:

    def __init__(self, root_genome: GenomePaths, root_pangenome: PathsPangenome):
        self.root_genome = root_genome
        self.root_pangenome = root_pangenome

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

    def build_ani_table(self):
        with open(self.root_genome.ani_kepangenomes_json) as fh:
            data_ani_clade = json.load(fh)
        with open(self.root_genome.ani_fitness_json) as fh:
            data_ani_fitness = json.load(fh)
        with open(self.root_genome.ani_phenotypes_json) as fh:
            data_ani_phenotypes = json.load(fh)

        conn = sqlite3.connect(str(self.root_pangenome.out_sqlite3_file))
        cur = conn.cursor()
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
        for genome_user, ani_matches in data_ani_clade.items():
            for genome_ani, (ani, af1, af2) in ani_matches.items():
                cur.execute(
                    "INSERT INTO ani (genome1, genome2, ani, af1, af2) VALUES (?, ?, ?, ?, ?)",
                    (f'user_{genome_user}', genome_ani, ani, af1, af2)
                )

        for genome_user, ani_matches in data_ani_fitness.items():
            for genome_ani, (ani, af1, af2) in ani_matches.items():
                cur.execute(
                    "INSERT INTO ani (genome1, genome2, ani, af1, af2) VALUES (?, ?, ?, ?, ?)",
                    (f'user_{genome_user}', genome_ani, ani, af1, af2)
                )

        for genome_user, ani_matches in data_ani_phenotypes.items():
            for genome_ani, (ani, af1, af2) in ani_matches.items():
                cur.execute(
                    "INSERT INTO ani (genome1, genome2, ani, af1, af2) VALUES (?, ?, ?, ?, ?)",
                    (f'user_{genome_user}', genome_ani, ani, af1, af2)
                )

        conn.commit()
        conn.close()

    def build_user_genome_features_table(self):
        conn = sqlite3.connect(str(self.root_pangenome.out_sqlite3_file))
        cur = conn.cursor()
        cur.execute("""
                CREATE TABLE IF NOT EXISTS user_feature (
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
        pass

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
        for feature in genome.feature:
            pass
        pass

    def build_phenotype_kegg_table(self):
        pass

