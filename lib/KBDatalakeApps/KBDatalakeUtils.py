import sys
import os
import re
import sqlite3
import json
import uuid
import shutil
import hashlib
from os import path

sys.path = ["/deps/KBUtilLib/src","/deps/cobrakbase","/deps/ModelSEEDpy","/deps/cb_annotation_ontology_api"] + sys.path

# Import utilities with error handling
from kbutillib import KBAnnotationUtils, KBReadsUtils, KBGenomeUtils, MSReconstructionUtils, MSFBAUtils, MSBiochemUtils

import pandas as pd
import cobra
from modelseedpy.core.msgenome import MSGenome, MSFeature
from modelseedpy import AnnotationOntology, MSPackageManager, MSMedia, MSModelUtil, MSBuilder, MSATPCorrection, MSGapfill, MSGrowthPhenotype, MSGrowthPhenotypes, ModelSEEDBiochem, MSExpression
from modelseedpy.core.mstemplate import MSTemplateBuilder
#from berdl import OntologyEnrichment

# BERDL query modules (Jose P. Faria)
"""
try:
    from berdl.berdl import QueryPangenomeBERDL, OntologyEnrichment
except ImportError:
    # Fallback for when berdl module is not installed
    QueryPangenomeBERDL = None
    OntologyEnrichment = None
"""

class KBDataLakeUtils(KBGenomeUtils, MSReconstructionUtils, MSFBAUtils):
    def __init__(self,reference_path,module_path,**kwargs):
        super().__init__(
                name="KBDataLakeUtils",
                **kwargs
        )
        self.module_path = module_path
        self.reference_path = reference_path

    def run_user_genome_to_tsv(self,genome_ref,output_filename):
        # Load genome object into object_hash
        self.build_genome_tsv(genome_ref,output_filename)

    def pipeline_process_arguments_into_user_genome_table(self,output_filename):
        """
        Pipeline step for processing arguments.
        Translates input reference list (genomes or genome sets) into a table
        of user genomes and saves as user_genomes.tsv in self.directory.
        """
        rows = []
        input_refs = self.app_parameters.get('input_refs', [])

        for ref in input_refs:
            obj = self.get_object(ref)
            if obj is None:
                print(f"Warning: Could not retrieve object {ref}, skipping")
                continue

            obj_type = obj.get("info", [None]*3)[2] or ""
            obj_data = obj.get("data", {})

            if "GenomeSet" in obj_type:
                # GenomeSet: extract individual genome refs
                genome_refs = []
                if "elements" in obj_data:
                    for _, elem in obj_data["elements"].items():
                        genome_refs.append(elem.get("ref", ""))
                elif "items" in obj_data:
                    for item in obj_data["items"]:
                        genome_refs.append(item.get("ref", ""))
                for gref in genome_refs:
                    if gref:
                        row = self._extract_genome_metadata(gref)
                        if row:
                            rows.append(row)
            elif "Genome" in obj_type:
                row = self._extract_genome_metadata(ref, obj=obj)
                if row:
                    rows.append(row)
            else:
                print(f"Warning: Object {ref} has unsupported type {obj_type}, skipping")

        # Create DataFrame and save
        columns = [
            'genome_id', 'species_name', 'taxonomy', 'genome_ref',
            'assembly_ref', 'genome_type', 'genome_source_id',
            'genome_source_name', 'num_contigs', 'num_proteins',
            'num_noncoding_genes'
        ]
        df = pd.DataFrame(rows, columns=columns)
        os.makedirs(output_filename, exist_ok=True)
        df.to_csv(output_filename, sep='\t', index=False)
        print(f"Saved {len(df)} genomes to {output_filename}")

    def _extract_genome_metadata(self, ref, obj=None):
        """Extract metadata from a genome object for the user genomes table."""
        if obj is None:
            obj = self.get_object(ref)
        if obj is None:
            return None

        info = obj.get("info", [])
        data = obj.get("data", {})

        genome_id = info[1] if len(info) > 1 else ref
        species_name = data.get("scientific_name", "Unknown")
        taxonomy = data.get("taxonomy", "")
        genome_ref = f"{info[6]}/{info[0]}/{info[4]}" if len(info) > 6 else ref
        assembly_ref = data.get("assembly_ref", "")
        genome_type = data.get("domain", "Unknown")
        genome_source_id = data.get("source_id", "")
        genome_source_name = data.get("source", "")
        num_contigs = data.get("num_contigs", 0)

        num_proteins = 0
        num_noncoding = 0
        if "features" in data:
            for ftr in data["features"]:
                if ftr.get("protein_translation"):
                    num_proteins += 1
                else:
                    num_noncoding += 1
        if "non_coding_features" in data:
            num_noncoding += len(data["non_coding_features"])

        return [
            genome_id, species_name, taxonomy, genome_ref,
            assembly_ref, genome_type, genome_source_id,
            genome_source_name, num_contigs, num_proteins,
            num_noncoding
        ]

    def pipeline_run_pangenome_kberdl_query(self):
        """
        Pipeline step for running pangenome query against KBase BERDL.

        Uses SKANI results to find related pangenome data for each user genome.
        Queries BERDL for:
        - Clade membership information
        - Gene clusters for matched clades
        - ANI matrices for matched genomes

        Results are saved to self.directory/berdl_pangenome/

        Author: Jose P. Faria (jplfaria@gmail.com)
        """

        raise NotImplementedError("AI hallucinated")

        if QueryPangenomeBERDL is None:
            print("Warning: BERDL pangenome module not available, skipping")
            return

        # Get token from environment or app parameters
        token = os.environ.get('KBASE_AUTH_TOKEN') or self.app_parameters.get('token', '')
        if not token:
            print("Warning: No BERDL token available, skipping pangenome query")
            return

        skani_file = os.path.join(self.directory, "skani", "skani_pangenome.tsv")
        if not os.path.exists(skani_file):
            print(f"Warning: SKANI pangenome results not found at {skani_file}, skipping")
            return

        output_dir = os.path.join(self.directory, "berdl_pangenome")
        os.makedirs(output_dir, exist_ok=True)

        print("Running BERDL pangenome queries...")

        try:
            # Initialize BERDL query client
            qp = QueryPangenomeBERDL(token=token)

            # Load SKANI results
            skani_df = pd.read_csv(skani_file, sep='\t')
            print(f"  Loaded {len(skani_df)} SKANI hits")

            # Get unique reference genomes from SKANI hits
            if 'reference_genome' not in skani_df.columns:
                print("  Warning: No reference_genome column in SKANI results")
                return

            reference_genomes = skani_df['reference_genome'].unique().tolist()
            print(f"  Found {len(reference_genomes)} unique reference genomes")

            # Query clade information for each reference genome
            clade_results = []
            for ref_genome in reference_genomes[:50]:  # Limit to avoid timeout
                try:
                    clade_id = qp.get_member_representative(ref_genome)
                    clade_results.append({
                        'reference_genome': ref_genome,
                        'gtdb_species_clade_id': clade_id
                    })
                except Exception as e:
                    print(f"  Warning: Could not get clade for {ref_genome}: {e}")

            if clade_results:
                clade_df = pd.DataFrame(clade_results)
                clade_path = os.path.join(output_dir, "reference_clades.tsv")
                clade_df.to_csv(clade_path, sep='\t', index=False)
                print(f"  Saved {len(clade_results)} clade mappings to {clade_path}")

                # Get unique clades and query their members
                unique_clades = clade_df['gtdb_species_clade_id'].unique().tolist()
                print(f"  Querying {len(unique_clades)} unique clades...")

                all_clade_members = []
                for clade_id in unique_clades[:20]:  # Limit
                    try:
                        members_df = qp.get_clade_members(clade_id)
                        if not members_df.empty:
                            members_df['query_clade'] = clade_id
                            all_clade_members.append(members_df)
                    except Exception as e:
                        print(f"  Warning: Could not get members for clade {clade_id}: {e}")

                if all_clade_members:
                    combined_members = pd.concat(all_clade_members, ignore_index=True)
                    members_path = os.path.join(output_dir, "clade_members.tsv")
                    combined_members.to_csv(members_path, sep='\t', index=False)
                    print(f"  Saved {len(combined_members)} clade members to {members_path}")

            print("BERDL pangenome queries complete")

        except Exception as e:
            print(f"Error in BERDL pangenome query: {e}")

    def pipeline_run_ontology_term_kberdl_query(self, filename_datalake_db: str):
        """
        Pipeline step for running ontology term query against KBase BERDL.

        Extracts ontology term IDs (GO, EC, KEGG, COG, PFAM, SO) from genome
        data and enriches them with labels and definitions from BERDL API.

        Input sources (checked in order):
        1. SQLite database (berdl_tables.db) - genome_features table
        2. TSV files in genomes/ directory

        Results are saved to:
        - self.directory/ontology_terms.tsv (all enriched terms)
        - Adds ontology_terms table to SQLite database if it exists

        Author: Jose P. Faria (jplfaria@gmail.com)
        """
        if OntologyEnrichment is None:
            print("Warning: BERDL ontology module not available, skipping")
            return

        # Get token from environment or app parameters
        # FIXME: wrong token this is KBase auth token not for BERDL. (Please check if this works)
        token = os.environ.get('KBASE_AUTH_TOKEN') or self.app_parameters.get('token', '')
        if not token:
            print("Warning: No BERDL token available, skipping ontology enrichment")
            return

        # Try to load genome data from multiple sources
        genome_dataframes = []
        data_source = None

        # Source 1: SQLite database
        if os.path.exists(filename_datalake_db):
            try:
                conn = sqlite3.connect(str(filename_datalake_db))
                # Check if genome_features table exists
                cursor = conn.cursor()
                cursor.execute("SELECT name FROM sqlite_master WHERE type='table' AND name='genome_features'")
                if cursor.fetchone():
                    genome_df = pd.read_sql_query("SELECT * FROM genome_features", conn)
                    if not genome_df.empty:
                        genome_dataframes.append(genome_df)
                        data_source = "SQLite database"
                        print(f"  Loading genome data from SQLite: {len(genome_df)} features")
                conn.close()
            except Exception as e:
                print(f"  Warning: Could not read from SQLite: {e}")

        # Source 2: TSV files in genomes/ directory
        if not genome_dataframes:
            genomes_dir = os.path.join(self.directory, "genomes")
            if os.path.exists(genomes_dir):
                genome_files = [f for f in os.listdir(genomes_dir) if f.endswith('.tsv')]
                for genome_file in genome_files:
                    filepath = os.path.join(genomes_dir, genome_file)
                    try:
                        genome_df = pd.read_csv(filepath, sep='\t')
                        genome_dataframes.append(genome_df)
                    except Exception as e:
                        print(f"  Warning: Could not read {genome_file}: {e}")
                if genome_dataframes:
                    data_source = f"TSV files ({len(genome_dataframes)} genomes)"

        if not genome_dataframes:
            print("Warning: No genome data found (checked SQLite and TSV files), skipping")
            return

        print(f"Running ontology term enrichment from {data_source}...")

        try:
            # Initialize enrichment client
            # FIXME: which token BERDL or KBase ?
            enricher = OntologyEnrichment(token=token)

            # Collect all unique ontology terms across all genomes
            all_terms = set()

            for genome_df in genome_dataframes:
                try:
                    terms_by_type = enricher.extract_ontology_terms(genome_df)
                    for terms in terms_by_type.values():
                        all_terms.update(terms)
                except Exception as e:
                    print(f"  Warning: Could not extract terms from dataframe: {e}")

            if not all_terms:
                print("  No ontology terms found in genomes")
                return

            print(f"  Found {len(all_terms)} unique ontology terms")

            # Enrich all terms
            enriched_df = enricher.enrich_terms(list(all_terms))

            # Save enriched terms to TSV
            output_path = os.path.join(self.directory, "ontology_terms.tsv")
            enriched_df.to_csv(output_path, sep='\t', index=False)
            print(f"  Saved {len(enriched_df)} enriched terms to {output_path}")

            # Also save to SQLite database if it exists
            if os.path.exists(filename_datalake_db):
                try:
                    conn = sqlite3.connect(str(filename_datalake_db))
                    enriched_df.to_sql('ontology_terms', conn, if_exists='replace', index=False)
                    conn.close()
                    print(f"  Added 'ontology_terms' table to SQLite database")
                except Exception as e:
                    print(f"  Warning: Could not save to SQLite: {e}")

            # Summary by ontology type
            if not enriched_df.empty:
                print("\n  Enrichment summary:")
                for prefix in ['GO:', 'EC:', 'KEGG:', 'COG:', 'PFAM:', 'SO:']:
                    count = len(enriched_df[enriched_df['identifier'].str.startswith(prefix)])
                    if count > 0:
                        with_label = len(enriched_df[(enriched_df['identifier'].str.startswith(prefix)) &
                                                      (enriched_df['label'] != '')])
                        print(f"    {prefix[:-1]}: {count} terms, {with_label} with labels")

            print("Ontology term enrichment complete")

        except Exception as e:
            print(f"Error in ontology enrichment: {e}")

    def pipeline_save_annotated_genomes(self):
        """
        Pipeline step for saving annotated genomes back to KBase.
        Uses annotation ontology API to save RAST annotations to genome objects.
        """
        pass

    def pipeline_phenotype_tables(self, phenotypes_dir):
        """Build summary phenotype tables from individual simulation results."""
        accuracy_rows = []
        gene_pheno_rows = []
        pheno_gap_rows = []
        gapfill_rxn_rows = []

        for result_file in os.listdir(phenotypes_dir):
            if not result_file.endswith('_phenosim.json'):
                continue

            result_path = os.path.join(phenotypes_dir, result_file)
            with open(result_path) as f:
                result = json.load(f)

            model_id = result.get('model_id', result_file.replace('_phenosim.json', ''))

            # Genome accuracy table
            if 'accuracy' in result:
                accuracy_rows.append({
                    'genome_id': model_id,
                    **result['accuracy']
                })

            # Gene phenotype reactions
            for gp in result.get('gene_phenotype_reactions', []):
                gene_pheno_rows.append({'genome_id': model_id, **gp})

            # Phenotype gaps
            for pg in result.get('phenotype_gaps', []):
                pheno_gap_rows.append({'genome_id': model_id, **pg})

            # Gapfilled reactions
            for gr in result.get('gapfilled_reactions', []):
                gapfill_rxn_rows.append({'genome_id': model_id, **gr})

        # Save tables
        if accuracy_rows:
            pd.DataFrame(accuracy_rows).to_csv(
                os.path.join(phenotypes_dir, 'genome_accuracy.tsv'), sep='\t', index=False)
        if gene_pheno_rows:
            pd.DataFrame(gene_pheno_rows).to_csv(
                os.path.join(phenotypes_dir, 'genome_gene_phenotype_reactions.tsv'), sep='\t', index=False)
        if pheno_gap_rows:
            pd.DataFrame(pheno_gap_rows).to_csv(
                os.path.join(phenotypes_dir, 'genome_phenotype_gaps.tsv'), sep='\t', index=False)
        if gapfill_rxn_rows:
            pd.DataFrame(gapfill_rxn_rows).to_csv(
                os.path.join(phenotypes_dir, 'gapfilled_reactions.tsv'), sep='\t', index=False)

        print(f"Built phenotype summary tables in {phenotypes_dir}")

    def build_phenotype_tables(self, phenotypes_dir, phenosim_directory, experiment_data_file=None, phenoset_file=None):
        """
        Pipeline step for building phenotype experimental data, genome phenotype, and gene phenotype tables.

        Reads individual per-genome phenosim JSON files from phenosim_directory.
        Each file is named {genome_id}.json and contains:
          - data: dict keyed by cpd_id with objective_value, class, reactions,
                  gfreactions, gapfill_count, fluxes, etc.
          - data.summary: dict with accuracy, CP, CN, FP, FN, P, N

        Creates three TSV files:
        1. model_performance.tsv - Model accuracy metrics per genome
        2. genome_phenotypes.tsv - Phenotype predictions per genome
        3. gene_phenotypes.tsv - Gene-phenotype associations from gapfilled reactions

        Args:
            phenotypes_dir: Directory to save output TSV files
            phenosim_directory: Directory containing per-genome phenosim JSON files
            experiment_data_file: Path to experimental_data.json (optional)
            phenoset_file: Path to full_phenotype_set.json (optional)
        """
        # Default file paths
        if phenoset_file is None:
            phenoset_file = os.path.join(self.module_path, "data/full_phenotype_set.json")
        if experiment_data_file is None:
            experiment_data_file = os.path.join(self.module_path, "data/experimental_data.json")

        # Load phenotype set for compound names
        cpd_names = {}
        if os.path.exists(phenoset_file):
            with open(phenoset_file) as f:
                phenoset_data = json.load(f)
            for pheno in phenoset_data.get("phenotypes", []):
                cpd_names[pheno["id"]] = pheno.get("name", pheno["id"])

        # Load experimental data indexed by genome
        genome_experiment_data = {}
        if os.path.exists(experiment_data_file):
            with open(experiment_data_file) as f:
                experiment_data = json.load(f)
            for _, item_val in experiment_data.items():
                genome_id = item_val["genome_id"]
                if genome_id not in genome_experiment_data:
                    genome_experiment_data[genome_id] = {}
                genome_experiment_data[genome_id][item_val["cpd_id"]] = {
                    "growth": item_val["signal"],
                    "source": item_val["experiment"],
                    "chebi_id": item_val.get("chebi_id", "")
                }

        # Collect per-genome phenosim files
        phenosim_files = []
        if os.path.isdir(phenosim_directory):
            phenosim_files = [f for f in os.listdir(phenosim_directory) if f.endswith('.json')]
        if not phenosim_files:
            print(f"No phenosim JSON files found in {phenosim_directory}")
            return

        print(f"Processing {len(phenosim_files)} phenosim files from {phenosim_directory}")

        # Create output directory
        os.makedirs(phenotypes_dir, exist_ok=True)

        model_performance_rows = []
        genome_phenotype_rows = []
        gene_phenotype_rows = []

        for phenosim_file in phenosim_files:
            genome_id = phenosim_file.replace('.json', '')
            filepath = os.path.join(phenosim_directory, phenosim_file)
            try:
                with open(filepath) as f:
                    file_data = json.load(f)
            except (json.JSONDecodeError, IOError) as e:
                print(f"  Warning: Could not read {phenosim_file}: {e}")
                continue

            data = file_data.get("data", {})
            summary = data.get("summary", {})
            exp_data = genome_experiment_data.get(genome_id, {})

            # ===== TABLE 1: Model Performance row =====
            positive_growth_count = 0
            negative_growth_count = 0
            positive_gaps = []
            negative_gaps = []

            for cpd_id, cpd_data in data.items():
                if cpd_id == "summary" or not isinstance(cpd_data, dict):
                    continue
                obj_value = cpd_data.get("objective_value", 0) or 0
                gap_count = cpd_data.get("gapfill_count", 0) or 0

                if obj_value > 0.01:
                    positive_growth_count += 1
                    if gap_count > 0:
                        positive_gaps.append(gap_count)
                else:
                    negative_growth_count += 1
                    if gap_count > 0:
                        negative_gaps.append(gap_count)

            model_performance_rows.append({
                "genome_id": genome_id,
                "taxonomy": "",
                "false_positives": summary.get("FP", 0) or 0,
                "false_negatives": summary.get("FN", 0) or 0,
                "true_positives": summary.get("CP", 0) or 0,
                "true_negatives": summary.get("CN", 0) or 0,
                "positive_growth": positive_growth_count,
                "negative_growth": negative_growth_count,
                "avg_positive_growth_gaps": round(sum(positive_gaps) / len(positive_gaps), 4) if positive_gaps else 0,
                "avg_negative_growth_gaps": round(sum(negative_gaps) / len(negative_gaps), 4) if negative_gaps else 0,
                "closest_user_genomes": ""
            })

            # ===== TABLE 2 & 3: Per-phenotype rows =====
            # Track gene associations per genome for Table 3
            gene_phenotype_map = {}  # gene_id -> {cpd_id -> {reactions, scores}}

            for cpd_id, cpd_data in data.items():
                if cpd_id == "summary" or not isinstance(cpd_data, dict):
                    continue

                obj_value = cpd_data.get("objective_value", 0) or 0
                pred_class = cpd_data.get("class", "")
                reactions = cpd_data.get("reactions", [])
                gfreactions = cpd_data.get("gfreactions", {})
                gap_count = cpd_data.get("gapfill_count", 0) or 0
                fluxes = cpd_data.get("fluxes", {})

                # Find closest experimental data for this phenotype
                closest_exp = ""
                if cpd_id in exp_data:
                    exp_info = exp_data[cpd_id]
                    closest_exp = f"{exp_info['source']}:{'growth' if exp_info['growth'] else 'no_growth'}"

                # Table 2 row
                genome_phenotype_rows.append({
                    "genome_id": genome_id,
                    "phenotype_id": cpd_id,
                    "phenotype_name": cpd_names.get(cpd_id, cpd_id),
                    "growth_prediction": "positive" if obj_value > 0.01 else "negative",
                    "class": pred_class,
                    "objective_value": round(obj_value, 6),
                    "gap_count": gap_count,
                    "gapfilled_reactions": ";".join(gfreactions.keys()) if gfreactions else "",
                    "reaction_count": cpd_data.get("reaction_count", len(reactions)),
                    "closest_experimental_data": closest_exp
                })

                # Table 3: extract gene associations from gapfilled reactions
                # gfreactions format: {rxn_id: [direction, gene_or_null]}
                if gfreactions:
                    for rxn_id, rxn_info in gfreactions.items():
                        gene_id = rxn_info[1] if isinstance(rxn_info, list) and len(rxn_info) > 1 and rxn_info[1] else None
                        if gene_id:
                            if gene_id not in gene_phenotype_map:
                                gene_phenotype_map[gene_id] = {}
                            if cpd_id not in gene_phenotype_map[gene_id]:
                                gene_phenotype_map[gene_id][cpd_id] = {
                                    "reactions": set(),
                                    "fluxes": []
                                }
                            gene_phenotype_map[gene_id][cpd_id]["reactions"].add(rxn_id)
                            # Use flux as a fitness proxy if available
                            if rxn_id in fluxes:
                                gene_phenotype_map[gene_id][cpd_id]["fluxes"].append(abs(fluxes[rxn_id]))

            # Convert gene map to Table 3 rows
            for gene_id, phenotype_dict in gene_phenotype_map.items():
                for cpd_id, gdata in phenotype_dict.items():
                    scores = gdata["fluxes"]
                    gene_phenotype_rows.append({
                        "genome_id": genome_id,
                        "gene_id": gene_id,
                        "phenotype_id": cpd_id,
                        "associated_reactions": ";".join(sorted(gdata["reactions"])),
                        "max_fitness_score": round(max(scores), 6) if scores else 0,
                        "min_fitness_score": round(min(scores), 6) if scores else 0,
                        "avg_fitness_score": round(sum(scores) / len(scores), 6) if scores else 0
                    })

        # Write TSV files
        if model_performance_rows:
            df = pd.DataFrame(model_performance_rows)
            df.to_csv(os.path.join(phenotypes_dir, "model_performance.tsv"), sep="\t", index=False)
            print(f"Saved model_performance.tsv with {len(df)} genomes")

        if genome_phenotype_rows:
            df = pd.DataFrame(genome_phenotype_rows)
            df.to_csv(os.path.join(phenotypes_dir, "genome_phenotypes.tsv"), sep="\t", index=False)
            print(f"Saved genome_phenotypes.tsv with {len(df)} records")

        if gene_phenotype_rows:
            df = pd.DataFrame(gene_phenotype_rows)
            df.to_csv(os.path.join(phenotypes_dir, "gene_phenotypes.tsv"), sep="\t", index=False)
            print(f"Saved gene_phenotypes.tsv with {len(df)} records")
        else:
            print("No gene-phenotype associations found in phenosim data")

        print(f"Built phenotype tables in {phenotypes_dir}")

    def build_model_tables(self, models_dir):
        pass

    def pipeline_build_sqllite_db(self):
        """
        Pipeline step for building the SQLite database from all output data.
        """
        db_path = os.path.join(self.directory, "berdl_tables.db")
        if os.path.exists(db_path):
            os.remove(db_path)

        conn = sqlite3.connect(db_path)

        # Table 1: genome - from user_genomes.tsv
        genomes_file = os.path.join(self.directory, "user_genomes.tsv")
        if os.path.exists(genomes_file):
            df = pd.read_csv(genomes_file, sep='\t')
            df.rename(columns={
                'genome_id': 'id',
                'taxonomy': 'gtdb_taxonomy',
            }, inplace=True)
            if 'ncbi_taxonomy' not in df.columns:
                df['ncbi_taxonomy'] = df.get('gtdb_taxonomy', '')
            if 'n_features' not in df.columns:
                df['n_features'] = df.get('num_proteins', 0) + df.get('num_noncoding_genes', 0)
            genome_cols = ['id', 'gtdb_taxonomy', 'ncbi_taxonomy', 'num_contigs', 'n_features']
            existing_cols = [c for c in genome_cols if c in df.columns]
            df[existing_cols].to_sql('genome', conn, if_exists='replace', index=False)
            print("Added 'genome' table to SQLite DB")

        # Table 2: genome_ani - from skani results
        skani_dir = os.path.join(self.directory, "skani")
        if os.path.exists(skani_dir):
            ani_rows = []
            for skani_file in os.listdir(skani_dir):
                if not skani_file.endswith('.tsv'):
                    continue
                kind = skani_file.replace('skani_', '').replace('.tsv', '')
                skani_df = pd.read_csv(os.path.join(skani_dir, skani_file), sep='\t')
                for _, row in skani_df.iterrows():
                    ani_rows.append({
                        'genome1': row.get('genome_id', ''),
                        'genome2': row.get('reference_genome', ''),
                        'ani': row.get('ani_percentage', 0.0),
                        'af1': 0.0,
                        'af2': 0.0,
                        'kind': kind,
                    })
            if ani_rows:
                pd.DataFrame(ani_rows).to_sql('genome_ani', conn, if_exists='replace', index=False)
                print("Added 'genome_ani' table to SQLite DB")

        # Table 3: genome_features - from genome TSV files
        genomes_dir = os.path.join(self.directory, "genomes")
        if os.path.exists(genomes_dir):
            all_features = []
            for genome_file in os.listdir(genomes_dir):
                if not genome_file.endswith('.tsv'):
                    continue
                genome_id = genome_file.replace('.tsv', '')
                gdf = pd.read_csv(os.path.join(genomes_dir, genome_file), sep='\t')
                gdf['genome_id'] = genome_id
                # Rename columns to match schema
                col_map = {
                    'gene_id': 'feature_id',
                    'contig': 'contig_id',
                    'dna_sequence': 'sequence',
                    'functions': 'rast_function',
                }
                gdf.rename(columns=col_map, inplace=True)

                # Compute sequence hash
                if 'sequence' in gdf.columns:
                    gdf['sequence_hash'] = gdf['sequence'].apply(
                        lambda x: hashlib.md5(str(x).encode()).hexdigest() if pd.notna(x) and x else '')

                # Compute length from start/end
                if 'start' in gdf.columns and 'end' in gdf.columns:
                    gdf['length'] = abs(gdf['end'] - gdf['start'])

                all_features.append(gdf)

            if all_features:
                features_df = pd.concat(all_features, ignore_index=True)
                # Select columns that exist
                desired_cols = [
                    'genome_id', 'contig_id', 'feature_id', 'length', 'start',
                    'end', 'strand', 'sequence', 'sequence_hash', 'rast_function',
                    'rast_functions', 'protein_translation', 'ontology_terms'
                ]
                existing = [c for c in desired_cols if c in features_df.columns]
                features_df[existing].to_sql('genome_features', conn, if_exists='replace', index=False)
                print("Added 'genome_features' table to SQLite DB")

        # Table 4: phenotype tables
        phenotypes_dir = os.path.join(self.directory, "phenotypes")
        if os.path.exists(phenotypes_dir):
            pheno_tables = {
                'genome_accuracy.tsv': 'genome_accuracy',
                'genome_gene_phenotype_reactions.tsv': 'genome_gene_phenotype_reactions',
                'genome_phenotype_gaps.tsv': 'genome_phenotype_gaps',
                'gapfilled_reactions.tsv': 'gapfilled_reactions',
            }
            for tsv_file, table_name in pheno_tables.items():
                tsv_path = os.path.join(phenotypes_dir, tsv_file)
                if os.path.exists(tsv_path):
                    pdf = pd.read_csv(tsv_path, sep='\t')
                    pdf.to_sql(table_name, conn, if_exists='replace', index=False)
                    print(f"Added '{table_name}' table to SQLite DB")

        conn.close()
        print(f"SQLite database built at {db_path}")

    def pipeline_save_genometables_workspace_object(self):
        """
        Pipeline step for saving the genome tables workspace object.
        Left blank for future implementation.
        """
        pass

    def pipeline_save_kbase_report(self):
        """
        Pipeline step for saving the KBase report.
        Generates an HTML report and saves it to KBase workspace.
        """
        # Prepare report directory
        report_dir = os.path.join(self.directory, "report")
        os.makedirs(report_dir, exist_ok=True)

        # Build summary message
        genomes_file = os.path.join(self.directory, "user_genomes.tsv")
        models_dir = os.path.join(self.directory, "models")

        message = "Datalake pipeline completed.\n"
        if os.path.exists(genomes_file):
            df = pd.read_csv(genomes_file, sep='\t')
            message += f"Genomes processed: {len(df)}\n"

        if os.path.exists(models_dir):
            model_count = len([f for f in os.listdir(models_dir) if f.endswith('_model.json')])
            message += f"Models built: {model_count}\n"

        db_path = os.path.join(self.directory, "berdl_tables.db")
        if os.path.exists(db_path):
            message += f"SQLite database: {db_path}\n"

        # Prepare file links for downloadable files
        file_links = []
        if os.path.exists(db_path):
            file_links.append({
                'path': db_path,
                'name': 'berdl_tables.db',
                'description': 'SQLite database with all analysis results'
            })

        warnings = []

        # Use save_report_to_kbase from KBCallbackUtils if available,
        # otherwise use report_client directly
        try:
            report_client = self.report_client()

            # Copy HTML template if available
            html_source = '/kb/module/data/html'
            if os.path.exists(html_source):
                output_directory = os.path.join(self.directory, str(uuid.uuid4()))
                shutil.copytree(html_source, output_directory)

                dfu = self.dfu_client()
                shock_id = dfu.file_to_shock({
                    'file_path': output_directory,
                    'pack': 'zip'
                })['shock_id']

                html_report = [{
                    'shock_id': shock_id,
                    'name': 'index.html',
                    'label': 'BERDL Tables',
                    'description': 'BERDL Table Viewer'
                }]

                report_params = {
                    'message': message,
                    'warnings': warnings,
                    'workspace_name': self.workspace_name,
                    'objects_created': getattr(self, 'obj_created', []),
                    'html_links': html_report,
                    'direct_html_link_index': 0,
                    'html_window_height': 800,
                    'file_links': file_links,
                }
            else:
                report_params = {
                    'message': message,
                    'warnings': warnings,
                    'workspace_name': self.workspace_name,
                    'objects_created': getattr(self, 'obj_created', []),
                    'file_links': file_links,
                }

            report_info = report_client.create_extended_report(report_params)
            self.report_name = report_info['name']
            self.report_ref = report_info['ref']
            print(f"Report saved: {self.report_name} ({self.report_ref})")

        except Exception as e:
            print(f"Warning: Failed to save KBase report: {e}")

def run_phenotype_simulation(model_filename,output_filename,max_phenotypes,kbversion):
    """
    Worker function for parallel phenotype simulation.
    This runs in a separate process with its own memory space.

    Args:
        work_item: Dictionary with model_id, model_path, phenotypes_dir

    Returns:
        Dictionary with success status and simulation results
    """
    # Create safe model ID
    genome_id = os.path.splitext(os.path.basename(model_filename))[0]

    # Load model
    model = cobra.io.load_json_model(model_filename)
    mdlutl = MSModelUtil(model)

    #Loading the phenotype set from the reference path
    filename = "/kb/module/data/full_phenotype_set.json"
    with open(filename) as f:
        phenoset_data = json.load(f)
    #Setting max phenotypes if specified in the work item
    if max_phenotypes is not None:
        phenoset_data["phenotypes"] = phenoset_data["phenotypes"][:max_phenotypes]
    #Instantiating the phenotype set
    phenoset = MSGrowthPhenotypes.from_dict(phenoset_data)

    # Create a minimal util instance for this worker
    class PhenotypeWorkerUtil(MSReconstructionUtils,MSFBAUtils,MSBiochemUtils):
        def __init__(self, kbversion):
            super().__init__(name="PhenotypeWorkerUtil", kbversion=kbversion)
    pheno_util = PhenotypeWorkerUtil(kbversion=kbversion)

    # Get template for gapfilling
    template = pheno_util.get_template(pheno_util.templates["gn"], None)

    # Retrieve ATP test conditions from the model
    atpcorrection = MSATPCorrection(mdlutl)
    atp_tests = atpcorrection.build_tests()

    # Create gapfiller with ATP test conditions
    gapfiller = MSGapfill(
        mdlutl,
        default_gapfill_templates=[template],
        default_target='bio1',
        minimum_obj=0.01,
        test_conditions=[atp_tests[0]]
    )
    pheno_util.set_media(gapfiller.gfmodelutl, "KBaseMedia/Carbon-Pyruvic-Acid")

    # Prefilter gapfilling database with ATP test conditions
    #print("Prefiltering gapfilling database...")
    #gapfiller.prefilter()
    #print("Prefiltering complete")

    # Filter out mass imbalanced (MI) reactions from the gapfill model
    mi_blocked_count = 0
    reaction_scores = {}
    for rxn in gapfiller.gfmodelutl.model.reactions:
        # Extract the base ModelSEED reaction ID (rxnXXXXX) from the reaction ID
        if rxn.id not in mdlutl.model.reactions and rxn.id in gapfiller.gfpkgmgr.getpkg("GapfillingPkg").gapfilling_penalties:
            reaction_scores[rxn.id] = {
                "<": 10 * gapfiller.gfpkgmgr.getpkg("GapfillingPkg").gapfilling_penalties[rxn.id].get("reverse", 1),
                ">": 10 * gapfiller.gfpkgmgr.getpkg("GapfillingPkg").gapfilling_penalties[rxn.id].get("forward", 1)
            }
        ms_rxn_id = pheno_util.reaction_id_to_msid(rxn.id)
        if ms_rxn_id:
            ms_rxn = pheno_util.get_reaction_by_id(ms_rxn_id)
            if ms_rxn and hasattr(ms_rxn, 'status') and ms_rxn.status and "MI" in ms_rxn.status:
                reaction_scores[rxn.id] = {"<":1000,">":1000}
                # Check if status contains "MI" (mass imbalanced)
                if rxn.id not in mdlutl.model.reactions:
                    #If reaction is not in model, set bounds to 0
                    rxn.lower_bound = 0
                    rxn.upper_bound = 0
                    mi_blocked_count += 1

    # Run simulations with gapfilling for zero-growth phenotypes
    # Note: test_conditions=None since we already ran prefilter
    results = phenoset.simulate_phenotypes(
        mdlutl,
        add_missing_exchanges=True,
        gapfill_negatives=True,
        msgapfill=gapfiller,
        test_conditions=None,
        ignore_experimental_data=True,
        annoont=None,
        growth_threshold=0.01,
        #reaction_scores=reaction_scores
    )

    output_dir = os.path.dirname(output_filename)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
    with open(output_filename, "w") as f:
        json.dump(results, f, indent=4, skipkeys=True)

    return {"success": True, "genome_id": genome_id}


def run_model_reconstruction(input_filename, output_filename, classifier_dir,kbversion):
    worker_util = MSReconstructionUtils(kbversion=kbversion)

    # Clear MSModelUtil cache for this process
    MSModelUtil.mdlutls.clear()

    # Create safe model ID
    genome_id = os.path.splitext(os.path.basename(input_filename))[0]
    print(f'genome_id: {genome_id}')
    # Load features from genome TSV
    gene_df = pd.read_csv(input_filename, sep='\t')

    # Create MSGenome from features
    genome = MSGenome()
    genome.id = genome_id
    genome.scientific_name = genome_id

    # Detect TSV format based on columns present
    columns = set(gene_df.columns)
    is_simple_format = 'id' in columns and 'functions' in columns and 'protein_translation' not in columns
    print(f'TSV format detected: {"simple (id, function)" if is_simple_format else "full genome TSV"}')

    ms_features = []
    for _, gene in gene_df.iterrows():
        if is_simple_format:
            # Simple format: columns are 'id' and 'function'
            # 'function' is pipe-delimited list of RAST descriptions
            gene_id_val = gene.get('id', '')
            if pd.notna(gene_id_val) and gene_id_val:
                feature = MSFeature(str(gene_id_val), '')  # No protein sequence in simple format
                func_col = gene.get('function', '')
                if pd.notna(func_col) and func_col:
                    # Split on pipe delimiter for multiple RAST functions
                    for func_desc in str(func_col).split('|'):
                        func_desc = func_desc.strip()
                        if func_desc:
                            feature.add_ontology_term('RAST', func_desc)
                ms_features.append(feature)
        else:
            # Full genome TSV format
            protein = gene.get('protein_translation', '')
            gene_id_val = gene.get('gene_id', '')
            if pd.notna(protein) and protein:
                feature = MSFeature(str(gene_id_val), str(protein))

                # Parse plain functions column for RAST descriptions
                # Format: "Threonine synthase (EC 4.2.3.1)" or "Func1;Func2"
                func_col = gene.get('functions', '')
                if pd.notna(func_col) and func_col:
                    for func_desc in str(func_col).split(';'):
                        func_desc = func_desc.strip()
                        if func_desc:
                            feature.add_ontology_term('RAST', func_desc)

                # Parse Annotation:SSO column for SSO IDs
                # Format: SSO:nnnnn:description|MSRXN:rxn1,rxn2;SSO:mmmmm:desc2|rxn3
                sso_col = gene.get('Annotation:SSO', '')
                if pd.notna(sso_col) and sso_col:
                    for entry in str(sso_col).split(';'):
                        entry = entry.strip()
                        if not entry:
                            continue
                        term_part = entry.split('|')[0]
                        parts = term_part.split(':')
                        if len(parts) >= 2 and parts[0] == 'SSO':
                            sso_id = parts[1]
                            feature.add_ontology_term('SSO', sso_id)

                ms_features.append(feature)

    genome.add_features(ms_features)

    # Load classifier
    genome_classifier = worker_util.get_classifier(classifier_dir)

    # Build the model
    current_output, mdlutl = worker_util.build_metabolic_model(
        genome=genome,
        genome_classifier=genome_classifier,
        model_id=genome_id,
        model_name=genome_id,
        gs_template="auto",
        atp_safe=True,
        load_default_medias=True,
        max_gapfilling=10,
        gapfilling_delta=0,
    )

    if mdlutl is None:
        return {
            'success': False,
            'error': f"Model build returned None: {current_output.get('Comments', ['Unknown'])}"
        }

    model = mdlutl.model

    # Gapfill if media specified
    gf_rxns = 0
    growth = 'NA'
    gapfill_media = worker_util.get_media("KBaseMedia/Carbon-Pyruvic-Acid")
    gf_output, _, _, _ = worker_util.gapfill_metabolic_model(
        mdlutl=mdlutl,
        genome=genome,
        media_objs=[gapfill_media],
        templates=[model.template],
        atp_safe=True,
        objective='bio1',
        minimum_objective=0.01,
        gapfilling_mode="Sequential",
    )
    gf_rxns = gf_output.get('GS GF', 0)
    growth = gf_output.get('Growth', 'Unknown')

    # Save model
    cobra.io.save_json_model(model, output_filename+"_cobra.json")

    genome_class = current_output.get('Class', 'Unknown')
    core_gf = current_output.get('Core GF', 0)

    # Get minimal and rich media for analysis
    minimal_media = gapfill_media  # Use the gapfill media as minimal
    rich_media = worker_util.get_media("KBaseMedia/Complete")

    # Collect gapfilled reactions by category
    core_gf_rxns = []
    minimal_gf_rxns = []
    rich_gf_rxns = []

    # Get gapfilled reactions from model attributes
    if hasattr(mdlutl, 'attributes') and 'gapfilling' in mdlutl.attributes:
        for gf_entry in mdlutl.attributes['gapfilling']:
            media_id = gf_entry.get('media', {}).get('id', '')
            for rxn_id in gf_entry.get('reactions', []):
                if 'core' in media_id.lower() or gf_entry.get('target', '') == 'core':
                    if rxn_id not in core_gf_rxns:
                        core_gf_rxns.append(rxn_id)
                else:
                    if rxn_id not in minimal_gf_rxns:
                        minimal_gf_rxns.append(rxn_id)

    # Run FVA in rich media to identify essential gapfilled reactions
    all_gf_rxns = list(set(core_gf_rxns + minimal_gf_rxns))
    if all_gf_rxns and rich_media:
        mdlutl.pkgmgr.getpkg("KBaseMediaPkg").build_package(rich_media)
        try:
            fva_result = cobra.flux_analysis.flux_variability_analysis(
                model,
                reaction_list=[model.reactions.get_by_id(rxn_id) for rxn_id in all_gf_rxns if rxn_id in model.reactions],
                fraction_of_optimum=0.001
            )
            for rxn_id in all_gf_rxns:
                if rxn_id in fva_result.index:
                    min_flux = fva_result.loc[rxn_id, 'minimum']
                    max_flux = fva_result.loc[rxn_id, 'maximum']
                    # If flux bounds don't include zero, reaction is essential
                    if min_flux > 1e-6 or max_flux < -1e-6:
                        if rxn_id not in rich_gf_rxns:
                            rich_gf_rxns.append(rxn_id)
        except Exception as e:
            print(f"Warning: Rich media FVA failed: {e}")

    # Run pFBA and FVA analysis in minimal media
    minimal_pfba_fluxes = {}
    minimal_fva_classes = {}
    mdlutl.pkgmgr.getpkg("KBaseMediaPkg").build_package(minimal_media)
    try:
        # pFBA in minimal media
        pfba_solution = cobra.flux_analysis.pfba(model)
        if pfba_solution.status == 'optimal':
            minimal_pfba_fluxes = {rxn.id: pfba_solution.fluxes[rxn.id] for rxn in model.reactions}

        # FVA in minimal media for classification
        fva_minimal = cobra.flux_analysis.flux_variability_analysis(model, fraction_of_optimum=0.001)
        for rxn in model.reactions:
            if rxn.id in fva_minimal.index:
                min_f = fva_minimal.loc[rxn.id, 'minimum']
                max_f = fva_minimal.loc[rxn.id, 'maximum']
                if abs(min_f) < 1e-6 and abs(max_f) < 1e-6:
                    minimal_fva_classes[rxn.id] = 'blocked'
                elif min_f > 1e-6:
                    minimal_fva_classes[rxn.id] = 'essential_forward'
                elif max_f < -1e-6:
                    minimal_fva_classes[rxn.id] = 'essential_reverse'
                elif min_f < -1e-6 and max_f > 1e-6:
                    minimal_fva_classes[rxn.id] = 'reversible'
                elif max_f > 1e-6:
                    minimal_fva_classes[rxn.id] = 'forward_only'
                elif min_f < -1e-6:
                    minimal_fva_classes[rxn.id] = 'reverse_only'
                else:
                    minimal_fva_classes[rxn.id] = 'variable'
    except Exception as e:
        print(f"Warning: Minimal media analysis failed: {e}")

    # Run pFBA and FVA analysis in rich media
    rich_pfba_fluxes = {}
    rich_fva_classes = {}
    if rich_media:
        mdlutl.pkgmgr.getpkg("KBaseMediaPkg").build_package(rich_media)
        try:
            # pFBA in rich media
            pfba_solution = cobra.flux_analysis.pfba(model)
            if pfba_solution.status == 'optimal':
                rich_pfba_fluxes = {rxn.id: pfba_solution.fluxes[rxn.id] for rxn in model.reactions}

            # FVA in rich media for classification
            fva_rich = cobra.flux_analysis.flux_variability_analysis(model, fraction_of_optimum=0.001)
            for rxn in model.reactions:
                if rxn.id in fva_rich.index:
                    min_f = fva_rich.loc[rxn.id, 'minimum']
                    max_f = fva_rich.loc[rxn.id, 'maximum']
                    if abs(min_f) < 1e-6 and abs(max_f) < 1e-6:
                        rich_fva_classes[rxn.id] = 'blocked'
                    elif min_f > 1e-6:
                        rich_fva_classes[rxn.id] = 'essential_forward'
                    elif max_f < -1e-6:
                        rich_fva_classes[rxn.id] = 'essential_reverse'
                    elif min_f < -1e-6 and max_f > 1e-6:
                        rich_fva_classes[rxn.id] = 'reversible'
                    elif max_f > 1e-6:
                        rich_fva_classes[rxn.id] = 'forward_only'
                    elif min_f < -1e-6:
                        rich_fva_classes[rxn.id] = 'reverse_only'
                    else:
                        rich_fva_classes[rxn.id] = 'variable'
        except Exception as e:
            print(f"Warning: Rich media analysis failed: {e}")

    # Build metabolite list as JSON dictionaries
    metabolites_list = []
    for met in model.metabolites:
        met_dict = {
            'id': met.id,
            'name': met.name,
            'formula': met.formula,
            'charge': met.charge,
            'compartment': met.compartment,
            'annotation': dict(met.annotation) if met.annotation else {}
        }
        metabolites_list.append(met_dict)

    # Build reaction list as JSON dictionaries
    reactions_list = []
    for rxn in model.reactions:
        rxn_dict = {
            'id': rxn.id,
            'name': rxn.name,
            'reaction': rxn.reaction,
            'lower_bound': rxn.lower_bound,
            'upper_bound': rxn.upper_bound,
            'gene_reaction_rule': rxn.gene_reaction_rule,
            'subsystem': rxn.subsystem if hasattr(rxn, 'subsystem') else '',
            'annotation': dict(rxn.annotation) if rxn.annotation else {},
            'metabolites': {met.id: coef for met, coef in rxn.metabolites.items()}
        }
        reactions_list.append(rxn_dict)

    # Build the comprehensive output data
    output_data = {
        'success': True,
        'model_info': {
            'model_id': model.id,
            'num_reactions': len(model.reactions),
            'num_metabolites': len(model.metabolites),
            'num_genes': len(model.genes),
            'genome_class': genome_class,
            'core_gapfill': core_gf,
            'gs_gapfill': gf_rxns,
            'growth': growth
        },
        'metabolites': metabolites_list,
        'reactions': reactions_list,
        'gapfilled_reactions': {
            'core': core_gf_rxns,
            'minimal_media': minimal_gf_rxns,
            'rich_media_essential': rich_gf_rxns
        },
        'flux_analysis': {
            'minimal_media': {
                'pfba_fluxes': minimal_pfba_fluxes,
                'fva_classes': minimal_fva_classes
            },
            'rich_media': {
                'pfba_fluxes': rich_pfba_fluxes,
                'fva_classes': rich_fva_classes
            }
        }
    }

    # Save to JSON file
    with open(output_filename + "_data.json", 'w') as f:
        json.dump(output_data, f, indent=2)
def generate_ontology_tables(
    clade_folder: str,
    reference_data_path: str = "/data/reference_data",
    genome_features_table: str = "genome_features",
    output_folder_name: str = "ontology_data"
) -> bool:
    """
    Generate ontology tables for a clade folder.

    This function reads genome features from a db.sqlite file, maps RAST
    annotations to SEED roles, extracts EC numbers, enriches all ontology
    terms, and saves three output tables.

    Args:
        clade_folder: Path to the clade folder (e.g., /path/to/pangenome/s__Escherichia_coli)
                      Must contain a db.sqlite file with genome_features table.
        reference_data_path: Path to directory containing reference files:
                            - seed.json (RAST → seed.role mapping)
                            - statements.parquet (labels, definitions, relationships)
                            - kegg_ko_definitions.parquet
                            - cog_definitions.parquet
                            Default: /data/reference_data
        genome_features_table: Name of the table in db.sqlite to read features from.
                              Default: genome_features
        output_folder_name: Name of the output folder to create.
                           Default: ontology_data

    Returns:
        True on success, False on failure.

    Output files (in clade_folder/output_folder_name/):
        - ontology_terms.tsv: All ontology terms with labels and definitions
        - ontology_definition.tsv: Ontology prefix definitions
        - ontology_relationships.tsv: Term relationships (is_a, enables_reaction)
    """
    import sqlite3
    import re
    import time
    from pathlib import Path
    import pyarrow.parquet as pq

    clade_path = Path(clade_folder)
    db_path = clade_path / "db.sqlite"
    output_path = clade_path / output_folder_name

    # Check if db.sqlite exists
    if not db_path.exists():
        print(f"Warning: db.sqlite not found in {clade_folder}, skipping ontology generation")
        return False

    print(f"\n{'='*70}")
    print(f"Generating ontology tables for: {clade_folder}")
    print(f"{'='*70}")

    try:
        # =====================================================================
        # STEP 1: Load genome features from SQLite
        # =====================================================================
        print(f"\n1. Loading genome features from {db_path}...")
        conn = sqlite3.connect(str(db_path))

        # Check if table exists
        cursor = conn.cursor()
        cursor.execute(f"SELECT name FROM sqlite_master WHERE type='table' AND name='{genome_features_table}'")
        if not cursor.fetchone():
            print(f"   Warning: Table '{genome_features_table}' not found in db.sqlite")
            conn.close()
            return False

        genome_df = pd.read_sql_query(f"SELECT * FROM {genome_features_table}", conn)
        conn.close()
        print(f"   Loaded {len(genome_df)} features")
        print(f"   Columns: {list(genome_df.columns)[:10]}...")

        # =====================================================================
        # STEP 2: Initialize RASTSeedMapper for RAST → seed.role mapping
        # =====================================================================
        print("\n2. Loading RAST → seed.role mapper...")

        ref_path = Path(reference_data_path)
        seed_json_path = ref_path / "seed.json"

        mapper = None
        if seed_json_path.exists():
            mapper = RASTSeedMapper(str(seed_json_path))
        else:
            print(f"   Warning: seed.json not found at {seed_json_path}")
            print(f"   RAST → seed.role mapping will be skipped")

        # =====================================================================
        # STEP 3: Extract ontology terms from genome features
        # =====================================================================
        print("\n3. Extracting ontology terms...")

        terms_by_type = {
            'GO': set(), 'EC': set(), 'KEGG': set(),
            'COG': set(), 'PFAM': set(), 'SO': set(), 'seed.role': set()
        }

        # Patterns for extracting existing term IDs from annotation columns
        patterns = {
            'GO': re.compile(r'GO:\d+'),
            'EC': re.compile(r'EC:[\d\.-]+'),
            'KEGG': re.compile(r'(?:KEGG:)?K\d{5}'),
            'COG': re.compile(r'COG:(?:COG\d+|[A-Z])'),
            'PFAM': re.compile(r'(?:PFAM:)?PF\d+(?:\.\d+)?'),
            'SO': re.compile(r'SO:\d+'),
            'seed.role': re.compile(r'seed\.role:\d+'),
        }

        # Pattern for extracting EC from RAST function strings like "enzyme (EC 1.1.1.1)"
        ec_in_rast_pattern = re.compile(r'\(EC[:\s]*([\d\.-]+)\)')

        # Track RAST functions for seed.role mapping
        rast_functions = set()
        seed_role_to_label = {}  # seed.role ID -> RAST function label

        # Find the RAST function column
        rast_col = None
        for col in ['rast_function', 'rast_functions', 'functions', 'Annotation:SSO']:
            if col in genome_df.columns:
                rast_col = col
                break

        if rast_col:
            print(f"   Using RAST function column: {rast_col}")

        # Extract terms from all columns
        for col in genome_df.columns:
            for _, row in genome_df.iterrows():
                value = str(row.get(col, ''))
                if not value or value == 'nan':
                    continue

                # Extract existing ontology term IDs
                for ont_type, pattern in patterns.items():
                    matches = pattern.findall(value)
                    for match in matches:
                        # Normalize prefixes
                        if ont_type == 'KEGG' and not match.startswith('KEGG:'):
                            match = f'KEGG:{match}'
                        elif ont_type == 'PFAM' and not match.startswith('PFAM:'):
                            match = f'PFAM:{match}'
                        terms_by_type[ont_type].add(match)

                # Extract EC from RAST function strings
                if col == rast_col:
                    ec_matches = ec_in_rast_pattern.findall(value)
                    for ec_num in ec_matches:
                        terms_by_type['EC'].add(f'EC:{ec_num}')

                    # Collect RAST functions for seed.role mapping
                    if value and value != 'nan':
                        # Split multi-function annotations
                        for separator in [' / ', ' @ ', '; ']:
                            if separator in value:
                                parts = value.split(separator)
                                for part in parts:
                                    part = part.strip()
                                    if part:
                                        rast_functions.add(part)
                        if not any(sep in value for sep in [' / ', ' @ ', '; ']):
                            rast_functions.add(value)

        # =====================================================================
        # STEP 4: Map RAST functions to seed.role IDs
        # =====================================================================
        if mapper and rast_functions:
            print(f"\n4. Mapping {len(rast_functions)} RAST functions to seed.role IDs...")

            mapped_count = 0
            for rast_func in rast_functions:
                # Get all matching seed.role IDs for this function
                mappings = mapper.map_all_annotations(rast_func)
                for matched_part, seed_id in mappings:
                    if seed_id:
                        terms_by_type['seed.role'].add(seed_id)
                        seed_role_to_label[seed_id] = matched_part
                        mapped_count += 1

            print(f"   Mapped {mapped_count} RAST functions to {len(terms_by_type['seed.role'])} unique seed.role IDs")
        else:
            print("\n4. Skipping RAST → seed.role mapping (no mapper or no RAST functions)")

        # Summary of extracted terms
        total_terms = sum(len(terms) for terms in terms_by_type.values())
        print(f"\n   Total unique terms: {total_terms}")
        for ont_type, terms in sorted(terms_by_type.items()):
            if terms:
                print(f"     {ont_type}: {len(terms)}")

        if total_terms == 0:
            print("   Warning: No ontology terms found in genome features")
            os.makedirs(output_path, exist_ok=True)
            pd.DataFrame(columns=['ontology_prefix', 'identifier', 'label', 'definition']).to_csv(
                output_path / 'ontology_terms.tsv', sep='\t', index=False)
            pd.DataFrame(columns=['ontology_prefix', 'definition']).to_csv(
                output_path / 'ontology_definition.tsv', sep='\t', index=False)
            pd.DataFrame(columns=['subject', 'predicate', 'object']).to_csv(
                output_path / 'ontology_relationships.tsv', sep='\t', index=False)
            return True

        # =====================================================================
        # STEP 5: Enrich terms from local parquet files
        # =====================================================================
        print("\n5. Enriching terms from local parquet files...")

        statements_path = ref_path / "statements.parquet"
        kegg_path = ref_path / "kegg_ko_definitions.parquet"
        cog_path = ref_path / "cog_definitions.parquet"

        enriched_terms = []
        statements_df = None

        # Collect all terms that need enrichment from statements.parquet
        berdl_terms = list(terms_by_type['GO'] | terms_by_type['EC'] |
                          terms_by_type['SO'] | terms_by_type['PFAM'] |
                          terms_by_type['seed.role'])

        if berdl_terms and statements_path.exists():
            print(f"   Loading statements.parquet...")
            statements_df = pq.read_table(statements_path).to_pandas()
            print(f"   Loaded {len(statements_df)} statements")

            # Filter to relevant subjects and predicates for labels/definitions
            mask = (
                statements_df['subject'].isin(berdl_terms) &
                statements_df['predicate'].isin(['rdfs:label', 'IAO:0000115'])
            )
            filtered = statements_df[mask]

            # Build lookup dict
            term_info = {}
            for _, row in filtered.iterrows():
                subj = row['subject']
                pred = row['predicate']
                val = row['value'] if 'value' in row else ''

                if subj not in term_info:
                    term_info[subj] = {'label': '', 'definition': ''}

                if pred == 'rdfs:label':
                    term_info[subj]['label'] = val
                elif pred == 'IAO:0000115':
                    term_info[subj]['definition'] = val

            # Add to enriched_terms
            for term_id in berdl_terms:
                prefix = term_id.split(':')[0]
                info = term_info.get(term_id, {'label': '', 'definition': ''})

                # For seed.role, use the RAST function as label if no label found
                label = info['label']
                if prefix == 'seed.role' and not label and term_id in seed_role_to_label:
                    label = seed_role_to_label[term_id]

                enriched_terms.append({
                    'ontology_prefix': prefix,
                    'identifier': term_id,
                    'label': label,
                    'definition': info['definition']
                })

            print(f"   Enriched {len([t for t in enriched_terms if t['label']])} terms with labels")

        # Enrich KEGG from kegg_ko_definitions.parquet
        kegg_terms = list(terms_by_type['KEGG'])
        if kegg_terms and kegg_path.exists():
            print(f"   Loading kegg_ko_definitions.parquet...")
            kegg_df = pq.read_table(kegg_path).to_pandas()
            ko_lookup = dict(zip(kegg_df['ko_id'], kegg_df['definition']))

            for ko_id in kegg_terms:
                k_num = ko_id.replace('KEGG:', '')
                definition = ko_lookup.get(k_num, '')
                label = re.sub(r'\s*\[EC:[^\]]+\]', '', definition).strip() if definition else ''

                enriched_terms.append({
                    'ontology_prefix': 'KEGG',
                    'identifier': ko_id,
                    'label': label,
                    'definition': definition
                })

        # Enrich COG from cog_definitions.parquet
        cog_terms = list(terms_by_type['COG'])
        if cog_terms and cog_path.exists():
            print(f"   Loading cog_definitions.parquet...")
            cog_df = pq.read_table(cog_path).to_pandas()
            cog_lookup = {row['cog_id']: row for _, row in cog_df.iterrows()}

            for cog_id in cog_terms:
                raw_id = cog_id.replace('COG:', '')
                info = cog_lookup.get(raw_id, {})

                enriched_terms.append({
                    'ontology_prefix': 'COG',
                    'identifier': cog_id,
                    'label': info.get('name', '') if isinstance(info, dict) else '',
                    'definition': info.get('pathway', '') if isinstance(info, dict) else ''
                })

        # =====================================================================
        # STEP 6: Extract relationships from statements.parquet
        # =====================================================================
        print("\n6. Extracting ontology relationships...")

        relationships = []
        all_term_ids = set()
        for terms in terms_by_type.values():
            all_term_ids.update(terms)

        seed_reaction_terms = set()

        if statements_df is not None:
            # Look for is_a (GO) and enables_reaction (seed.role -> seed.reaction)
            relevant_predicates = {
                'rdfs:subClassOf',  # is_a hierarchy
                '<https://modelseed.org/ontology/enables_reaction>',  # seed.role → reaction
            }

            # Filter for our terms and relevant predicates
            mask = (
                statements_df['subject'].isin(all_term_ids) &
                statements_df['predicate'].isin(relevant_predicates)
            )
            rel_df = statements_df[mask]

            # Clean predicates and add to relationships
            predicate_labels = {
                'rdfs:subClassOf': 'is_a',
                '<https://modelseed.org/ontology/enables_reaction>': 'enables_reaction',
            }

            for _, row in rel_df.iterrows():
                subj = row['subject']
                pred = row['predicate']
                obj = row['object']

                # Skip self-referential or blank nodes
                if subj == obj or str(obj).startswith('_:'):
                    continue

                # Skip EC and SO parent hierarchy (not useful per team decision)
                if pred == 'rdfs:subClassOf':
                    if subj.startswith('EC:') or subj.startswith('SO:'):
                        continue

                # Track seed.reaction terms for backfill
                if str(obj).startswith('seed.reaction:'):
                    seed_reaction_terms.add(obj)

                clean_pred = predicate_labels.get(pred, pred)
                relationships.append({
                    'subject': subj,
                    'predicate': clean_pred,
                    'object': obj
                })

            print(f"   Found {len(relationships)} relationships")
            print(f"   Found {len(seed_reaction_terms)} seed.reaction terms")

            # Backfill seed.reaction terms into enriched_terms
            if seed_reaction_terms:
                print(f"   Backfilling seed.reaction term labels...")
                mask = (
                    statements_df['subject'].isin(seed_reaction_terms) &
                    statements_df['predicate'].isin(['rdfs:label', 'IAO:0000115'])
                )
                rxn_filtered = statements_df[mask]

                rxn_info = {}
                for _, row in rxn_filtered.iterrows():
                    subj = row['subject']
                    pred = row['predicate']
                    val = row['value'] if 'value' in row else ''

                    if subj not in rxn_info:
                        rxn_info[subj] = {'label': '', 'definition': ''}

                    if pred == 'rdfs:label':
                        rxn_info[subj]['label'] = val
                    elif pred == 'IAO:0000115':
                        rxn_info[subj]['definition'] = val

                for rxn_id in seed_reaction_terms:
                    info = rxn_info.get(rxn_id, {'label': '', 'definition': ''})
                    enriched_terms.append({
                        'ontology_prefix': 'seed.reaction',
                        'identifier': rxn_id,
                        'label': info['label'],
                        'definition': info['definition']
                    })

        # =====================================================================
        # STEP 7: Add EC column to ontology terms
        # =====================================================================
        print("\n7. Adding EC column to ontology terms...")

        # Load KEGG KO -> EC mapping from reference file
        kegg_ec_mapping_path = ref_path / "kegg_ko_ec_mapping.tsv"
        ko_to_ec = {}

        if kegg_ec_mapping_path.exists():
            print(f"   Loading KEGG KO -> EC mapping...")
            with open(kegg_ec_mapping_path, 'r') as f:
                for line in f:
                    line = line.strip()
                    if '\t' in line:
                        ec_raw, ko_raw = line.split('\t')
                        ec_id = ec_raw.replace('ec:', 'EC:')
                        ko_id = ko_raw.replace('ko:', 'KEGG:')
                        if ko_id not in ko_to_ec:
                            ko_to_ec[ko_id] = []
                        ko_to_ec[ko_id].append(ec_id)
            print(f"   Loaded {len(ko_to_ec)} KEGG KO -> EC mappings")
        else:
            print(f"   Warning: kegg_ko_ec_mapping.tsv not found at {kegg_ec_mapping_path}")

        # Patterns for extracting EC and TC from labels
        ec_label_pattern = re.compile(r'\(EC\s*([\d\.-]+)\)')
        tc_label_pattern = re.compile(r'\(TC\s*([\d\.\w]+)\)')

        kegg_ec_count = 0
        seed_ec_count = 0
        seed_tc_count = 0
        ec_copy_count = 0

        for term in enriched_terms:
            ec_values = []
            prefix = term['ontology_prefix']
            identifier = term['identifier']
            label = term.get('label', '')

            if prefix == 'KEGG':
                # KEGG KO: lookup from mapping file
                if identifier in ko_to_ec:
                    ec_values.extend(ko_to_ec[identifier])
                    kegg_ec_count += 1

            elif prefix == 'seed.role':
                # seed.role: extract EC and TC from label
                if label:
                    ec_matches = ec_label_pattern.findall(label)
                    if ec_matches:
                        ec_values.extend(['EC:' + m for m in ec_matches])
                        seed_ec_count += 1

                    tc_matches = tc_label_pattern.findall(label)
                    if tc_matches:
                        ec_values.extend(['TC:' + m for m in tc_matches])
                        seed_tc_count += 1

            elif prefix == 'EC':
                # EC terms: copy identifier itself
                ec_values.append(identifier)
                ec_copy_count += 1

            # Join multiple values with pipe
            term['ec'] = '|'.join(ec_values) if ec_values else ''

        print(f"   KEGG KO with EC: {kegg_ec_count}")
        print(f"   seed.role with EC: {seed_ec_count}")
        print(f"   seed.role with TC: {seed_tc_count}")
        print(f"   EC terms copied: {ec_copy_count}")

        total_with_ec = sum(1 for t in enriched_terms if t.get('ec'))
        print(f"   Total terms with ec column: {total_with_ec}")

        # =====================================================================
        # STEP 8: Create ontology definitions
        # =====================================================================
        print("\n8. Creating ontology definitions...")

        ontology_definitions = {
            'GO': 'Gene Ontology - standardized vocabulary for gene and protein functions',
            'EC': 'Enzyme Commission numbers - classification of enzymes by reaction type',
            'SO': 'Sequence Ontology - vocabulary for sequence features',
            'PFAM': 'Protein Families database - protein domain families',
            'KEGG': 'KEGG Orthologs - ortholog groups linking genes across species',
            'COG': 'Clusters of Orthologous Groups - protein functional categories',
            'seed.role': 'SEED Role Ontology - functional roles from RAST annotation',
            'seed.reaction': 'SEED Reaction Ontology - biochemical reactions from ModelSEED',
        }

        # Only include definitions for prefixes we actually have terms for
        present_prefixes = set(t['ontology_prefix'] for t in enriched_terms)
        definition_rows = [
            {'ontology_prefix': prefix, 'definition': desc}
            for prefix, desc in ontology_definitions.items()
            if prefix in present_prefixes
        ]

        # =====================================================================
        # STEP 9: Save output files
        # =====================================================================
        print(f"\n9. Saving output to {output_path}...")

        os.makedirs(output_path, exist_ok=True)

        # Save ontology_terms.tsv
        terms_df = pd.DataFrame(enriched_terms)
        terms_df = terms_df.drop_duplicates(subset=['identifier'])
        # Sort by ontology_prefix, then by identifier for proper ordering
        terms_df = terms_df.sort_values(['ontology_prefix', 'identifier']).reset_index(drop=True)
        terms_path = output_path / 'ontology_terms.tsv'
        terms_df.to_csv(terms_path, sep='\t', index=False)
        print(f"   Saved {len(terms_df)} terms to ontology_terms.tsv")

        # Summary by prefix
        for prefix in terms_df['ontology_prefix'].unique():
            count = len(terms_df[terms_df['ontology_prefix'] == prefix])
            print(f"     {prefix}: {count}")

        # Save ontology_definition.tsv
        defs_df = pd.DataFrame(definition_rows)
        defs_path = output_path / 'ontology_definition.tsv'
        defs_df.to_csv(defs_path, sep='\t', index=False)
        print(f"   Saved {len(defs_df)} definitions to ontology_definition.tsv")

        # Save ontology_relationships.tsv
        rels_df = pd.DataFrame(relationships)
        if not rels_df.empty:
            rels_df = rels_df.drop_duplicates()
        rels_path = output_path / 'ontology_relationships.tsv'
        rels_df.to_csv(rels_path, sep='\t', index=False)
        print(f"   Saved {len(rels_df)} relationships to ontology_relationships.tsv")

        if not rels_df.empty:
            print(f"   By predicate:")
            for pred in rels_df['predicate'].unique():
                count = len(rels_df[rels_df['predicate'] == pred])
                print(f"     {pred}: {count}")

        print(f"\n{'='*70}")
        print(f"Ontology table generation complete!")
        print(f"{'='*70}")

        return True

    except Exception as e:
        import traceback
        print(f"\nError generating ontology tables: {e}")
        print(traceback.format_exc())
        return False

class RASTSeedMapper:
    """
    Maps RAST annotations to SEED role ontology identifiers.
    
    This mapper handles:
    - Direct exact matches
    - Multi-function annotations with various separators (/, @, ;)
    - Different SEED ontology formats (URL-based and clean IDs)
    
    Usage:
        mapper = RASTSeedMapper("/data/reference_data/seed.json")
        seed_id = mapper.map_annotation("Alcohol dehydrogenase")
        # Returns: "seed.role:0000000001234"
    """
    
    def __init__(self, seed_ontology_path: str):
        """
        Initialize the mapper with a SEED ontology file.
        
        Args:
            seed_ontology_path: Path to SEED ontology JSON file (seed.json)
        """
        self.seed_mapping = {}
        self.multi_func_separators = [' / ', ' @ ', '; ']
        self._load_seed_ontology(seed_ontology_path)
    
    def _load_seed_ontology(self, path: str) -> None:
        """Load SEED ontology from JSON file."""
        import json
        from pathlib import Path
        
        path = Path(path)
        if not path.exists():
            raise FileNotFoundError(f"Ontology file not found: {path}")
        
        with open(path, 'r', encoding='utf-8') as f:
            data = json.load(f)
        
        # Extract nodes from JSON-LD format
        graphs = data.get("graphs", [])
        if not graphs:
            print("Warning: No graphs found in ontology file")
            return
            
        nodes = graphs[0].get("nodes", [])
        
        for node in nodes:
            label = node.get("lbl")
            node_id = node.get("id")
            
            if not label or not node_id:
                continue
                
            # Parse different ID formats
            seed_role_id = self._parse_seed_role_id(node_id)
            if seed_role_id:
                self.seed_mapping[label] = seed_role_id
                
        print(f"    Loaded {len(self.seed_mapping)} SEED role mappings")
        
    def _parse_seed_role_id(self, raw_id: str) -> str:
        """Parse SEED role ID from various formats."""
        if not raw_id:
            return None
            
        # URL format with Role parameter
        if "Role=" in raw_id:
            try:
                role_number = raw_id.split("Role=")[-1]
                return f"seed.role:{role_number}"
            except IndexError:
                return None
                
        # Already in clean format
        if raw_id.startswith("seed.role:"):
            return raw_id
            
        # OBO-style IDs (e.g., seed.role_0000000001234)
        if '_' in raw_id and 'seed.role_' in raw_id:
            ontology_part = raw_id.split('/')[-1]
            return ontology_part.replace("_", ":", 1)
            
        return None
    
    def split_multi_function(self, annotation: str) -> list:
        """Split multi-function annotations into individual components."""
        if not annotation:
            return []
            
        parts = [annotation]
        for separator in self.multi_func_separators:
            new_parts = []
            for part in parts:
                split_parts = part.split(separator)
                new_parts.extend(p.strip() for p in split_parts if p.strip())
            parts = new_parts
            
        return parts
    
    def map_annotation(self, annotation: str) -> str:
        """
        Map a RAST annotation to its SEED role ID.
        
        Args:
            annotation: RAST annotation string
            
        Returns:
            seed.role ID if found, None otherwise
        """
        if not annotation:
            return None
            
        # Try direct match first
        if annotation in self.seed_mapping:
            return self.seed_mapping[annotation]
            
        # Try splitting multi-function annotations
        parts = self.split_multi_function(annotation)
        
        if len(parts) > 1:
            for part in parts:
                if part in self.seed_mapping:
                    return self.seed_mapping[part]
                    
        return None
    
    def map_all_annotations(self, annotation: str) -> list:
        """
        Map a RAST annotation to ALL matching SEED role IDs.
        
        For multi-function annotations like "Thioredoxin / Glutaredoxin",
        returns all matching roles.
        
        Args:
            annotation: RAST annotation string
            
        Returns:
            List of tuples (matched_part, seed_role_id)
        """
        if not annotation:
            return []
        
        results = []
        
        # Try direct match first
        if annotation in self.seed_mapping:
            results.append((annotation, self.seed_mapping[annotation]))
        
        # Try splitting multi-function annotations
        parts = self.split_multi_function(annotation)
        
        for part in parts:
            if part in self.seed_mapping and part != annotation:
                results.append((part, self.seed_mapping[part]))
        
        return results
