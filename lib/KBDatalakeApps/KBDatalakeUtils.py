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

global_kbversion = "prod"

class KBDataLakeUtils(KBGenomeUtils, MSReconstructionUtils, MSFBAUtils):
    def __init__(self,kbversion,**kwargs):
        super().__init__(
                name="KBDataLakeUtils",
                kbversion=kbversion,
                **kwargs
        )
        global global_kbversion
        global_kbversion = kbversion

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

def run_phenotype_simulation(model_filename,output_filename,max_phenotypes=5):
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
    reference_path = "/data/reference_data"

    # Load model
    model = cobra.io.load_json_model(model_filename)
    mdlutl = MSModelUtil(model)

    #Loading the phenotype set from the reference path
    filename = reference_path + "/phenotypes/full_phenotype_set.json"
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
    pheno_util = PhenotypeWorkerUtil(kbversion=global_kbversion)

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

def run_model_reconstruction(input_filename,output_filename):
    worker_util = MSReconstructionUtils(kbversion=global_kbversion)

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
    genome_classifier = worker_util.get_classifier()

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
