import sys
import os
import re
import sqlite3
import json
import uuid
import shutil
from os import path
from concurrent.futures import ProcessPoolExecutor, as_completed
# Add the parent directory to the sys.path
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
script_path = os.path.abspath(__file__)
script_dir = os.path.dirname(script_path)
base_dir = os.path.dirname(os.path.dirname(script_dir))
folder_name = os.path.basename(script_dir)

sys.path = ["/deps/KBUtilLib/src","/deps/cobrakbase","/deps/ModelSEEDpy"] + sys.path

# Import utilities with error handling
from kbutillib import KBModelUtils, KBReadsUtils, KBGenomeUtils, SKANIUtils, MSReconstructionUtils, MSFBAUtils, MSBiochemUtils

import hashlib
import pandas as pd
import cobra
from modelseedpy import AnnotationOntology, MSPackageManager, MSTemplateBuilder, MSMedia, MSModelUtil, MSBuilder, MSATPCorrection, MSGapfill, MSGrowthPhenotype, MSGrowthPhenotypes, ModelSEEDBiochem, MSExpression
from modelseedpy.helpers import get_template

# BERDL query modules (Jose P. Faria)
try:
    from berdl.berdl import QueryPangenomeBERDL, OntologyEnrichment
except ImportError:
    # Fallback for when berdl module is not installed
    QueryPangenomeBERDL = None
    OntologyEnrichment = None


class KBDataLakeUtils(KBReadsUtils, KBGenomeUtils, SKANIUtils, MSReconstructionUtils, MSFBAUtils, MSBiochemUtils):
    def __init__(self, directory, reference_path, worker_count, parameters, **kwargs):
        super().__init__(
                name="KBDataLakeUtils",
                **kwargs
        )
        self.directory = directory
        self.reference_path = reference_path
        self.worker_count = worker_count
        self.app_parameters = parameters
        self.workspace_name = self.app_parameters['workspace_name']
        os.makedirs(self.directory , exist_ok=True)

    def run_full_pipeline(self):
        """
        Run the full pipeline for modeling analysis.
        """
        # Step 1: Process input arguments into user genome table
        self.pipeline_process_arguments_into_user_genome_table()

        # Step 2: Download assemblies and genome genes (independent, can run in parallel via threads)
        self.pipeline_download_user_genome_assmemblies()
        self.pipeline_download_user_genome_genes()

        # Step 3: Run SKANI analysis on downloaded assemblies
        self.pipeline_run_skani_analysis()

        # Step 4: Annotate genomes with RAST
        self.pipeline_annotate_user_genome_with_rast()

        # Step 5: Build metabolic models (parallelized internally via ProcessPoolExecutor)
        self.pipeline_run_moddeling_analysis()

        # Step 6: Run phenotype simulations (parallelized internally via ProcessPoolExecutor)
        self.pipeline_run_phenotype_simulations()

        # Step 7: Build SQLite database from all output data
        self.pipeline_build_sqllite_db()

        # Step 8: Save outputs to KBase workspace
        self.pipeline_save_annotated_genomes()
        self.pipeline_save_models_to_kbase()
        self.pipeline_save_genometables_workspace_object()

        # Step 9: Generate and save report
        self.pipeline_save_kbase_report()

    def pipeline_process_arguments_into_user_genome_table(self):
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
        os.makedirs(self.directory, exist_ok=True)
        output_path = os.path.join(self.directory, "user_genomes.tsv")
        df.to_csv(output_path, sep='\t', index=False)
        print(f"Saved {len(df)} genomes to {output_path}")

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

    def pipeline_download_user_genome_assmemblies(self):
        """
        Pipeline step for downloading genome assemblies.
        Reads user_genomes.tsv and downloads one assembly per genome to
        self.directory/assemblies/<genome_id>.fasta.  If multiple genomes
        share the same assembly_ref the file is copied rather than
        re-downloaded.
        """
        genomes_file = os.path.join(self.directory, "user_genomes.tsv")
        df = pd.read_csv(genomes_file, sep='\t')

        assemblies_dir = os.path.join(self.directory, "assemblies")
        os.makedirs(assemblies_dir, exist_ok=True)

        downloaded = {}  # assembly_ref -> first local fasta path

        for _, row in df.iterrows():
            genome_id = row['genome_id']
            assembly_ref = row.get('assembly_ref', '')
            if not assembly_ref or pd.isna(assembly_ref):
                print(f"Warning: No assembly_ref for {genome_id}, skipping")
                continue

            target_path = os.path.join(assemblies_dir, f"{genome_id}.fasta")

            if assembly_ref in downloaded:
                # Same assembly already downloaded for another genome – copy it
                shutil.copy2(downloaded[assembly_ref], target_path)
                print(f"  Copied assembly for {genome_id} (shared assembly_ref)")
                continue

            print(f"Downloading assembly for {genome_id}...")
            try:
                assembly_set = self.download_assembly([assembly_ref], assemblies_dir)
                # Rename the assembly-id-named file to genome_id.fasta
                for _name, assembly in assembly_set.assemblies.items():
                    src = assembly.fasta_file
                    if os.path.exists(src):
                        shutil.move(src, target_path)
                    break
                downloaded[assembly_ref] = target_path
                print(f"  Saved: {target_path}")
            except Exception as e:
                print(f"  Warning: Failed to download assembly for {genome_id}: {e}")

        # Clean up the metadata file that download_assembly creates
        meta_file = os.path.join(assemblies_dir, "assemblies_metadata.json")
        if os.path.exists(meta_file):
            os.remove(meta_file)

        print(f"Downloaded {len(downloaded)} unique assemblies for {len(df)} genomes")

    def pipeline_run_skani_analysis(self):
        """
        Pipeline step for running SKANI analysis.
        Runs SKANI against pangenome, fitness, and phenotype sketch databases.
        Results saved to self.directory/skani as TSV files.
        """
        assemblies_dir = os.path.join(self.directory, "assemblies")
        skani_dir = os.path.join(self.directory, "skani")
        os.makedirs(skani_dir, exist_ok=True)

        # Collect all FASTA files from assemblies directory
        fasta_files = []
        for f in os.listdir(assemblies_dir):
            if f.endswith(('.fasta', '.fa', '.fna')):
                fasta_files.append(os.path.join(assemblies_dir, f))

        if not fasta_files:
            print("No FASTA files found in assemblies directory")
            return

        # Run SKANI against each sketch database
        databases = ['pangenome', 'fitness', 'phenotype']
        for db_name in databases:
            print(f"Running SKANI against {db_name} database...")
            try:
                results = self.query_genomes(
                    query_fasta=fasta_files,
                    database_name=db_name,
                    threads=self.worker_count
                )

                # Build output table: genome_id, reference_genome, ani_percentage
                rows = []
                for query_id, hits in results.items():
                    for hit in hits:
                        rows.append({
                            'genome_id': query_id,
                            'reference_genome': hit.get('reference', ''),
                            'ani_percentage': hit.get('ani', 0.0)
                        })

                output_df = pd.DataFrame(rows)
                output_path = os.path.join(skani_dir, f"skani_{db_name}.tsv")
                output_df.to_csv(output_path, sep='\t', index=False)
                print(f"  Saved {len(rows)} hits to {output_path}")

            except Exception as e:
                print(f"  Warning: SKANI analysis against {db_name} failed: {e}")

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
        if QueryPangenomeBERDL is None:
            print("Warning: BERDL pangenome module not available, skipping")
            return

        # Get BERDL service account token (Boris's forever token)
        token = self.get_token(namespace="berdl")
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

    def pipeline_download_user_genome_genes(self):
        """
        Pipeline step for downloading genome genes and annotations.
        Creates a TSV table per genome in self.directory/genomes/<genome_id>.tsv.

        Uses KBAnnotationUtils.process_object to standardise features and
        ontology terms.  Each distinct ontology type discovered on a genome's
        features becomes its own column headed ``Annotation:<type>`` (e.g.
        ``Annotation:SSO``, ``Annotation:KO``).  Individual terms within the
        column are semicolon-separated and formatted as::

            TERMID:description|rxn1,rxn2

        where the reaction IDs are ModelSEED reaction mappings (omitted when
        no mapping exists).
        """
        genomes_file = os.path.join(self.directory, "user_genomes.tsv")
        df = pd.read_csv(genomes_file, sep='\t')

        genomes_dir = os.path.join(self.directory, "genomes")
        os.makedirs(genomes_dir, exist_ok=True)

        for _, row in df.iterrows():
            genome_id = row['genome_id']
            genome_ref = row['genome_ref']

            print(f"Downloading genes for {genome_id}...")
            try:
                # Load genome object into object_hash
                obj = self.load_kbase_gene_container(genome_ref, localname=genome_id)
                genome_data = obj["data"]
                obj_type = (
                    obj.get("info", [None] * 3)[2]
                    if "info" in obj
                    else "KBaseGenomes.Genome"
                )

                # Use KBAnnotationUtils.process_object to standardise
                # features, aliases, and ontology terms in self.ftrhash
                self.process_object({"object": genome_data, "type": obj_type})

                # Discover all cleaned ontology types across this genome
                all_ontology_types = set()
                for ftr in self.ftrhash.values():
                    for raw_tag in ftr.get("ontology_terms", {}):
                        all_ontology_types.add(self.clean_tag(raw_tag))
                sorted_ont_types = sorted(all_ontology_types)

                # Build rows – only genes and noncoding features
                gene_rows = []
                for ftr_id, ftr in self.ftrhash.items():
                    ftr_type = self.ftrtypes.get(ftr_id, 'Unknown')
                    if ftr_type not in ('gene', 'noncoding'):
                        continue
                    # Aliases (upgrade_feature already normalised to [[src, val], ...])
                    aliases = ftr.get("aliases", [])
                    alias_str = (
                        ";".join(f"{a[0]}:{a[1]}" for a in aliases if len(a) >= 2)
                        if aliases else ""
                    )

                    # Location
                    locations = ftr.get("location", [])
                    contig = locations[0][0] if locations else ""
                    start = locations[0][1] if locations else 0
                    strand = locations[0][2] if locations else "+"
                    length = locations[0][3] if locations else 0
                    end = start + length if strand == "+" else start - length

                    # Functions (upgrade_feature converted 'function' → 'functions' list)
                    functions = ftr.get("functions", [])
                    if isinstance(functions, list):
                        functions_str = ";".join(functions)
                    else:
                        functions_str = str(functions) if functions else ""

                    row_data = {
                        'gene_id': ftr.get('id', ''),
                        'aliases': alias_str,
                        'contig': contig,
                        'start': start,
                        'end': end,
                        'strand': strand,
                        'type': ftr_type,
                        'functions': functions_str,
                        'protein_translation': ftr.get('protein_translation', ''),
                        'dna_sequence': ftr.get('dna_sequence', ''),
                    }

                    # Add a column per ontology type
                    for ont_type in sorted_ont_types:
                        terms_list = []
                        for raw_tag, terms_dict in ftr.get("ontology_terms", {}).items():
                            if self.clean_tag(raw_tag) != ont_type:
                                continue
                            for raw_term in terms_dict:
                                cleaned = self.clean_term(raw_term, raw_tag, ont_type)
                                name = self.get_term_name(ont_type, cleaned)
                                rxns = self.translate_term_to_modelseed(cleaned)
                                entry = f"{cleaned}:{name}"
                                if rxns:
                                    entry += "|" + ",".join(rxns)
                                terms_list.append(entry)
                        row_data[f"Annotation:{ont_type}"] = (
                            ";".join(terms_list) if terms_list else ""
                        )

                    gene_rows.append(row_data)

                gene_df = pd.DataFrame(gene_rows)
                output_path = os.path.join(genomes_dir, f"{genome_id}.tsv")
                gene_df.to_csv(output_path, sep='\t', index=False)
                ont_msg = ", ".join(sorted_ont_types) if sorted_ont_types else "none"
                print(f"  Saved {len(gene_rows)} features for {genome_id} "
                      f"(ontologies: {ont_msg})")

            except Exception as e:
                print(f"  Warning: Failed to download genes for {genome_id}: {e}")

    def pipeline_annotate_user_genome_with_rast(self):
        """
        Pipeline step for annotating genomes with RAST.
        Submits protein sequences for RAST annotation, translates the
        returned function strings to SSO terms, and populates the
        ``Annotation:SSO`` column in each genome TSV file.

        Skips any genome whose TSV already contains a non-empty
        ``Annotation:SSO`` column (e.g. from existing KBase annotations).
        """
        genomes_dir = os.path.join(self.directory, "genomes")
        genome_files = [f for f in os.listdir(genomes_dir) if f.endswith('.tsv')]

        rast_client = self.rast_client()

        for genome_file in genome_files:
            genome_id = genome_file.replace('.tsv', '')
            filepath = os.path.join(genomes_dir, genome_file)

            try:
                df = pd.read_csv(filepath, sep='\t')

                # Skip if Annotation:SSO already has data
                if 'Annotation:SSO' in df.columns:
                    non_empty = (
                        df['Annotation:SSO']
                        .fillna('')
                        .astype(str)
                        .str.strip()
                        .ne('')
                        .sum()
                    )
                    if non_empty > 0:
                        print(f"  Skipping {genome_id}: Annotation:SSO already "
                              f"populated ({non_empty} entries)")
                        continue

                print(f"Annotating {genome_id} with RAST...")

                # Collect protein sequences for annotation
                proteins = []
                protein_indices = []
                for idx, row in df.iterrows():
                    protein = row.get('protein_translation', '')
                    if pd.notna(protein) and protein:
                        proteins.append(str(protein))
                        protein_indices.append(idx)

                if not proteins:
                    if 'Annotation:SSO' not in df.columns:
                        df['Annotation:SSO'] = ''
                    df.to_csv(filepath, sep='\t', index=False)
                    continue

                # Call RAST annotate_proteins
                result = rast_client.annotate_proteins({'proteins': proteins})
                functions_list = result.get('functions', [])

                # Translate RAST functions → SSO terms and populate column
                sso_col = [''] * len(df)
                for i, idx in enumerate(protein_indices):
                    if i >= len(functions_list):
                        continue
                    funcs = functions_list[i]
                    if isinstance(funcs, str):
                        funcs = [funcs]
                    if not isinstance(funcs, list):
                        continue

                    sso_entries = []
                    for func_str in funcs:
                        # Split multi-role function strings the same way
                        # KBAnnotationUtils.upgrade_feature does
                        roles = re.split(r"\s*;\s+|\s+[\@\/]\s+", func_str)
                        for role in roles:
                            role = role.strip()
                            if not role:
                                continue
                            sso_id = self.translate_rast_function_to_sso(role)
                            if sso_id is None:
                                continue
                            name = self.get_term_name("SSO", sso_id)
                            rxns = self.translate_term_to_modelseed(sso_id)
                            entry = f"{sso_id}:{name}"
                            if rxns:
                                entry += "|" + ",".join(rxns)
                            sso_entries.append(entry)
                    sso_col[idx] = ";".join(sso_entries)

                df['Annotation:SSO'] = sso_col
                df.to_csv(filepath, sep='\t', index=False)
                annotated = sum(1 for x in sso_col if x)
                print(f"  Added RAST/SSO annotations for {annotated}/"
                      f"{len(protein_indices)} proteins in {genome_id}")

            except Exception as e:
                print(f"  Warning: RAST annotation failed for {genome_id}: {e}")

    def pipeline_run_moddeling_analysis(self):
        """
        Pipeline step for building metabolic models for a list of genomes.
        Runs in parallel using ProcessPoolExecutor.
        """
        models_dir = os.path.join(self.directory, "models")
        os.makedirs(models_dir, exist_ok=True)

        genome_dir = os.path.join(self.directory, "genomes")
        genome_files = [f for f in os.listdir(genome_dir) if f.endswith('.tsv')]
        genome_ids = [f.replace('.tsv', '') for f in genome_files]

        print(f"\nBuilding {len(genome_ids)} models with {self.worker_count} workers")

        # Prepare work items with all data needed by the worker
        work_items = []
        for genome_id in genome_ids:
            genome_tsv = os.path.join(genome_dir, f"{genome_id}.tsv")
            work_items.append({
                'genome_id': genome_id,
                'genome_tsv': genome_tsv,
                'models_dir': models_dir,
                'gapfill_media_ref': 'KBaseMedia/Carbon-Pyruvic-Acid'
            })

        errors = []
        completed = 0

        # Use ProcessPoolExecutor for true parallelism
        with ProcessPoolExecutor(max_workers=self.worker_count) as executor:
            future_to_genome = {
                executor.submit(_build_single_model_worker, item): item['genome_id']
                for item in work_items
            }

            for future in as_completed(future_to_genome):
                genome_id = future_to_genome[future]
                completed += 1
                try:
                    result = future.result()
                    if result.get('success'):
                        info = result.get('model_info', {})
                        print(f"[{completed}/{len(genome_ids)}] {genome_id}: "
                              f"{info.get('num_reactions', 0)} rxns, "
                              f"{info.get('num_genes', 0)} genes, "
                              f"class={info.get('genome_class', 'N/A')}")
                    else:
                        errors.append((genome_id, result.get('error', 'Unknown error')))
                        print(f"[{completed}/{len(genome_ids)}] {genome_id}: FAILED - {result.get('error', 'Unknown')}")
                except Exception as e:
                    errors.append((genome_id, str(e)))
                    print(f"[{completed}/{len(genome_ids)}] {genome_id}: EXCEPTION - {str(e)}")

        if errors:
            print(f"\n{len(errors)} models failed to build:")
            for gid, err in errors:
                print(f"  {gid}: {err}")

    def pipeline_run_ontology_term_kberdl_query(self):
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

        # Get BERDL service account token (Boris's forever token)
        token = self.get_token(namespace="berdl")
        if not token:
            print("Warning: No BERDL token available, skipping ontology enrichment")
            return

        # Try to load genome data from multiple sources
        genome_dataframes = []
        data_source = None

        # Source 1: SQLite database
        db_path = os.path.join(self.directory, "berdl_tables.db")
        if os.path.exists(db_path):
            try:
                conn = sqlite3.connect(db_path)
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
            if os.path.exists(db_path):
                try:
                    conn = sqlite3.connect(db_path)
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
        genomes_file = os.path.join(self.directory, "user_genomes.tsv")
        genomes_dir = os.path.join(self.directory, "genomes")
        df = pd.read_csv(genomes_file, sep='\t')

        suffix = self.app_parameters.get('suffix', '.datalake')
        saved_refs = []

        for _, row in df.iterrows():
            genome_id = row['genome_id']
            genome_ref = row['genome_ref']
            genome_tsv = os.path.join(genomes_dir, f"{genome_id}.tsv")

            if not os.path.exists(genome_tsv):
                print(f"Warning: No genome TSV found for {genome_id}, skipping save")
                continue

            try:
                gene_df = pd.read_csv(genome_tsv, sep='\t')

                # Build annotations dict from RAST column
                # Format: {gene_id: {ontology: {term: {"type": "RAST"}}}}
                annotations = {}
                for _, gene_row in gene_df.iterrows():
                    gene_id = gene_row.get('gene_id', '')
                    rast_funcs = gene_row.get('rast_functions', '')
                    if pd.notna(rast_funcs) and rast_funcs:
                        annotations[gene_id] = {}
                        annotations[gene_id]['SSO'] = {}
                        for func in str(rast_funcs).split(';'):
                            func = func.strip()
                            if func:
                                annotations[gene_id]['SSO'][func] = {"type": "RAST"}

                if annotations:
                    # Load object info hash for add_annotations_to_object
                    self.object_to_proteins(genome_ref)
                    result = self.add_annotations_to_object(genome_ref, suffix, annotations)
                    saved_ref = result.get('output_ref', '')
                    if saved_ref:
                        saved_refs.append(saved_ref)
                    print(f"Saved annotated genome for {genome_id}: {saved_ref}")

            except Exception as e:
                print(f"Warning: Failed to save annotated genome {genome_id}: {e}")

        # Create GenomeSet with all saved genomes
        genome_set_name = self.app_parameters.get('genome_set_name', 'datalake_genomes')
        if saved_refs:
            genome_set_data = {
                'description': f'Genome set from datalake pipeline with {len(saved_refs)} genomes',
                'elements': {
                    f"genome_{i}": {"ref": ref}
                    for i, ref in enumerate(saved_refs)
                }
            }
            self.set_ws(self.workspace_name)
            params = {
                "id": self.ws_id,
                "objects": [{
                    "data": genome_set_data,
                    "name": genome_set_name,
                    "type": "KBaseSearch.GenomeSet",
                    "meta": {},
                    "provenance": self.provenance(),
                }]
            }
            self.ws_client().save_objects(params)
            print(f"Saved GenomeSet '{genome_set_name}' with {len(saved_refs)} genomes")

    def pipeline_save_models_to_kbase(self):
        """
        Pipeline step for saving metabolic models to the KBase workspace.
        """
        genomes_file = os.path.join(self.directory, "user_genomes.tsv")
        models_dir = os.path.join(self.directory, "models")
        df = pd.read_csv(genomes_file, sep='\t')

        genome_ids = set(df['genome_id'].tolist())

        for model_file in os.listdir(models_dir):
            if not model_file.endswith('_model.json'):
                continue

            model_id = model_file.replace('_model.json', '')
            # Convert safe ID back to check against genome_ids
            original_id = model_id.replace('_', '.')
            if model_id not in genome_ids and original_id not in genome_ids:
                # Also check with underscore-replaced version
                found = False
                for gid in genome_ids:
                    if gid.replace('.', '_') == model_id:
                        found = True
                        break
                if not found:
                    continue

            model_path = os.path.join(models_dir, model_file)
            try:
                model = cobra.io.load_json_model(model_path)
                mdlutl = MSModelUtil(model)
                mdlutl.wsid = model_id
                self.save_model(mdlutl, workspace=self.workspace_name, objid=model_id)
                print(f"Saved model {model_id} to workspace {self.workspace_name}")
            except Exception as e:
                print(f"Warning: Failed to save model {model_id}: {e}")

    def pipeline_run_phenotype_simulations(self):
        """
        Pipeline step for running phenotype simulations.
        Runs simulations in parallel and saves results to self.directory/phenotypes.
        """
        models_dir = os.path.join(self.directory, "models")
        phenotypes_dir = os.path.join(self.directory, "phenotypes")
        os.makedirs(phenotypes_dir, exist_ok=True)

        model_files = [f for f in os.listdir(models_dir) if f.endswith('_model.json')]
        if not model_files:
            print("No models found for phenotype simulation")
            return

        print(f"\nRunning phenotype simulations for {len(model_files)} models "
              f"with {self.worker_count} workers")

        # Prepare work items
        work_items = []
        for model_file in model_files:
            model_id = model_file.replace('_model.json', '')
            work_items.append({
                'genome_id': model_id,
                'directory': self.directory,
                'max_phenotypes': 5,
                'reference_path': self.reference_path
            })

        errors = []
        completed = 0

        with ProcessPoolExecutor(max_workers=self.worker_count) as executor:
            future_to_model = {
                executor.submit(_simulate_phenotypes_worker, item): item['genome_id']
                for item in work_items
            }

            for future in as_completed(future_to_model):
                model_id = future_to_model[future]
                completed += 1
                try:
                    result = future.result()
                    if result.get('success'):
                        print(f"[{completed}/{len(model_files)}] {model_id}: "
                              f"simulated {result.get('num_phenotypes', 0)} phenotypes")
                    else:
                        errors.append((model_id, result.get('error', 'Unknown')))
                        print(f"[{completed}/{len(model_files)}] {model_id}: FAILED - {result.get('error')}")
                except Exception as e:
                    errors.append((model_id, str(e)))
                    print(f"[{completed}/{len(model_files)}] {model_id}: EXCEPTION - {str(e)}")

        # Build summary tables from individual phenotype results
        self._build_phenotype_tables(phenotypes_dir)

    def _build_phenotype_tables(self, phenotypes_dir):
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


def _build_single_model_worker(work_item):
    """
    Worker function for parallel model building.
    This runs in a separate process with its own memory space.

    Args:
        work_item: Dictionary with genome_id, genome_tsv, models_dir, gapfill_media_ref

    Returns:
        Dictionary with success status, model_info or error message
    """
    import os
    import sys
    import pandas as pd

    # Re-setup paths for worker process
    sys.path = [
        "/deps/KBUtilLib/src",
        "/deps/cobrakbase",
        "/deps/ModelSEEDpy",
    ] + sys.path

    try:
        import cobra
        from modelseedpy.core.msgenome import MSGenome, MSFeature
        from modelseedpy.core.msmodelutl import MSModelUtil
        from kbutillib import MSReconstructionUtils

        # Create a minimal util instance for this worker
        class WorkerUtil(MSReconstructionUtils):
            def __init__(self):
                super().__init__(name="WorkerUtil")

        worker_util = WorkerUtil()

        genome_id = work_item['genome_id']
        genome_tsv = work_item['genome_tsv']
        models_dir = work_item['models_dir']
        gapfill_media_ref = work_item.get('gapfill_media_ref')

        # Clear MSModelUtil cache for this process
        MSModelUtil.mdlutls.clear()

        # Create safe model ID
        safe_genome_id = genome_id.replace('.', '_')
        model_path = os.path.join(models_dir, f'{safe_genome_id}_model.json')

        # Load features from genome TSV
        gene_df = pd.read_csv(genome_tsv, sep='\t')

        # Create MSGenome from features
        genome = MSGenome()
        genome.id = safe_genome_id
        genome.scientific_name = genome_id

        ms_features = []
        for _, gene in gene_df.iterrows():
            protein = gene.get('protein_translation', '')
            gene_id = gene.get('gene_id', '')
            if pd.notna(protein) and protein:
                feature = MSFeature(gene_id, str(protein))
                # Parse Annotation:SSO column
                # Format: SSO:nnnnn:description|rxn1,rxn2;SSO:mmmmm:desc2|rxn3
                sso_col = gene.get('Annotation:SSO', '')
                if pd.notna(sso_col) and sso_col:
                    for entry in str(sso_col).split(';'):
                        entry = entry.strip()
                        if not entry:
                            continue
                        term_part = entry.split('|')[0]
                        parts = term_part.split(':')
                        if len(parts) >= 2 and parts[0] == 'SSO':
                            sso_id = parts[0] + ':' + parts[1]
                            feature.add_ontology_term('SSO', sso_id)
                            # Extract description for classifier
                            if len(parts) >= 3:
                                description = ':'.join(parts[2:])
                                if description:
                                    feature.add_ontology_term('RAST', description)
                ms_features.append(feature)

        genome.add_features(ms_features)

        # Load classifier
        genome_classifier = worker_util.get_classifier()

        # Build the model
        current_output, mdlutl = worker_util.build_metabolic_model(
            genome=genome,
            genome_classifier=genome_classifier,
            model_id=safe_genome_id,
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
        if gapfill_media_ref:
            gapfill_media = worker_util.get_media(gapfill_media_ref)
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
        cobra.io.save_json_model(model, model_path)

        genome_class = current_output.get('Class', 'Unknown')
        core_gf = current_output.get('Core GF', 0)

        return {
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
            }
        }

    except Exception as e:
        import traceback
        return {
            'success': False,
            'error': f"{str(e)}\n{traceback.format_exc()}"
        }

def _simulate_phenotypes_worker(work_item):
    """
    Worker function for parallel phenotype simulation.
    This runs in a separate process with its own memory space.

    Args:
        work_item: Dictionary with model_id, model_path, phenotypes_dir

    Returns:
        Dictionary with success status and simulation results
    """
    import os
    import sys
    import json
    import cobra
    from modelseedpy import MSGrowthPhenotypes, MSATPCorrection,MSGapfill,MSModelUtil

    sys.path = [
        "/deps/KBUtilLib/src",
        "/deps/cobrakbase",
        "/deps/ModelSEEDpy",
    ] + sys.path

    try:
        genome_id = work_item['genome_id']
        phenotypes_dir = work_item['directory'] + '/phenotypes/'
        reference_path = work_item['reference_path']

        # Load model
        model = cobra.io.load_json_model(work_item['directory'] + '/models/' + genome_id + '_model.json')
        mdlutl = MSModelUtil(model)

        #Loading the phenotype set from the reference path
        filename = reference_path + "/phenotypes/full_phenotype_set.json"
        with open(filename) as f:
            phenoset_data = json.load(f)
        #Setting max phenotypes if specified in the work item
        if "max_phenotypes" in work_item:
            phenoset_data["phenotypes"] = phenoset_data["phenotypes"][:work_item["max_phenotypes"]]
        #Instantiating the phenotype set
        phenoset = MSGrowthPhenotypes.from_dict(phenoset_data)

        # Create a minimal util instance for this worker
        class PhenotypeWorkerUtil(MSReconstructionUtils,MSFBAUtils,MSBiochemUtils):
            def __init__(self):
                super().__init__(name="PhenotypeWorkerUtil")
        pheno_util = PhenotypeWorkerUtil()

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

        os.makedirs(phenotypes_dir, exist_ok=True)
        with open(phenotypes_dir+"/"+genome_id+".json", "w") as f:
            json.dump(results, f, indent=4, skipkeys=True)

        return {"success": True, "genome_id": genome_id}

    except Exception as e:
        import traceback
        return {
            'success': False,
            'model_id': work_item.get('model_id', 'unknown'),
            'error': f"{str(e)}\n{traceback.format_exc()}"
        }


# =============================================================================
# ONTOLOGY TABLE GENERATION FUNCTION
# =============================================================================
# Standalone function for generating ontology tables from a clade folder.
# Called by the pipeline coordinator to generate ontology data for each clade.
#
# Author: Jose P. Faria (jplfaria@gmail.com)
# Date: February 2026
# =============================================================================



# =============================================================================
# ONTOLOGY TABLE GENERATION FUNCTION
# =============================================================================
# Standalone function for generating ontology tables from a clade folder.
# Called by the pipeline coordinator to generate ontology data for each clade.
#
# This function:
# 1. Reads genome features from db.sqlite
# 2. Maps RAST functions to seed.role IDs using RASTSeedMapper
# 3. Extracts EC numbers from RAST function strings
# 4. Extracts existing ontology terms (GO, KEGG, COG, PFAM, SO)
# 5. Gets seed.role → seed.reaction relationships
# 6. Enriches all terms with labels and definitions
# 7. Saves ontology_terms.tsv, ontology_definition.tsv, ontology_relationships.tsv
#
# Author: Jose P. Faria (jplfaria@gmail.com)
# Date: February 2026
# =============================================================================

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
