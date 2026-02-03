import sys
import os
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


class KBDataLakeUtils(KBReadsUtils, KBGenomeUtils, SKANIUtils, MSReconstructionUtils, MSFBAUtils, MSBiochemUtils):
    def __init__(self, directory, worker_count, parameters, **kwargs):
        super().__init__(
            name="KBDataLakeUtils",
            **kwargs
        )
        self.directory = directory
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
        """
        #TODO: Jose should fill in his code for kberdle query using API here using SKANI results as input

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

                # Build rows
                gene_rows = []
                for ftr_id, ftr in self.ftrhash.items():
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
                        'type': self.ftrtypes.get(ftr_id, 'Unknown'),
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
        Submits protein sequences for RAST annotation and adds results to genome TSV files.
        """
        genomes_dir = os.path.join(self.directory, "genomes")
        genome_files = [f for f in os.listdir(genomes_dir) if f.endswith('.tsv')]

        rast_client = self.rast_client()

        for genome_file in genome_files:
            genome_id = genome_file.replace('.tsv', '')
            filepath = os.path.join(genomes_dir, genome_file)
            print(f"Annotating {genome_id} with RAST...")

            try:
                df = pd.read_csv(filepath, sep='\t')

                # Collect protein sequences for annotation
                proteins = []
                protein_indices = []
                for idx, row in df.iterrows():
                    protein = row.get('protein_translation', '')
                    if pd.notna(protein) and protein:
                        proteins.append(str(protein))
                        protein_indices.append(idx)

                if not proteins:
                    df['rast_functions'] = ''
                    df.to_csv(filepath, sep='\t', index=False)
                    continue

                # Call RAST annotate_proteins
                result = rast_client.annotate_proteins({'proteins': proteins})
                functions_list = result.get('functions', [])

                # Map RAST annotations back to DataFrame
                rast_col = [''] * len(df)
                for i, idx in enumerate(protein_indices):
                    if i < len(functions_list):
                        funcs = functions_list[i]
                        if isinstance(funcs, list):
                            rast_col[idx] = ';'.join(funcs)
                        elif isinstance(funcs, str):
                            rast_col[idx] = funcs

                df['rast_functions'] = rast_col
                df.to_csv(filepath, sep='\t', index=False)
                print(f"  Added RAST annotations for {len(protein_indices)} proteins in {genome_id}")

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
        """
        #TODO: Jose should fill in his code for ontology term query using API here using genome TSV files as input

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
                'model_id': model_id,
                'model_path': os.path.join(models_dir, model_file),
                'phenotypes_dir': phenotypes_dir,
            })

        errors = []
        completed = 0

        with ProcessPoolExecutor(max_workers=self.worker_count) as executor:
            future_to_model = {
                executor.submit(_simulate_phenotypes_worker, item): item['model_id']
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
                # Add RAST annotations if present
                rast_funcs = gene.get('rast_functions', '')
                if pd.notna(rast_funcs) and rast_funcs:
                    for func in str(rast_funcs).split(';'):
                        if func.strip():
                            feature.add_ontology_term('RAST', func.strip())
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

    sys.path = [
        "/deps/KBUtilLib/src",
        "/deps/cobrakbase",
        "/deps/ModelSEEDpy",
    ] + sys.path

    try:
        import cobra
        from modelseedpy import MSModelUtil

        model_id = work_item['model_id']
        model_path = work_item['model_path']
        phenotypes_dir = work_item['phenotypes_dir']

        # Load model
        model = cobra.io.load_json_model(model_path)
        MSModelUtil(model)

        # Build result structure
        result = {
            'model_id': model_id,
            'success': True,
            'accuracy': {},
            'gene_phenotype_reactions': [],
            'phenotype_gaps': [],
            'gapfilled_reactions': [],
            'num_phenotypes': 0,
        }

        # Extract gapfilled reactions info
        for rxn in model.reactions:
            if hasattr(rxn, 'annotation'):
                sbo = rxn.annotation.get('sbo', '')
                if sbo == 'SBO:0000176':
                    result['gapfilled_reactions'].append({
                        'reaction_id': rxn.id,
                        'reaction_name': rxn.name,
                        'gpr': rxn.gene_reaction_rule,
                        'lower_bound': rxn.lower_bound,
                        'upper_bound': rxn.upper_bound,
                    })

        # Save results
        output_path = os.path.join(phenotypes_dir, f'{model_id}_phenosim.json')
        with open(output_path, 'w') as f:
            json.dump(result, f, indent=2)

        return result

    except Exception as e:
        import traceback
        return {
            'success': False,
            'model_id': work_item.get('model_id', 'unknown'),
            'error': f"{str(e)}\n{traceback.format_exc()}"
        }
