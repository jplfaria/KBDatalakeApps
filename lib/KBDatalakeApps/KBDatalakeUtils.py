import sys
import os
import re
import io
import sqlite3
import json
import uuid
import shutil
import hashlib
from os import path
from contextlib import redirect_stdout

sys.path = ["/deps/KBUtilLib/src","/deps/cobrakbase","/deps/ModelSEEDpy","/deps/cb_annotation_ontology_api"] + sys.path

# Import utilities with error handling
from kbutillib import KBAnnotationUtils, KBReadsUtils, KBGenomeUtils, MSReconstructionUtils, MSFBAUtils, MSBiochemUtils

import pandas as pd
import cobra
from modelseedpy.core.msgenome import MSGenome, MSFeature
from modelseedpy import AnnotationOntology, MSPackageManager, MSMedia, MSModelUtil, MSBuilder, MSATPCorrection, MSGapfill, MSGrowthPhenotype, MSGrowthPhenotypes, ModelSEEDBiochem, MSExpression
from modelseedpy.core.mstemplate import MSTemplateBuilder
# from KBDatalakeApps.ontology_enrichment import OntologyEnrichment

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

    def pipeline_save_annotated_genomes(self):
        """
        Pipeline step for saving annotated genomes back to KBase.
        Uses annotation ontology API to save RAST annotations to genome objects.
        """
        pass

    def build_phenotype_tables(self, output_dir, phenosim_directory,
                               experiment_data_file=None, phenoset_file=None,
                               fitness_mapping_dir=None, fitness_genomes_dir=None,
                               model_data_dir=None, reference_phenosim_dir=None,
                               fitness_threshold=-0.5, essential_genes_file=None):
        """
        Pipeline step for building phenotype tables.

        Reads individual per-genome phenosim JSON files from phenosim_directory.
        Each file is named {genome_id}.json and contains:
          - details: parallel arrays (Phenotype, Class, Simulated objective,
                     Observed objective, Transports missing)
          - data: dict keyed by cpd_id with objective_value, class, reactions,
                  gfreactions, gapfill_count, fluxes, etc.
          - data.summary: dict with accuracy, CP, CN, FP, FN, P, N

        Creates three TSV files:
        1. model_performance.tsv - Model accuracy metrics per genome
        2. genome_phenotypes.tsv - Phenotype predictions per genome (with transports)
        3. gene_phenotypes.tsv - Gene-phenotype associations from gapfilling,
           model predictions, and fitness data

        Args:
            output_dir: Directory to save output TSV files
            phenosim_directory: Directory containing per-genome phenosim JSON files
            experiment_data_file: Path to experimental_data.json (optional)
            phenoset_file: Path to full_phenotype_set.json (optional)
            fitness_mapping_dir: Directory containing {genome_id}_fitness.parquet
                                 files mapping user genes to fitness data (optional)
            fitness_genomes_dir: Directory containing fitness genome JSON files
                                 (e.g., ANA3.json) with condition metadata (optional)
            model_data_dir: Directory containing *_data.json model output files
                            with minimal media fluxes for model predictions (optional)
            reference_phenosim_dir: Directory containing pre-computed phenosim JSON
                                    files for experimental/reference genomes (optional).
                                    These are added to model_performance and
                                    genome_phenotypes tables with source="experiment".
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
        # Normalize genome IDs: replace dots with underscores to match phenosim filenames
        # (e.g., "106654.22" in experimental_data.json -> "106654_22" in phenosim files)
        genome_experiment_data = {}
        if os.path.exists(experiment_data_file):
            with open(experiment_data_file) as f:
                experiment_data = json.load(f)
            for _, item_val in experiment_data.items():
                genome_id = item_val["genome_id"].replace(".", "_")
                if genome_id not in genome_experiment_data:
                    genome_experiment_data[genome_id] = {}
                genome_experiment_data[genome_id][item_val["cpd_id"]] = {
                    "growth": item_val["signal"],
                    "source": item_val["experiment"],
                    "chebi_id": item_val.get("chebi_id", "")
                }

        # === Load fitness condition mappings from fitness genome JSONs ===
        # Build (fitness_genome_id, set_id) -> msid for carbon_source conditions
        setid_to_msid = {}
        if fitness_genomes_dir and os.path.isdir(fitness_genomes_dir):
            for fg_file in os.listdir(fitness_genomes_dir):
                if not fg_file.endswith('.json'):
                    continue
                fg_id = fg_file.replace('.json', '')
                fg_path = os.path.join(fitness_genomes_dir, fg_file)
                try:
                    with open(fg_path) as f:
                        fg_data = json.load(f)
                    for gene_data in fg_data.get("genes", {}).values():
                        for set_id, fitness_entry in gene_data.get("fitness", {}).items():
                            if fitness_entry.get("type") == "carbon_source" and "msid" in fitness_entry:
                                key = (fg_id, set_id)
                                if key not in setid_to_msid:
                                    setid_to_msid[key] = fitness_entry["msid"]
                except (json.JSONDecodeError, IOError) as e:
                    print(f"  Warning: Could not read fitness genome {fg_file}: {e}")
            if setid_to_msid:
                n_orgs = len(set(k[0] for k in setid_to_msid))
                print(f"Loaded {len(setid_to_msid)} condition-to-compound mappings from {n_orgs} fitness genomes")

        # Pre-build the condition mapping DataFrame for efficient joins
        mapping_df = None
        phenotypes_with_fitness_data = set()  # compound IDs that have fitness experiments
        if setid_to_msid:
            mapping_rows = [(fg_id, sid, msid) for (fg_id, sid), msid in setid_to_msid.items()]
            mapping_df = pd.DataFrame(mapping_rows, columns=['fitness_genome_id', 'set_id', 'msid'])
            phenotypes_with_fitness_data = set(mapping_df['msid'])

        # === Load model data for gene-reaction mapping ===
        model_rxn_to_genes = {}   # genome_id -> {rxn_id: set(gene_ids)}
        if model_data_dir and os.path.isdir(model_data_dir):
            for mf in os.listdir(model_data_dir):
                if not mf.endswith('_data.json'):
                    continue
                mf_path = os.path.join(model_data_dir, mf)
                try:
                    with open(mf_path) as f:
                        mdata = json.load(f)
                except (json.JSONDecodeError, IOError):
                    continue
                if not mdata.get('success', False):
                    continue
                mid = mdata.get('model_info', {}).get('model_id', mf.replace('_data.json', ''))
                rxn_genes = {}
                for rxn in mdata.get('reactions', []):
                    gene_rule = rxn.get('gene_reaction_rule', '')
                    if gene_rule:
                        genes = set(re.findall(r'\b(?!and\b|or\b)[A-Za-z]\w*', gene_rule))
                        if genes:
                            rxn_genes[rxn['id']] = genes
                model_rxn_to_genes[mid] = rxn_genes
            if model_rxn_to_genes:
                print(f"Loaded model data for {len(model_rxn_to_genes)} genomes")

        # === Load essential gene IDs for essentiality fraction ===
        essential_gene_ids = None
        if essential_genes_file and os.path.exists(essential_genes_file):
            ess_df = pd.read_csv(essential_genes_file)
            essential_gene_ids = set(ess_df['locusId'].astype(str))
            print(f"Loaded {len(essential_gene_ids)} essential gene IDs from {essential_genes_file}")

        # Collect phenosim files from user directory and reference directory
        # Each entry: (filepath, genome_id, source)
        phenosim_items = []
        phenosim_genome_ids = set()

        if os.path.isdir(phenosim_directory):
            for f in os.listdir(phenosim_directory):
                if f.endswith('_phenosim.json'):
                    gid = f.replace('_phenosim.json', '')
                    source = "user" if gid.startswith("user_") else "pangenome"
                    phenosim_items.append((os.path.join(phenosim_directory, f), gid, source))
                    phenosim_genome_ids.add(gid)

        if reference_phenosim_dir and os.path.isdir(reference_phenosim_dir):
            ref_count = 0
            for f in os.listdir(reference_phenosim_dir):
                if f.endswith('.json'):
                    gid = f.replace('.json', '')
                    if gid not in phenosim_genome_ids:
                        phenosim_items.append((os.path.join(reference_phenosim_dir, f), gid, "experiment"))
                        phenosim_genome_ids.add(gid)
                        ref_count += 1
            if ref_count:
                print(f"Found {ref_count} reference phenosim files from {reference_phenosim_dir}")

        if not phenosim_items:
            print(f"No phenosim JSON files found in {phenosim_directory}")
            return

        print(f"Processing {len(phenosim_items)} total phenosim files ({sum(1 for _,_,s in phenosim_items if s=='user')} user, {sum(1 for _,_,s in phenosim_items if s=='pangenome')} pangenome, {sum(1 for _,_,s in phenosim_items if s=='experiment')} experiment)")
        os.makedirs(output_dir, exist_ok=True)

        model_performance_rows = []
        genome_phenotype_rows = []
        gene_phenotype_rows = []

        for filepath, genome_id, genome_source in phenosim_items:
            try:
                with open(filepath) as f:
                    file_data = json.load(f)
            except (json.JSONDecodeError, IOError) as e:
                print(f"  Warning: Could not read {filepath}: {e}")
                continue

            details = file_data.get("details", {})
            data = file_data.get("data", {})
            exp_data = genome_experiment_data.get(genome_id, {})

            # === Build phenotype_data from details arrays ===
            phenotype_data = {}  # cpd_id -> {class, simulated_objective, ...}
            phenotype_list = details.get("Phenotype", [])
            class_list = details.get("Class", [])
            sim_obj_list = details.get("Simulated objective", [])
            obs_obj_list = details.get("Observed objective", [])
            transport_list = details.get("Transports missing", [])

            for i, cpd_id in enumerate(phenotype_list):
                phenotype_data[cpd_id] = {
                    "class": class_list[i] if i < len(class_list) else "",
                    "simulated_objective": sim_obj_list[i] if i < len(sim_obj_list) else 0,
                    "observed_objective": obs_obj_list[i] if i < len(obs_obj_list) else 0,
                    "transports_added": transport_list[i] if i < len(transport_list) else "",
                }

            # Enrich with per-compound data (gapfilling results, fluxes)
            for cpd_id, cpd_data in data.items():
                if cpd_id == "summary" or not isinstance(cpd_data, dict):
                    continue
                if cpd_id not in phenotype_data:
                    phenotype_data[cpd_id] = {
                        "class": cpd_data.get("class", ""),
                        "simulated_objective": cpd_data.get("objective_value", 0) or 0,
                        "observed_objective": cpd_data.get("experimental_value", 0) or 0,
                        "transports_added": ";".join(cpd_data.get("missing_transports", [])),
                    }
                pd_entry = phenotype_data[cpd_id]
                pd_entry["objective_value"] = cpd_data.get("objective_value", 0) or 0
                pd_entry["gap_count"] = cpd_data.get("gapfill_count", 0) or 0
                pd_entry["gfreactions"] = cpd_data.get("gfreactions") or {}
                pd_entry["fluxes"] = cpd_data.get("fluxes") or {}
                pd_entry["reactions"] = cpd_data.get("reactions") or []
                pd_entry["reaction_count"] = cpd_data.get("reaction_count", len(pd_entry["reactions"]))

            # Enrich experimental genomes with observed data from experiment_data_file
            # Sets observed_objective and recomputes class (P/N -> CP/CN/FP/FN)
            if genome_source == "experiment" and exp_data:
                for cpd_id, pd_entry in phenotype_data.items():
                    if cpd_id in exp_data:
                        growth = exp_data[cpd_id]["growth"]
                        pd_entry["observed_objective"] = growth
                        sim_cls = pd_entry["class"]
                        if sim_cls in ("P", "CP", "FP"):
                            pd_entry["class"] = "CP" if growth else "FP"
                        elif sim_cls in ("N", "CN", "FN"):
                            pd_entry["class"] = "FN" if growth else "CN"

            # ===== TABLE 1: Model Performance =====
            # Use class from details: CP/CN/FP/FN (with exp data) or P/N (without)
            cp = cn = fp = fn = 0
            positive_growth_count = negative_growth_count = 0
            positive_gaps = []
            negative_gaps = []

            for cpd_id, pd_entry in phenotype_data.items():
                cls = pd_entry["class"]
                gap_count = pd_entry.get("gap_count", 0)
                if cls in ("P", "CP", "FP"):
                    positive_growth_count += 1
                    if gap_count:
                        positive_gaps.append(gap_count)
                else:
                    negative_growth_count += 1
                    if gap_count:
                        negative_gaps.append(gap_count)
                if cls == "CP":
                    cp += 1
                elif cls == "CN":
                    cn += 1
                elif cls == "FP":
                    fp += 1
                elif cls == "FN":
                    fn += 1

            total_compared = cp + cn + fp + fn
            accuracy = round((cp + cn) / total_compared, 4) if total_compared > 0 else 0

            model_performance_rows.append({
                "genome_id": genome_id,
                "taxonomy": "",
                "false_positives": fp,
                "false_negatives": fn,
                "true_positives": cp,
                "true_negatives": cn,
                "accuracy": accuracy,
                "positive_growth": positive_growth_count,
                "negative_growth": negative_growth_count,
                "avg_positive_growth_gaps": round(sum(positive_gaps) / len(positive_gaps), 4) if positive_gaps else 0,
                "avg_negative_growth_gaps": round(sum(negative_gaps) / len(negative_gaps), 4) if negative_gaps else 0,
                "closest_user_genomes": "",
                "source": genome_source
            })

            # ===== TABLE 2: Genome Phenotypes (from phenotype_data with transports) =====
            for cpd_id, pd_entry in phenotype_data.items():
                # Only populate closest_experimental_data for user genomes
                closest_exp = ""
                if genome_source == "user" and cpd_id in exp_data:
                    exp_info = exp_data[cpd_id]
                    closest_exp = f"{exp_info['source']}:{'growth' if exp_info['growth'] else 'no_growth'}"

                genome_phenotype_rows.append({
                    "genome_id": genome_id,
                    "phenotype_id": cpd_id,
                    "phenotype_name": cpd_names.get(cpd_id, cpd_id),
                    "class": pd_entry["class"],
                    "simulated_objective": round(pd_entry["simulated_objective"], 6) if pd_entry["simulated_objective"] else 0,
                    "observed_objective": pd_entry["observed_objective"],
                    "gap_count": pd_entry.get("gap_count", 0),
                    "gapfilled_reactions": ";".join(pd_entry.get("gfreactions", {}).keys()),
                    "reaction_count": pd_entry.get("reaction_count", 0),
                    "transports_added": pd_entry.get("transports_added", ""),
                    "closest_experimental_data": closest_exp,
                    "source": genome_source
                })

            # ===== TABLE 3: Gene-Phenotype Associations =====
            # Collect from three sources: gapfill, model predictions, fitness
            # gene_pheno_map: gene_id -> cpd_id -> {association data}
            gene_pheno_map = {}

            def _ensure_gene_pheno(gene_id, cpd_id):
                if gene_id not in gene_pheno_map:
                    gene_pheno_map[gene_id] = {}
                if cpd_id not in gene_pheno_map[gene_id]:
                    gene_pheno_map[gene_id][cpd_id] = {
                        "model_pred_fluxes": [],
                        "sources": set()
                    }
                return gene_pheno_map[gene_id][cpd_id]

            # Source 1: Gapfilled reactions
            for cpd_id, pd_entry in phenotype_data.items():
                gfreactions = pd_entry.get("gfreactions", {})
                if not gfreactions:
                    continue
                for rxn_id, rxn_info in gfreactions.items():
                    gene_id = (rxn_info[1]
                               if isinstance(rxn_info, list) and len(rxn_info) > 1 and rxn_info[1]
                               else None)
                    if gene_id:
                        entry = _ensure_gene_pheno(gene_id, cpd_id)
                        entry["sources"].add("gapfill")

            # Source 2: Model predictions - reactions with flux for the phenotype
            rxn_to_genes = model_rxn_to_genes.get(genome_id, {})
            if rxn_to_genes:
                for cpd_id, pd_entry in phenotype_data.items():
                    pheno_fluxes = pd_entry.get("fluxes", {})
                    if not pheno_fluxes:
                        continue
                    for rxn_id, pheno_flux in pheno_fluxes.items():
                        if abs(pheno_flux) < 1e-6:
                            continue
                        genes = rxn_to_genes.get(rxn_id, set())
                        for gene_id in genes:
                            entry = _ensure_gene_pheno(gene_id, cpd_id)
                            entry["sources"].add("model_prediction")
                            entry["model_pred_fluxes"].append(abs(pheno_flux))

            # Source 3: Fitness data from mapping parquet
            gene_fitness_for_genome = {}  # gene_id -> cpd_id -> {max, min, avg, count}
            gene_essentiality_fraction = {}  # gene_id -> fraction of mapped fitness genes that are essential
            genes_with_fitness_match = set()  # genes matched to fitness genome clusters
            if fitness_mapping_dir and os.path.isdir(fitness_mapping_dir):
                fitness_parquet = os.path.join(
                    fitness_mapping_dir, f"{genome_id}_fitness.parquet")
                if os.path.exists(fitness_parquet):
                    fitness_df = pd.read_parquet(fitness_parquet)
                    # All genes in the parquet are matched to fitness genome clusters
                    genes_with_fitness_match = set(fitness_df['feature_id'].unique())
                    # Extract and remove match-only marker rows
                    match_mask = fitness_df['set_id'] == 'fitness_genome_match'
                    fitness_df = fitness_df[~match_mask]
                    # Compute essentiality fraction per gene from fitness gene mappings
                    if essential_gene_ids is not None and not fitness_df.empty:
                        for fid, grp in fitness_df.groupby('feature_id'):
                            mapped_genes = set(grp['fitness_feature_id'])
                            n_essential = sum(1 for g in mapped_genes if g in essential_gene_ids)
                            gene_essentiality_fraction[fid] = round(n_essential / len(mapped_genes), 6)
                    # Join with condition mapping to get msid (compound ID)
                    if mapping_df is not None and not fitness_df.empty:
                        fitness_with_msid = fitness_df.merge(
                            mapping_df, on=['fitness_genome_id', 'set_id'], how='inner')
                        if not fitness_with_msid.empty:
                            # Aggregate by (feature_id, msid)
                            agg = (fitness_with_msid
                                   .groupby(['feature_id', 'msid'])['value']
                                   .agg(['max', 'min', 'mean', 'count'])
                                   .reset_index())
                            for _, row in agg.iterrows():
                                fid = row['feature_id']
                                if fid not in gene_fitness_for_genome:
                                    gene_fitness_for_genome[fid] = {}
                                gene_fitness_for_genome[fid][row['msid']] = {
                                    'max': row['max'],
                                    'min': row['min'],
                                    'avg': row['mean'],
                                    'count': int(row['count'])
                                }
                            print(f"  Loaded fitness data for {len(gene_fitness_for_genome)} genes in {genome_id}")
            # Merge fitness into gene_pheno_map (with threshold gating)
            for gene_id, cpd_fitness in gene_fitness_for_genome.items():
                for cpd_id, stats in cpd_fitness.items():
                    meets_threshold = stats['min'] <= fitness_threshold
                    if gene_id in gene_pheno_map and cpd_id in gene_pheno_map[gene_id]:
                        # Gene already associated via gapfill/model — always add stats
                        entry = gene_pheno_map[gene_id][cpd_id]
                        entry["fitness_stats"] = stats
                        if meets_threshold:
                            entry["sources"].add("fitness")
                    elif meets_threshold:
                        # New association — only create if threshold met
                        entry = _ensure_gene_pheno(gene_id, cpd_id)
                        entry["sources"].add("fitness")
                        entry["fitness_stats"] = stats

            # Build reverse mapping: gene -> all reactions in the model
            gene_to_all_reactions = {}
            for rxn_id, genes in rxn_to_genes.items():
                for gene_id in genes:
                    if gene_id not in gene_to_all_reactions:
                        gene_to_all_reactions[gene_id] = set()
                    gene_to_all_reactions[gene_id].add(rxn_id)

            # Convert gene_pheno_map to Table 3 rows
            for gene_id, pheno_dict in gene_pheno_map.items():
                for cpd_id, assoc in pheno_dict.items():
                    fitness_stats = assoc.get("fitness_stats", {})
                    model_fluxes = assoc["model_pred_fluxes"]
                    # Fitness values: None when no data, numeric when present
                    if fitness_stats:
                        fitness_max = round(fitness_stats.get('max', 0), 6)
                        fitness_min = round(fitness_stats.get('min', 0), 6)
                        fitness_avg = round(fitness_stats.get('avg', 0), 6)
                        fitness_count = fitness_stats.get('count', 0)
                    else:
                        fitness_max = None
                        fitness_min = None
                        fitness_avg = None
                        fitness_count = None
                    # Fitness match status column
                    has_fitness_score = (gene_id in gene_fitness_for_genome
                                        and cpd_id in gene_fitness_for_genome[gene_id])
                    if has_fitness_score:
                        fitness_match = "has_score"
                    elif cpd_id not in phenotypes_with_fitness_data:
                        fitness_match = "no_fitness_data_for_phenotype"
                    elif gene_id not in genes_with_fitness_match:
                        fitness_match = "no_fitness_ortholog"
                    else:
                        fitness_match = "no_score_for_gene_phenotype"
                    # All reactions associated with this gene in the model
                    mp_rxns = ";".join(sorted(
                        r.replace("_c0", "") for r in gene_to_all_reactions.get(gene_id, set())))
                    gene_phenotype_rows.append({
                        "genome_id": genome_id,
                        "gene_id": gene_id,
                        "phenotype_id": cpd_id,
                        "phenotype_name": cpd_names.get(cpd_id, cpd_id),
                        "association_sources": ";".join(sorted(assoc["sources"])),
                        "model_pred_reactions": mp_rxns,
                        "model_pred_max_flux": round(max(model_fluxes), 6) if model_fluxes else 0,
                        "fitness_match": fitness_match,
                        "fitness_max": fitness_max,
                        "fitness_min": fitness_min,
                        "fitness_avg": fitness_avg,
                        "fitness_count": fitness_count,
                        "essentiality_fraction": gene_essentiality_fraction.get(gene_id)
                    })

        # Add experiment-only genomes to Table 2
        # These are genomes in experimental_data.json that have no phenosim data
        # phenosim_genome_ids already built above when collecting phenosim_items
        exp_only_count = 0
        for exp_genome_id, exp_cpds in genome_experiment_data.items():
            if exp_genome_id in phenosim_genome_ids:
                continue
            for cpd_id, exp_info in exp_cpds.items():
                genome_phenotype_rows.append({
                    "genome_id": exp_genome_id,
                    "phenotype_id": cpd_id,
                    "phenotype_name": cpd_names.get(cpd_id, cpd_id),
                    "class": "",
                    "simulated_objective": 0,
                    "observed_objective": exp_info["growth"],
                    "gap_count": 0,
                    "gapfilled_reactions": "",
                    "reaction_count": 0,
                    "transports_added": "",
                    "closest_experimental_data": f"{exp_info['source']}:{'growth' if exp_info['growth'] else 'no_growth'}",
                    "source": "experiment"
                })
                exp_only_count += 1
        if exp_only_count:
            print(f"Added {exp_only_count} experiment-only phenotype records from {len(genome_experiment_data) - len(phenosim_genome_ids & set(genome_experiment_data.keys()))} genomes")

        # Write TSV files
        if model_performance_rows:
            df = pd.DataFrame(model_performance_rows)
            df.to_csv(os.path.join(output_dir, "model_performance.tsv"), sep="\t", index=False)
            print(f"Saved model_performance.tsv with {len(df)} genomes")

        if genome_phenotype_rows:
            df = pd.DataFrame(genome_phenotype_rows)
            df.to_csv(os.path.join(output_dir, "genome_phenotypes.tsv"), sep="\t", index=False)
            print(f"Saved genome_phenotypes.tsv with {len(df)} records")

        if gene_phenotype_rows:
            df = pd.DataFrame(gene_phenotype_rows)
            df.to_csv(os.path.join(output_dir, "gene_phenotypes.tsv"), sep="\t", index=False)
            print(f"Saved gene_phenotypes.tsv with {len(df)} records")
        else:
            print("No gene-phenotype associations found")

        print(f"Built phenotype tables in {output_dir}")

    def add_model_data_to_genome_table(self, database_path=None, model_data_dir=None):
        """
        Add model-related data to the genome table in the SQLite database,
        or write to TSV if database_path is None.

        Reads *_data.json files produced by run_model_reconstruction and adds
        columns: model_reactions, model_metabolites, model_genes, genome_class,
        energy_gapfill, gs_gapfill, pyruvate_yield.

        Args:
            database_path: Path to the SQLite database. If None, writes a
                           genome_model_data.tsv file to model_data_dir instead.
            model_data_dir: Directory containing *_data.json files.
                            If None, uses self.directory + "/models/"
        """
        if model_data_dir is None:
            model_data_dir = os.path.join(self.directory, "models")

        data_files = []
        if os.path.isdir(model_data_dir):
            data_files = [f for f in os.listdir(model_data_dir) if f.endswith('_data.json')]

        if not data_files:
            print(f"No model data files found in {model_data_dir}")
            return

        # Collect rows from all model data files
        rows = []
        for data_file in data_files:
            filepath = os.path.join(model_data_dir, data_file)
            try:
                with open(filepath) as f:
                    data = json.load(f)
            except (json.JSONDecodeError, IOError) as e:
                print(f"Warning: Could not read {data_file}: {e}")
                continue

            if not data.get('success', False):
                continue

            model_info = data.get('model_info', {})
            genome_id = model_info.get('model_id', data_file.replace('_data.json', ''))

            # Parse pyruvate yield from growth field (format: "Carbon-Pyruvic-Acid:0.1428...")
            growth_str = str(model_info.get('growth', ''))
            pyruvate_yield = 0.0
            if ':' in growth_str:
                try:
                    pyruvate_yield = float(growth_str.split(':')[1])
                except (ValueError, IndexError):
                    pyruvate_yield = 0.0

            rows.append({
                'genome_id': genome_id,
                'model_reactions': model_info.get('num_reactions', 0),
                'model_metabolites': model_info.get('num_metabolites', 0),
                'model_genes': model_info.get('num_genes', 0),
                'genome_class': model_info.get('genome_class', ''),
                'energy_gapfill': model_info.get('core_gapfill', 0),
                'gs_gapfill': model_info.get('gs_gapfill', 0),
                'pyruvate_yield': pyruvate_yield
            })

        if not rows:
            print("No valid model data found")
            return

        # TSV output mode
        if database_path is None:
            df = pd.DataFrame(rows)
            tsv_path = os.path.join(model_data_dir, "genome_model_data.tsv")
            df.to_csv(tsv_path, sep='\t', index=False)
            print(f"Saved genome model data for {len(df)} genomes to {tsv_path}")
            return

        # SQLite output mode
        conn = sqlite3.connect(database_path)
        cursor = conn.cursor()

        # Add missing columns to genome table
        cursor.execute("PRAGMA table_info(genome)")
        existing_cols = {row[1] for row in cursor.fetchall()}

        new_columns = {
            'model_reactions': 'INTEGER',
            'model_metabolites': 'INTEGER',
            'model_genes': 'INTEGER',
            'genome_class': 'TEXT',
            'energy_gapfill': 'INTEGER',
            'gs_gapfill': 'INTEGER',
            'pyruvate_yield': 'REAL'
        }

        for col_name, col_type in new_columns.items():
            if col_name not in existing_cols:
                cursor.execute(f"ALTER TABLE genome ADD COLUMN {col_name} {col_type}")

        updated_count = 0
        for row in rows:
            cursor.execute("""
                UPDATE genome SET
                    model_reactions = ?,
                    model_metabolites = ?,
                    model_genes = ?,
                    genome_class = ?,
                    energy_gapfill = ?,
                    gs_gapfill = ?,
                    pyruvate_yield = ?
                WHERE id = ?
            """, (
                row['model_reactions'],
                row['model_metabolites'],
                row['model_genes'],
                row['genome_class'],
                row['energy_gapfill'],
                row['gs_gapfill'],
                row['pyruvate_yield'],
                row['genome_id']
            ))

            if cursor.rowcount > 0:
                updated_count += 1
                print(f"  Updated genome '{row['genome_id']}' with model data")
            else:
                print(f"  Warning: Genome '{row['genome_id']}' not found in genome table")

        conn.commit()
        conn.close()
        print(f"Model data added to genome table for {updated_count} genomes in {database_path}")

    def build_model_tables(self, database_path=None, model_path=None, phenoset_file=None, data_files: list = None):
        """
        Create genome_reactions table and update feature tables with reaction data.

        Reads model data JSON file(s) and creates a genome_reactions table in the
        SQLite database. Also updates genome_features and pan_genome_features tables
        with per-gene reaction associations and flux/class data.

        Args:
            database_path: Path to the SQLite database. If None, writes TSV files
                           to the same directory as model_path instead.
            model_path: Path to a *_data.json file or directory containing them.
                        If None, uses self.directory + "/models/"
            data_files: List of data files to process overrides model_path
        """
        if data_files is None:
            if model_path is None:
                model_path = os.path.join(self.directory, "models")

            # Collect model data files otherwise use the data files provided
            if os.path.isdir(model_path):
                data_files = [os.path.join(model_path, f)
                              for f in os.listdir(model_path) if f.endswith('_data.json')]
            elif os.path.isfile(model_path):
                data_files = [model_path]

        if not data_files:
            print(f"No model data files found at {model_path}")
            return

        all_reaction_rows = []
        # (genome_id, gene_id) -> aggregated reaction data
        all_gene_data = {}

        for filepath in data_files:
            try:
                with open(filepath) as f:
                    data = json.load(f)
            except (json.JSONDecodeError, IOError) as e:
                print(f"Warning: Could not read {filepath}: {e}")
                continue

            if not data.get('success', False):
                continue

            model_info = data.get('model_info', {})
            genome_id = model_info.get('model_id', '')

            # Build metabolite ID -> name lookup
            met_names = {}
            for met in data.get('metabolites', []):
                met_names[met['id']] = met.get('name', met['id'])

            # Build gapfilling status lookup per reaction
            gf_status = {}
            gapfilled = data.get('gapfilled_reactions', {})
            for rxn_id in gapfilled.get('core', []):
                gf_status[rxn_id] = 'core'
            for rxn_id in gapfilled.get('rich_media_essential', []):
                if rxn_id not in gf_status:
                    gf_status[rxn_id] = 'rich'
            for rxn_id in gapfilled.get('minimal_media', []):
                if rxn_id not in gf_status:
                    gf_status[rxn_id] = 'minimal'

            # Get flux analysis data
            flux = data.get('flux_analysis', {})
            minimal_fluxes = flux.get('minimal_media', {}).get('pfba_fluxes', {})
            minimal_classes = flux.get('minimal_media', {}).get('fva_classes', {})
            rich_fluxes = flux.get('rich_media', {}).get('pfba_fluxes', {})
            rich_classes = flux.get('rich_media', {}).get('fva_classes', {})

            # Process each reaction
            for rxn in data.get('reactions', []):
                raw_rxn_id = rxn['id']
                rxn_id = re.sub(r'_[a-z]\d+$', '', raw_rxn_id)

                # Build equation with metabolite names by replacing IDs
                equation_ids = rxn.get('reaction', '')
                equation_names = equation_ids
                for met_id, met_name in met_names.items():
                    equation_names = equation_names.replace(met_id, met_name)

                # Determine directionality from bounds
                lb = rxn.get('lower_bound', 0)
                ub = rxn.get('upper_bound', 0)
                if lb < 0 and ub > 0:
                    directionality = 'reversible'
                elif lb >= 0 and ub > 0:
                    directionality = 'forward'
                elif lb < 0 and ub <= 0:
                    directionality = 'reverse'
                else:
                    directionality = 'blocked'

                all_reaction_rows.append({
                    'genome_id': genome_id,
                    'reaction_id': rxn_id,
                    'genes': rxn.get('gene_reaction_rule', ''),
                    'equation_names': equation_names,
                    'equation_ids': equation_ids,
                    'directionality': directionality,
                    'upper_bound': ub,
                    'lower_bound': lb,
                    'gapfilling_status': gf_status.get(raw_rxn_id, gf_status.get(rxn_id, 'none')),
                    'rich_media_flux': rich_fluxes.get(raw_rxn_id, 0.0),
                    'rich_media_class': rich_classes.get(raw_rxn_id, ''),
                    'minimal_media_flux': minimal_fluxes.get(raw_rxn_id, 0.0),
                    'minimal_media_class': minimal_classes.get(raw_rxn_id, '')
                })

                # Parse gene IDs from gene_reaction_rule
                gene_rule = rxn.get('gene_reaction_rule', '')
                if gene_rule:
                    tokens = gene_rule.replace('(', ' ').replace(')', ' ').split()
                    gene_ids = set(t for t in tokens if t not in ('and', 'or'))
                    for gene_id in gene_ids:
                        key = (genome_id, gene_id)
                        if key not in all_gene_data:
                            all_gene_data[key] = {
                                'reactions': [],
                                'rich_fluxes': [],
                                'rich_classes': [],
                                'minimal_fluxes': [],
                                'minimal_classes': []
                            }
                        all_gene_data[key]['reactions'].append(rxn_id)
                        all_gene_data[key]['rich_fluxes'].append(abs(rich_fluxes.get(raw_rxn_id, 0.0)))
                        all_gene_data[key]['rich_classes'].append(rich_classes.get(raw_rxn_id, ''))
                        all_gene_data[key]['minimal_fluxes'].append(abs(minimal_fluxes.get(raw_rxn_id, 0.0)))
                        all_gene_data[key]['minimal_classes'].append(minimal_classes.get(raw_rxn_id, ''))

            print(f"Processed {len(data.get('reactions', []))} reactions from {os.path.basename(filepath)}")

        # FVA class priority: essential > variable > blocked
        class_priority = {
            'essential_forward': 3, 'essential_reverse': 3,
            'forward_only': 2, 'reverse_only': 2, 'reversible': 2,
            'blocked': 1, '': 0
        }

        def most_constrained_class(classes):
            best = ''
            best_priority = 0
            for c in classes:
                p = class_priority.get(c, 0)
                if p > best_priority:
                    best_priority = p
                    best = c
            if best in ('essential_forward', 'essential_reverse'):
                return 'essential'
            elif best in ('forward_only', 'reverse_only', 'reversible'):
                return 'variable'
            elif best == 'blocked':
                return 'blocked'
            return best

        # Build per-gene aggregated rows
        gene_rows = []
        for (genome_id, gene_id), gdata in all_gene_data.items():
            gene_rows.append({
                'genome_id': genome_id,
                'gene_id': gene_id,
                'reaction': ';'.join(sorted(set(gdata['reactions']))),
                'rich_media_flux': max(gdata['rich_fluxes']) if gdata['rich_fluxes'] else 0.0,
                'rich_media_class': most_constrained_class(gdata['rich_classes']),
                'minimal_media_flux': max(gdata['minimal_fluxes']) if gdata['minimal_fluxes'] else 0.0,
                'minimal_media_class': most_constrained_class(gdata['minimal_classes'])
            })

        # === Build media compositions table ===
        media_rows = []
        def _get_cpd_name(cpd_id):
            try:
                return self.biochem_db.compounds.get_by_id(cpd_id).name
            except Exception:
                return cpd_id

        # KBase workspace media: minimal and rich
        for media_name, media_label in [("Carbon-Pyruvic-Acid", "minimal"), ("AuxoMedia", "rich")]:
            try:
                media_obj = self.get_object(media_name, ws="KBaseMedia")
                if media_obj:
                    for mc in media_obj["data"]["mediacompounds"]:
                        cpd_id = mc["compound_ref"].split("/")[-1]
                        media_rows.append({
                            "media_id": media_label,
                            "compound_id": cpd_id,
                            "max_uptake": mc["maxFlux"],
                            "compound_name": _get_cpd_name(cpd_id)
                        })
            except Exception as e:
                print(f"Warning: Could not retrieve KBase media {media_name}: {e}")

        # Phenotype set media formulations
        if phenoset_file and os.path.exists(phenoset_file):
            with open(phenoset_file) as f:
                phenoset_data = json.load(f)
            for pheno in phenoset_data.get("phenotypes", []):
                media_id = pheno.get("base_media_name", pheno["id"])
                for cpd_id, cpd_data in pheno.get("base_media", {}).items():
                    media_rows.append({
                        "media_id": media_id,
                        "compound_id": cpd_id,
                        "max_uptake": abs(cpd_data.get("lower_bound", 0)),
                        "compound_name": _get_cpd_name(cpd_id)
                    })

        # TSV output mode
        if database_path is None:
            output_dir = model_path if os.path.isdir(model_path) else os.path.dirname(model_path)

            if all_reaction_rows:
                df = pd.DataFrame(all_reaction_rows)
                rxn_path = os.path.join(output_dir, "genome_reactions.tsv")
                df.to_csv(rxn_path, sep='\t', index=False)
                print(f"Saved genome_reactions.tsv with {len(df)} rows to {rxn_path}")

            if gene_rows:
                df = pd.DataFrame(gene_rows)
                gene_path = os.path.join(output_dir, "gene_reaction_data.tsv")
                df.to_csv(gene_path, sep='\t', index=False)
                print(f"Saved gene_reaction_data.tsv with {len(df)} rows to {gene_path}")

            if media_rows:
                df = pd.DataFrame(media_rows)
                media_path = os.path.join(output_dir, "media_compositions.tsv")
                df.to_csv(media_path, sep='\t', index=False)
                print(f"Saved media_compositions.tsv with {len(df)} rows to {media_path}")

            return

        # SQLite output mode
        conn = sqlite3.connect(database_path)

        # Write genome_reactions table
        if all_reaction_rows:
            df = pd.DataFrame(all_reaction_rows)
            df.to_sql('genome_reactions', conn, if_exists='replace', index=False)
            print(f"Saved genome_reactions table with {len(df)} rows")

        # Update genome_features and pan_genome_features tables
        cursor = conn.cursor()
        for table_name in ['genome_features', 'pan_genome_features']:
            cursor.execute("SELECT name FROM sqlite_master WHERE type='table' AND name=?", (table_name,))
            if not cursor.fetchone():
                continue

            # Add new columns if they don't exist
            cursor.execute(f"PRAGMA table_info([{table_name}])")
            existing_cols = {row[1] for row in cursor.fetchall()}

            new_cols = {
                'reaction': 'TEXT',
                'rich_media_flux': 'REAL',
                'rich_media_class': 'TEXT',
                'minimal_media_flux': 'REAL',
                'minimal_media_class': 'TEXT'
            }

            for col_name, col_type in new_cols.items():
                if col_name not in existing_cols:
                    cursor.execute(f"ALTER TABLE [{table_name}] ADD COLUMN [{col_name}] {col_type}")

            # Determine the feature ID column name
            cursor.execute(f"PRAGMA table_info([{table_name}])")
            cols = [row[1] for row in cursor.fetchall()]
            feature_col = 'feature_id' if 'feature_id' in cols else 'id'
            has_genome_id = 'genome_id' in cols

            # Update each gene's reaction data
            update_count = 0
            for grow in gene_rows:
                if has_genome_id:
                    cursor.execute(f"""
                        UPDATE [{table_name}] SET
                            reaction = ?,
                            rich_media_flux = ?,
                            rich_media_class = ?,
                            minimal_media_flux = ?,
                            minimal_media_class = ?
                        WHERE [{feature_col}] = ? AND genome_id = ?
                    """, (grow['reaction'], grow['rich_media_flux'], grow['rich_media_class'],
                          grow['minimal_media_flux'], grow['minimal_media_class'],
                          grow['gene_id'], grow['genome_id']))
                else:
                    cursor.execute(f"""
                        UPDATE [{table_name}] SET
                            reaction = ?,
                            rich_media_flux = ?,
                            rich_media_class = ?,
                            minimal_media_flux = ?,
                            minimal_media_class = ?
                        WHERE [{feature_col}] = ?
                    """, (grow['reaction'], grow['rich_media_flux'], grow['rich_media_class'],
                          grow['minimal_media_flux'], grow['minimal_media_class'],
                          grow['gene_id']))

                update_count += cursor.rowcount

            conn.commit()
            print(f"Updated {update_count} rows in {table_name} with reaction data")

        if media_rows:
            df = pd.DataFrame(media_rows)
            df.to_sql('media_compositions', conn, if_exists='replace', index=False)
            print(f"Saved media_compositions table with {len(df)} rows")

        conn.close()
        print(f"Model tables built in {database_path}")

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

def run_phenotype_simulation(model_filename,output_filename,data_path,max_phenotypes,kbversion,token):
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

    # Suppress verbose stdout from modelseedpy library calls
    with redirect_stdout(io.StringIO()):
        # Load model
        model = cobra.io.load_json_model(model_filename)
        mdlutl = MSModelUtil(model)

        #Loading the phenotype set from the reference path
        filename = data_path + "/full_phenotype_set.json"
        with open(filename) as f:
            phenoset_data = json.load(f)
        #Setting max phenotypes if specified in the work item
        if max_phenotypes is not None:
            phenoset_data["phenotypes"] = phenoset_data["phenotypes"][:max_phenotypes]
        #Instantiating the phenotype set
        phenoset = MSGrowthPhenotypes.from_dict(phenoset_data)

        # Create a minimal util instance for this worker
        class PhenotypeWorkerUtil(MSReconstructionUtils,MSFBAUtils,MSBiochemUtils):
            def __init__(self, kbversion, token):
                super().__init__(name="PhenotypeWorkerUtil", kb_version=kbversion, token=token)
        pheno_util = PhenotypeWorkerUtil(kbversion=kbversion, token=token)

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


def run_model_reconstruction2(genome_id: str, genome: MSGenome, output_filename, classifier_dir, kb_version):
    print('run_model_reconstruction2', genome_id, len(genome.features), output_filename)

    # Diagnostic: count features with RAST terms
    rast_count = sum(1 for f in genome.features if 'RAST' in f.ontology_terms and f.ontology_terms['RAST'])
    print(f'  Genome {genome_id}: {len(genome.features)} features, {rast_count} with RAST annotations')
    if rast_count == 0:
        return {
            'success': False,
            'error': f"No RAST annotations found. Check that the TSV 'functions' column is being read correctly."
        }

    worker_util = MSReconstructionUtils(kbversion=kb_version, token=token)

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
            'error': f"Model build returned None: {current_output.get('Comments', ['Unknown']) if current_output else ['Unknown']}"
        }

    model = mdlutl.model

    # Save model before gapfilling so it's available even if gapfilling fails
    cobra.io.save_json_model(model, output_filename+"_cobra.json")

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
    if gf_output is not None:
        gf_rxns = gf_output.get('GS GF', 0)
        growth = gf_output.get('Growth', 'Unknown')
    else:
        gf_output, _, _, _ = worker_util.gapfill_metabolic_model(
            mdlutl=mdlutl,
            genome=genome,
            media_objs=[gapfill_media],
            templates=[model.template],
            atp_safe=False,
            objective='bio1',
            minimum_objective=0.01,
            gapfilling_mode="Sequential",
        )
        if gf_output is not None:
            gf_rxns = gf_output.get('GS GF', 0)
            growth = gf_output.get('Growth', 'Unknown')
        else:
            print(f"  Warning: Gapfilling returned None for {genome_id}")

    # Save model (re-save after gapfilling if it succeeded)
    cobra.io.save_json_model(model, output_filename+"_cobra.json")

    genome_class = current_output.get('Class', 'Unknown')
    core_gf = current_output.get('Core GF', 0)

    # Get minimal and rich media for analysis
    minimal_media = gapfill_media  # Use the gapfill media as minimal
    rich_media = worker_util.get_media("KBaseMedia/AuxoMedia")

    # Collect gapfilled reactions by category
    core_gf_rxns = []
    minimal_gf_rxns = []
    rich_gf_rxns = []

    # Identify gapfilled reactions: no gene association and not
    # biomass, exchange, demand, or sink reactions
    all_gf_rxn_ids = {
        rxn.id for rxn in model.reactions
        if not rxn.gene_reaction_rule
        and not any(rxn.id.startswith(p) for p in ('bio', 'EX', 'DM', 'SK'))
    }

    # Categorize using integrated gapfilling solutions from MSModelUtil
    categorized = set()
    if hasattr(mdlutl, 'integrated_gapfillings'):
        for gf_entry in mdlutl.integrated_gapfillings:
            media_obj = gf_entry.get('media', None)
            media_id = media_obj.id if media_obj and hasattr(media_obj, 'id') else ''
            rxn_ids = list(gf_entry.get('new', {}).keys()) + list(gf_entry.get('reversed', {}).keys())
            is_core = 'atp' in media_id.lower()
            for rxn_id in rxn_ids:
                if rxn_id in all_gf_rxn_ids and rxn_id not in categorized:
                    categorized.add(rxn_id)
                    if is_core:
                        core_gf_rxns.append(rxn_id)
                    else:
                        minimal_gf_rxns.append(rxn_id)

    # Any gapfilled reactions not in integrated_gapfillings default to minimal
    for rxn_id in sorted(all_gf_rxn_ids - categorized):
        minimal_gf_rxns.append(rxn_id)

    print(f"  Gapfilled reactions: {len(all_gf_rxn_ids)} total, {len(core_gf_rxns)} core, {len(minimal_gf_rxns)} minimal media")

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

    return output_data


def run_model_reconstruction(input_filename, output_filename, classifier_dir, kbversion, token):
    worker_util = MSReconstructionUtils(kb_version=kbversion, token=token)

    # Clear MSModelUtil cache for this process
    MSModelUtil.mdlutls.clear()

    # Create safe model ID
    genome_id = os.path.splitext(os.path.basename(input_filename))[0]
    print(f'genome_id: {genome_id}')
    # Load features from genome TSV
    gene_df = pd.read_csv(input_filename, sep='\t', index_col=0).to_dict()['functions']

    # Create MSGenome from features
    genome = MSGenome()
    genome.id = genome_id
    genome.scientific_name = genome_id

    ms_features = []
    for feature_id, rast_str in gene_df.items():
        feature = MSFeature(feature_id, "")
        if not pd.isna(rast_str):
            for v in rast_str.split('; '):
                feature.add_ontology_term('RAST', v)
        ms_features.append(feature)
    genome.add_features(ms_features)

    # Diagnostic: count features with RAST terms
    rast_count = sum(1 for f in genome.features if 'RAST' in f.ontology_terms and f.ontology_terms['RAST'])
    print(f'  Genome {genome_id}: {len(genome.features)} features, {rast_count} with RAST annotations')
    if rast_count == 0:
        return {
            'success': False,
            'error': f"No RAST annotations found in {input_filename}. Check that the TSV 'functions' column is being read correctly."
        }

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
            'error': f"Model build returned None: {current_output.get('Comments', ['Unknown']) if current_output else ['Unknown']}"
        }

    model = mdlutl.model

    # Save model before gapfilling so it's available even if gapfilling fails
    cobra.io.save_json_model(model, output_filename+"_cobra.json")

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
    if gf_output is not None:
        gf_rxns = gf_output.get('GS GF', 0)
        growth = gf_output.get('Growth', 'Unknown')
    else:
        gf_output, _, _, _ = worker_util.gapfill_metabolic_model(
            mdlutl=mdlutl,
            genome=genome,
            media_objs=[gapfill_media],
            templates=[model.template],
            atp_safe=False,
            objective='bio1',
            minimum_objective=0.01,
            gapfilling_mode="Sequential",
        )
        if gf_output is not None:
            gf_rxns = gf_output.get('GS GF', 0)
            growth = gf_output.get('Growth', 'Unknown')
        else:
            print(f"  Warning: Gapfilling returned None for {genome_id}")

    # Save model (re-save after gapfilling if it succeeded)
    cobra.io.save_json_model(model, output_filename+"_cobra.json")

    genome_class = current_output.get('Class', 'Unknown')
    core_gf = current_output.get('Core GF', 0)

    # Get minimal and rich media for analysis
    minimal_media = gapfill_media  # Use the gapfill media as minimal
    rich_media = worker_util.get_media("KBaseMedia/AuxoMedia")

    # Collect gapfilled reactions by category
    core_gf_rxns = []
    minimal_gf_rxns = []
    rich_gf_rxns = []

    # Identify gapfilled reactions: no gene association and not
    # biomass, exchange, demand, or sink reactions
    all_gf_rxn_ids = {
        rxn.id for rxn in model.reactions
        if not rxn.gene_reaction_rule
        and not any(rxn.id.startswith(p) for p in ('bio', 'EX', 'DM', 'SK'))
    }

    # Categorize using integrated gapfilling solutions from MSModelUtil
    categorized = set()
    if hasattr(mdlutl, 'integrated_gapfillings'):
        for gf_entry in mdlutl.integrated_gapfillings:
            media_obj = gf_entry.get('media', None)
            media_id = media_obj.id if media_obj and hasattr(media_obj, 'id') else ''
            rxn_ids = list(gf_entry.get('new', {}).keys()) + list(gf_entry.get('reversed', {}).keys())
            is_core = 'atp' in media_id.lower()
            for rxn_id in rxn_ids:
                if rxn_id in all_gf_rxn_ids and rxn_id not in categorized:
                    categorized.add(rxn_id)
                    if is_core:
                        core_gf_rxns.append(rxn_id)
                    else:
                        minimal_gf_rxns.append(rxn_id)

    # Any gapfilled reactions not in integrated_gapfillings default to minimal
    for rxn_id in sorted(all_gf_rxn_ids - categorized):
        minimal_gf_rxns.append(rxn_id)

    print(f"  Gapfilled reactions: {len(all_gf_rxn_ids)} total, {len(core_gf_rxns)} core, {len(minimal_gf_rxns)} minimal media")

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

    return output_data

def generate_ontology_tables(
    input_database: str,
    reference_data_path: str = "/data/",
    source_tables: list = ["genome_features"],
) -> bool:
    """
    Generate ontology tables from one or more database tables.

    Reads features from one or more tables in a SQLite database, maps RAST
    annotations to SEED roles, extracts EC numbers, enriches all ontology
    terms, and writes results back to the same database.

    Args:
        input_database: Path to the SQLite database file.
        reference_data_path: Path to directory containing reference files:
                            - seed.json (RAST -> seed.role mapping)
                            - statements.parquet (labels, definitions, relationships)
                            - kegg_ko_definitions.parquet
                            - cog_definitions.parquet
                            Default: /data/
        source_tables: List of table names in the database to read features from
                       and harvest ontology terms. Each table is processed
                       independently. Default: ["genome_features"]

    Returns:
        True on success, False on failure.

    Output tables (written to input_database):
        - ontology_terms: All ontology terms with labels and definitions
        - ontology_definitions: Ontology prefix definitions
        - ontology_relationships: Term relationships (is_a, enables_reaction)
    """
    import sqlite3
    import re
    import time
    from pathlib import Path
    import pyarrow.parquet as pq

    db_path = Path(input_database)

    # Check if db.sqlite exists
    if not db_path.exists():
        print(f"Warning: db.sqlite not found in {input_database}, skipping ontology generation")
        return False

    print(f"\n{'='*70}")
    print(f"Generating ontology tables for: {input_database}")
    print(f"{'='*70}")

    try:
        # =====================================================================
        # STEP 1: Initialize term collection and load RAST mapper
        # =====================================================================
        ref_path = Path(reference_data_path) / "reference_data/berdl_db/run_20250819_020438/parquet_files"
        seed_json_path = ref_path / "seed.json"

        print("\n1. Loading RAST → seed.role mapper...")
        mapper = None
        if seed_json_path.exists():
            mapper = RASTSeedMapper(str(seed_json_path))
        else:
            print(f"   Warning: seed.json not found at {seed_json_path}")
            print(f"   RAST → seed.role mapping will be skipped")

        terms_by_type = {
            'GO': set(), 'EC': set(), 'KEGG': set(),
            'COG': set(), 'PFAM': set(), 'SO': set(), 'seed.role': set()
        }
        rast_functions = set()
        seed_role_to_label = {}

        # Patterns for extracting ontology term IDs from cell values
        patterns = {
            'GO': re.compile(r'GO:\d+'),
            'EC': re.compile(r'EC:[\d\.-]+'),
            'KEGG': re.compile(r'(?:KEGG:)?K\d{5}'),
            'COG': re.compile(r'COG:(?:COG\d+|[A-Z])'),
            'PFAM': re.compile(r'(?:PFAM:)?PF\d+(?:\.\d+)?'),
            'SO': re.compile(r'SO:\d+'),
            'seed.role': re.compile(r'seed\.role:\d+'),
        }
        ec_in_rast_pattern = re.compile(r'\(EC[:\s]*([\d\.-]+)\)')
        rast_col_candidates = ['rast_function', 'rast_functions', 'functions', 'Annotation:SSO']

        # =====================================================================
        # STEP 2: Load each table and harvest ontology terms
        # =====================================================================
        print(f"\n2. Extracting ontology terms from tables: {source_tables}")
        conn = sqlite3.connect(str(db_path))
        cursor = conn.cursor()
        tables_read = 0

        for table_name in source_tables:
            cursor.execute("SELECT name FROM sqlite_master WHERE type='table' AND name=?", (table_name,))
            if not cursor.fetchone():
                print(f"   Warning: Table '{table_name}' not found in database, skipping")
                continue

            table_df = pd.read_sql_query(f"SELECT * FROM [{table_name}]", conn)
            tables_read += 1
            print(f"   Processing '{table_name}': {len(table_df)} rows, columns: {list(table_df.columns)[:10]}...")

            # Detect RAST function column for this table
            rast_col = None
            for candidate in rast_col_candidates:
                if candidate in table_df.columns:
                    rast_col = candidate
                    break
            if rast_col:
                print(f"     RAST function column: {rast_col}")

            # Extract terms from every column of every row
            for col in table_df.columns:
                for _, row in table_df.iterrows():
                    value = str(row.get(col, ''))
                    if not value or value == 'nan':
                        continue

                    # Match ontology term IDs
                    for ont_type, pattern in patterns.items():
                        matches = pattern.findall(value)
                        for match in matches:
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
                            for separator in [' / ', ' @ ', '; ']:
                                if separator in value:
                                    parts = value.split(separator)
                                    for part in parts:
                                        part = part.strip()
                                        if part:
                                            rast_functions.add(part)
                            if not any(sep in value for sep in [' / ', ' @ ', '; ']):
                                rast_functions.add(value)

        conn.close()

        if tables_read == 0:
            print(f"   Error: No valid tables found in database")
            return False

        # =====================================================================
        # STEP 3: Map RAST functions to seed.role IDs
        # =====================================================================
        if mapper and rast_functions:
            print(f"\n3. Mapping {len(rast_functions)} RAST functions to seed.role IDs...")

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
            print("\n3. Skipping RAST → seed.role mapping (no mapper or no RAST functions)")

        # Summary of extracted terms
        total_terms = sum(len(terms) for terms in terms_by_type.values())
        print(f"\n   Total unique terms: {total_terms}")
        for ont_type, terms in sorted(terms_by_type.items()):
            if terms:
                print(f"     {ont_type}: {len(terms)}")

        if total_terms == 0:
            print("   Warning: No ontology terms found, writing empty tables")
            conn = sqlite3.connect(str(db_path))
            for tbl in ['ontology_terms', 'ontology_definitions', 'ontology_relationships']:
                conn.execute(f"DROP TABLE IF EXISTS [{tbl}]")
            pd.DataFrame(columns=['ontology_prefix', 'identifier', 'label', 'definition']).to_sql(
                'ontology_terms', conn, if_exists='replace', index=False)
            pd.DataFrame(columns=['ontology_prefix', 'definition']).to_sql(
                'ontology_definitions', conn, if_exists='replace', index=False)
            pd.DataFrame(columns=['subject', 'predicate', 'object']).to_sql(
                'ontology_relationships', conn, if_exists='replace', index=False)
            conn.close()
            return True

        # =====================================================================
        # STEP 4: Enrich terms from local parquet files
        # =====================================================================
        print("\n4. Enriching terms from local parquet files...")

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
        # STEP 5: Extract relationships from statements.parquet
        # =====================================================================
        print("\n5. Extracting ontology relationships...")

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
        # STEP 6: Add EC column to ontology terms
        # =====================================================================
        print("\n6. Adding EC column to ontology terms...")

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

        # Extract GO -> EC mapping from statements.parquet (oio:hasDbXref to EC:)
        go_to_ec = {}
        if statements_df is not None:
            print(f"   Extracting GO -> EC mappings from statements.parquet...")
            # Filter for GO terms with hasDbXref predicate
            go_ec_mask = (
                statements_df['subject'].str.startswith('GO:', na=False) &
                (statements_df['predicate'] == 'oio:hasDbXref')
            )
            go_dbxref_df = statements_df[go_ec_mask]

            # Extract EC numbers from object or value column
            ec_pattern = re.compile(r'EC:[\d\.\-]+')
            for _, row in go_dbxref_df.iterrows():
                go_id = row['subject']
                # Check both object and value columns for EC reference
                obj_val = str(row.get('object', '')) + ' ' + str(row.get('value', ''))
                ec_matches = ec_pattern.findall(obj_val)
                if ec_matches:
                    if go_id not in go_to_ec:
                        go_to_ec[go_id] = []
                    go_to_ec[go_id].extend(ec_matches)

            # Deduplicate EC values per GO term
            for go_id in go_to_ec:
                go_to_ec[go_id] = list(set(go_to_ec[go_id]))

            print(f"   Found {len(go_to_ec)} GO terms with EC cross-references")

        # Patterns for extracting EC and TC from labels
        ec_label_pattern = re.compile(r'\(EC\s*([\d\.-]+)\)')
        tc_label_pattern = re.compile(r'\(TC\s*([\d\.\w]+)\)')

        kegg_ec_count = 0
        go_ec_count = 0
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

            elif prefix == 'GO':
                # GO terms: lookup EC from hasDbXref extracted above
                if identifier in go_to_ec:
                    ec_values.extend(go_to_ec[identifier])
                    go_ec_count += 1

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
        print(f"   GO terms with EC: {go_ec_count}")
        print(f"   seed.role with EC: {seed_ec_count}")
        print(f"   seed.role with TC: {seed_tc_count}")
        print(f"   EC terms copied: {ec_copy_count}")

        total_with_ec = sum(1 for t in enriched_terms if t.get('ec'))
        print(f"   Total terms with ec column: {total_with_ec}")

        # =====================================================================
        # STEP 7: Create ontology definitions
        # =====================================================================
        print("\n7. Creating ontology definitions...")

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
        # STEP 8: Write output tables to SQLite
        # =====================================================================
        print(f"\n8. Writing ontology tables to {db_path}...")

        conn = sqlite3.connect(str(db_path))
        for tbl in ['ontology_terms', 'ontology_definitions', 'ontology_relationships']:
            conn.execute(f"DROP TABLE IF EXISTS [{tbl}]")

        # Write ontology_terms
        terms_df = pd.DataFrame(enriched_terms)
        terms_df = terms_df.drop_duplicates(subset=['identifier'])
        terms_df = terms_df.sort_values(['ontology_prefix', 'identifier']).reset_index(drop=True)
        terms_df.to_sql('ontology_terms', conn, if_exists='replace', index=False)
        print(f"   Saved {len(terms_df)} rows to ontology_terms")

        for prefix in terms_df['ontology_prefix'].unique():
            count = len(terms_df[terms_df['ontology_prefix'] == prefix])
            print(f"     {prefix}: {count}")

        # Write ontology_definitions
        defs_df = pd.DataFrame(definition_rows)
        defs_df.to_sql('ontology_definitions', conn, if_exists='replace', index=False)
        print(f"   Saved {len(defs_df)} rows to ontology_definitions")

        # Write ontology_relationships
        rels_df = pd.DataFrame(relationships)
        if not rels_df.empty:
            rels_df = rels_df.drop_duplicates()
        rels_df.to_sql('ontology_relationships', conn, if_exists='replace', index=False)
        print(f"   Saved {len(rels_df)} rows to ontology_relationships")

        if not rels_df.empty:
            print(f"   By predicate:")
            for pred in rels_df['predicate'].unique():
                count = len(rels_df[rels_df['predicate'] == pred])
                print(f"     {pred}: {count}")

        conn.close()

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
