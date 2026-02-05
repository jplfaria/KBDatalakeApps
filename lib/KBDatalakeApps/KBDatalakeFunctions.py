import sys
import os
import re
import sqlite3
import json
import uuid
import shutil
from os import path

sys.path = ["/deps/KBUtilLib/src","/deps/cobrakbase","/deps/ModelSEEDpy","/deps/cb_annotation_ontology_api"] + sys.path

# Import utilities with error handling
from kbutillib import KBAnnotationUtils, MSReconstructionUtils,MSFBAUtils,MSBiochemUtils

import pandas as pd
import cobra
from modelseedpy import AnnotationOntology, MSPackageManager, MSTemplateBuilder, MSMedia, MSModelUtil, MSBuilder, MSATPCorrection, MSGapfill, MSGrowthPhenotype, MSGrowthPhenotypes, ModelSEEDBiochem, MSExpression

@staticmethod
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

def genome_to_model(input_filename,output_filename):
    worker_util = MSReconstructionUtils()

    # Clear MSModelUtil cache for this process
    MSModelUtil.mdlutls.clear()

    # Create safe model ID
    genome_id = os.path.splitext(os.path.basename(input_filename))[0]

    # Load features from genome TSV
    gene_df = pd.read_csv(input_filename, sep='\t')

    # Create MSGenome from features
    genome = MSGenome()
    genome.id = genome_id
    genome.scientific_name = genome_id

    ms_features = []
    for _, gene in gene_df.iterrows():
        protein = gene.get('protein_translation', '')
        gene_id = gene.get('gene_id', '')
        if pd.notna(protein) and protein:
            feature = MSFeature(gene_id, str(protein))
            # Parse Annotation:SSO column
            # Format: SSO:nnnnn:description|rxn1,rxn2;SSO:mmmmm:desc2|rxn3
            sso_col = gene.get('functions', '')
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

def run_user_genome_to_tsv(genome_ref,output_filename):
    # Load genome object into object_hash
    anno_util = KBAnnotationUtils()
    anno_util.build_genome_tsv(genome_ref,output_filename)
