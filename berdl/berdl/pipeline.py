import os
import argparse
from pathlib import Path
import json
from berdl.genome_paths import GenomePaths
from berdl.prep_genome_set import BERDLPreGenome
from cobrakbase import KBaseAPI


def main(input_params):
    print(os.environ)
    #  setup clients/methods
    kbase = KBaseAPI(input_params['_ctx']['token'], dev=True)
    #kbase = KBaseAPI(input_params['_ctx']['token'], config={
    #    'workspace-url': input_params['_config']['workspace-url']})

    #  extract genomes
    genomes = []
    for ref in input_params['input_refs']:
        kbase_input_object = kbase.get_from_ws(ref)
        kbase_input_object_type = kbase_input_object.info.type

        print('input_object is:', kbase_input_object_type)
        if kbase_input_object_type == 'KBaseGenomes.Genome':
            genomes.append(kbase.get_from_ws(str(ref)))
        elif kbase_input_object_type == 'wololo2':
            pass
        else:
            pass
            # raise ValueError('')

    paths = GenomePaths(root=Path(input_params['_config']['scratch']).resolve())
    berdl_prep_genomes = BERDLPreGenome(kbase, paths)
    user_to_clade, ani_clade, df_ani_fitness, df_ani_phenotype = berdl_prep_genomes.run(genomes)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Run BERDL genome pipeline"
    )
    parser.add_argument(
        "input_params",
        help="Path to input params JSON file"
    )
    #  read input params
    args = parser.parse_args()
    filename_input_params = args.input_params

    if not os.path.exists(filename_input_params):
        raise FileNotFoundError(
            f"Input params file not found: {filename_input_params}"
        )

    with open(filename_input_params, 'r') as fh:
        input_params = json.load(fh)

    main(input_params)
