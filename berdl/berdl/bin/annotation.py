import os
import argparse
from pathlib import Path
from berdl.query.query_ontology_local import QueryOntologyLocal
from modelseedpy import MSGenome


def main(filename_faa):
    path_faa = Path(filename_faa).resolve()
    if not path_faa.exists() or not path_faa.is_file():
        raise ValueError(f'invalid faa path: {path_faa}')

    path_output_folder = path_faa.parent
    genome = MSGenome.from_fasta(str(path_faa))
    query_ontology = QueryOntologyLocal()

    genome_id = path_faa.name[:-4]

    annotations = query_ontology.get_annotations(genome)
    for method, df in annotations.items():
        output = path_output_folder / f'{genome_id}_annotation_{method}.tsv'
        df.write_csv(output, separator='\t')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Run BERDL annotation pipeline"
    )
    parser.add_argument(
        "filename_faa",
        help="Path to Genome Protein Fasta"
    )
    #  read input params
    args = parser.parse_args()
    input_filename_faa = args.filename_faa

    if not os.path.exists(input_filename_faa):
        raise FileNotFoundError(
            f"FASTA file not found: {input_filename_faa}"
        )

    main(input_filename_faa)
