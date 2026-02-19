from pathlib import Path
from annotation.annotation import run_rast, run_kofam, run_psortb, run_bakta
import pandas as pd


def task_kofam(filename_faa: Path, client):
    print(f'start task_kofam: {filename_faa}')
    if filename_faa.exists() and filename_faa.is_file():
        _parent = filename_faa.parent
        output_annotation = _parent / f'{filename_faa.name[:-4]}_annotation_kofam.tsv'
        # run task
        run_kofam(client, filename_faa, output_annotation)

        print(f'done task_kofam: {output_annotation}')
        return output_annotation
    else:
        raise ValueError(f"invalid filename_faa: {filename_faa}")


def task_rast(filename_faa: Path, client):
    print(f'start task_rast: {filename_faa}')
    if filename_faa.exists() and filename_faa.is_file():
        _parent = filename_faa.parent
        output_annotation = _parent / f'{filename_faa.name[:-4]}_rast.tsv'
        # run task
        run_rast(client, filename_faa, output_annotation)

        print(f'done task_rast: {output_annotation}')
        return output_annotation
    else:
        raise ValueError(f"invalid filename_faa: {filename_faa}")


def task_psortb(filename_faa: Path, org_arg, client):
    print(f'start task_psortb: {filename_faa}')
    if filename_faa.exists() and filename_faa.is_file():
        _parent = filename_faa.parent
        output_annotation = _parent / f'{filename_faa.name[:-4]}_annotation_psortb.tsv'
        # run task
        out = run_psortb(client, org_arg, filename_faa, output_annotation)

        print(f'done task_psortb: {output_annotation}')
        print(out)
        return output_annotation
    else:
        raise ValueError(f"invalid filename_faa: {filename_faa}")


def task_bakta(filename_faa: Path, client):
    print(f'start task_bakta: {filename_faa}')
    if filename_faa.exists() and filename_faa.is_file():
        _parent = filename_faa.parent
        output_annotation = _parent / f'{filename_faa.name[:-4]}_annotation_bakta.tsv'
        # run task
        out = run_bakta(client, filename_faa, output_annotation)

        print(f'done task_bakta: {output_annotation}')
        print(out)
        return output_annotation
    else:
        raise ValueError(f"invalid filename_faa: {filename_faa}")
    

def run_RAST_annotation(input_filepath, output_filename, rast_client):
        sequence_hash = {}
        current_id = None
        current_seq = []
        with open(input_filepath) as f:
            for line in f:
                line = line.strip()
                if line.startswith(">"):
                    if current_id:
                        sequence_hash[current_id] = "".join(current_seq)
                    current_id = line[1:].split()[0]
                    current_seq = []
                elif line:
                    current_seq.append(line)
        if current_id:
            sequence_hash[current_id] = "".join(current_seq)
        proteins = []
        ids = []
        for id, sequence in sequence_hash.items():
            proteins.append(sequence)
            ids.append(id)

        annotate_protein_params = {'proteins': proteins}
        print('annotate_protein_params:', annotate_protein_params)

        result = rast_client.annotate_proteins(annotate_protein_params)
        print('rast annotation result', result)
        functions_list = result.get('functions', [])
        records = []
        for id, functions in zip(ids, functions_list):
            records.append({
                'id': id,
                'functions': functions
            })
        df = pd.DataFrame.from_records(records)
        df.to_csv(output_filename, sep='\t', index=False)
