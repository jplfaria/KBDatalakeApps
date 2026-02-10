from pathlib import Path
from annotation.annotation import run_rast, run_kofam, run_psortb, run_bakta


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
        out = run_psortb(client, filename_faa, org_arg, output_annotation)

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