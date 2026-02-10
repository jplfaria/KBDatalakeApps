import os
from pathlib import Path
from modelseedpy import MSGenome


def collect_ontology(_doc):
    ontology = []
    ontology.append(['bakta_product', _doc['product']])
    for db_xref in _doc.get('db_xrefs', []):
        if db_xref.startswith('UniRef:UniRef50'):
            ontology.append(['uniref_50', db_xref])  # not defined in berld ontology
        elif db_xref.startswith('UniRef:UniRef90'):
            ontology.append(['uniref_90', db_xref])  # not defined in berld ontology
        elif db_xref.startswith('UniRef:UniRef100'):
            ontology.append(['uniref_100', db_xref])  # not defined in berld ontology
        elif db_xref.startswith('SO:'):
            ontology.append(['SO', db_xref])
        elif db_xref.startswith('EC:'):
            ontology.append(['EC', db_xref])
        elif db_xref.startswith('KEGG:'):
            ontology.append(['KEGG', db_xref])  # not defined in berld ontology
        elif db_xref.startswith('GO:'):
            ontology.append(['GO', db_xref])
        elif db_xref.startswith('COG:'):
            ontology.append(['COG', db_xref])
        elif db_xref.startswith('PFAM:'):
            ontology.append(['PFAM', db_xref])
        elif db_xref.startswith('UniRef:UniRef100:'):
            ontology.append(['uniref_100', db_xref])
        else:
            ontology.append(['others', db_xref])
    return ontology


def parse_kofam(data: str):
    annotation = {}
    for line in data.split('\n'):
        if line:
            _parts = line.strip().split()
            if len(_parts) == 2:
                feature_id, ko = _parts
                if feature_id not in annotation:
                    annotation[feature_id] = set()
                annotation[feature_id].add(ko)

    return annotation


def parse_bakta(data: dict):
    annotation = {}
    features = data.get("features", [])
    for doc in features:
        i = doc['id']
        annotation[i] = {}

        feature_ontology = collect_ontology(doc)
        for on, v in feature_ontology:
            if on not in annotation[i]:
                annotation[i][on] = set()
            annotation[i][on].add(v)
    return annotation


def parse_psortb(data: str):
    annotation = {}
    lines = data.split('\n')
    h = lines[0].strip().split('\t')
    print(h)
    # skip header
    for line in lines[1:]:
        if line:
            r = line.strip().split('\t')
            print(r)
            d = {h[i]: r[i].strip() for i in range(len(h))}
            annotation[d['SeqID']] = d

    return annotation


def run_rast(client_rast, genome_file_input, output_file):
    print(f'run_rast {genome_file_input} -> {output_file}')
    genome = MSGenome.from_fasta(str(genome_file_input))
    proteins = {f.id: f.seq for f in genome.features if f.seq}
    l_sequences = []
    l_feature_id = []
    for i, s in proteins.items():
        l_sequences.append(s)
        l_feature_id.append(i)

    result = client_rast.annotate_proteins({'proteins': l_sequences})
    annotation = {}
    for i in range(len(l_feature_id)):
        annotation[l_feature_id[i]] = result['functions'][i]
    print('write: ', str(output_file))
    with open(str(output_file), 'w') as fh:
        fh.write('feature_id\tRAST\n')
        for feature_id, list_rast in annotation.items():
            _str = '; '.join(list_rast)
            fh.write(f'{feature_id}\t{_str}\n')


def run_kofam(client_kofam, genome_file_input, output_file):
    print(f'run_kofam {genome_file_input} -> {output_file}')
    genome = MSGenome.from_fasta(str(genome_file_input))
    proteins = {f.id: f.seq for f in genome.features if f.seq}

    result = client_kofam.annotate_proteins(proteins)
    annotation = parse_kofam(result)
    print('write: ', str(output_file))
    with open(str(output_file), 'w') as fh:
        fh.write('feature_id\tKO\n')
        for feature_id, ko_set in annotation.items():
            ko_str = '; '.join(ko_set)
            fh.write(f'{feature_id}\t{ko_str}\n')


def run_psortb(client, org_flag, genome_file_input, output_file):
    """
    org_flag: -n -p -a
    """
    print(f'run_psortb {org_flag} {genome_file_input} -> {output_file}')
    genome = MSGenome.from_fasta(str(genome_file_input))
    proteins = {f.id: f.seq for f in genome.features if f.seq}

    result = client.annotate_proteins(proteins, org_flag)
    annotation = parse_psortb(result)
    print('write: ', str(output_file))
    with open(str(output_file), 'w') as fh:
        fh.write('feature_id\tprimary_localization_psortb\tsecondary_localization_psortb\n')
        for feature_id, d in annotation.items():
            pri_loc = d.get('Final_Localization', '')
            sec_loc = d.get('Secondary_Localization', '')
            fh.write(f'{feature_id}\t{pri_loc}\tsec_loc\n')
    pass


def run_bakta(client, genome_file_input, output_file):
    print(f'run_bakta {genome_file_input} -> {output_file}')
    genome = MSGenome.from_fasta(str(genome_file_input))
    proteins = {f.id: f.seq for f in genome.features if f.seq}

    result = client.annotate_proteins(proteins)
    annotation = parse_bakta(result)
    print('write: ', str(output_file))
    keys = set()
    for v in annotation.values():
        keys |= v
    keys = sorted(keys)
    with open(str(output_file), 'w') as fh:
        header = "\t".join(keys)
        fh.write(f'feature_id\t{header}\n')
        for feature_id, ontology_set in annotation.items():
            values = []
            for k in keys:
                values.append('; '.join(ontology_set.get(k, [])))
            _str_values = '\t'.join(values)
            fh.write(f'{feature_id}\t{_str_values}\n')


def test_annotation(client_kofam, client_bakta, client_psortb, client_rast):
    import time
    proteins = {
        # "tRNA:Cm32/Um32 methyltransferase"
        "Test3.CDS.1": "LFILTATGNMSLCGLKKECLIAASELVTCRE",
        # Aspartokinase (EC 2.7.2.4);Homoserine dehydrogenase (EC 1.1.1.3)
        "Test3.CDS.2": "MRVLKFGGTSVANAERFLRVADILESNARQGQVATVLSAPAKITNHLVAMIEKTISGQDALPNISDAERIFAELLTGLAAAQPGFPLAQLKTFVDQEFAQIKHVLHGISLLGQCPDSINAALICRGEKMSIAIMAGVLEARGHNVTVIDPVEKLLAVGHYLESTVDIAESTRRIAASRIPADHMVLMAGFTAGNEKGELVVLGRNGSDYSAAVLAACLRADCCEIWTDVDGVYTCDPRQVPDARLLKSMSYQEAMELSYFGAKVLHPRTITPIAQFQIPCLIKNTGNPQAPGTLIGASRDEDELPVKGISNLNNMAMFSVSGPGMKGMVGMAARVFAAMSRARISVVLITQSSSEYSISFCVPQSDCVRAERAMQEEFYLELKEGLLEPLAVTERLAIISVVGDGMRTLRGISAKFFAALARANINIVAIAQGSSERSISVVVNNDDATTGVRVTHQMLFNTDQVIEVFVIGVGGVGGALLEQLKRQQSWLKNKHIDLRVCGVANSKALLTNVHGLNLENWQEELAQAKEPFNLGRLIRLVKEYHLLNPVIVDCTSSQAVADQYADFLREGFHVVTPNKKANTSSMDYYHQLRYAAEKSRRKFLYDTNVGAGLPVIENLQNLLNAGDELMKFSGILSGSLSYIFGKLDEGMSFSEATTLAREMGYTEPDPRDDLSGMDVARKLLILARETGRELELADIEIEPVLPAEFNAEGDVAAFMANLSQLDDLFAARVAKARDEGKVLRYVGNIDEDGVCRVKIAEVDGNDPLFKVKNGENALAFYSHYYQPLPLVLRGYGAGNDVTAAGVFADLLRTLSWKLGV",
    }
    try:
        print(f"test kb_kofam annotation")
        start_time = time.perf_counter()
        result = client_kofam.annotate_proteins(proteins)
        end_time = time.perf_counter()
        print(f"Execution time: {end_time - start_time} seconds")
        print(f'received results of type {type(result)} and size {len(result)}')
        print(result)
        print('parse', parse_kofam(result))
    except Exception as ex:
        print(f'nope {ex}')

    try:
        print(f"test kb_bakta annotation")
        start_time = time.perf_counter()
        result = client_bakta.annotate_proteins(proteins)
        end_time = time.perf_counter()
        print(f"Execution time: {end_time - start_time} seconds")
        print(f'received results of type {type(result)} and size {len(result)}')
        print(result)
        print('parse', parse_bakta(result))
    except Exception as ex:
        print(f'nope {ex}')

    try:
        print(f"test kb_psortb annotation")
        start_time = time.perf_counter()
        result = client_psortb.annotate_proteins(proteins, "-n")
        end_time = time.perf_counter()
        print(f"Execution time: {end_time - start_time} seconds")
        print(f'received results of type {type(result)} and size {len(result)}')
        print(result)
        print('parse', parse_psortb(result))
    except Exception as ex:
        print(f'nope {ex}')

    try:
        print(f"test RAST_SDK annotation")
        start_time = time.perf_counter()
        l_sequences = []
        l_feature_id = []
        for i, s in proteins.items():
            l_sequences.append(s)
            l_feature_id.append(i)
        result = client_rast.annotate_proteins({'proteins': l_sequences})
        end_time = time.perf_counter()
        print(f"Execution time: {end_time - start_time} seconds")
        print(f'received results of type {type(result)} and size {len(result)}')
        print(result)
        annotation = {}
        for i in range(len(l_feature_id)):
            annotation[l_feature_id[i]] = result['functions'][i]
        print('parse', annotation)
    except Exception as ex:
        print(f'nope {ex}')
