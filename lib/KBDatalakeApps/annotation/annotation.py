import os


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


def parse_bakta(data: str):
    pass


def parse_psortb(data: str):
    annotation = {}
    lines = data.split('\n')
    h = lines[0].strip()
    # skip header
    for line in lines[1:]:
        r = line.strip().split('\t')
        d = {h[i]: r[i].strip() for i in range(len(h))}
        annotation[d['SeqID']] = d

    return annotation


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
        print('parse', parse_psortb(result))
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
            annotation[l_feature_id[i]] = result[i]
    except Exception as ex:
        print(f'nope {ex}')
