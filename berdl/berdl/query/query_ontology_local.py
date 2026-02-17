from pathlib import Path
import polars as pl
from berdl.query.query_ontology import QueryOntologyABC
from berdl.hash_seq import ProteinSequence
from modelseedpy.core.msgenome import MSFeature, MSGenome


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


class QueryOntologyLocal(QueryOntologyABC):

    def __init__(self, root=Path('/data/reference_data/berdl_db')):
        self.root = root
        self.ldf_annotation_bakta = pl.scan_parquet(root / 'annotation/bakta/*.parquet')
        self.ldf_annotation_kofam = pl.scan_parquet(root / 'annotation/kofam/*.parquet')

    def get_protein_ontology(self, seq: str):
        protein = ProteinSequence(seq)
        h = protein.hash_value

        annotation_bakta = self.ldf_annotation_bakta.filter(pl.col("_id") == h).collect()
        annotation_kofam = self.ldf_annotation_kofam.filter(pl.col("_id") == h).collect()

        return {
            'bakta': annotation_bakta,
            'kofam': annotation_kofam,
        }

    def get_annotation_bakta(self, hash_collection):
        res = self.ldf_annotation_bakta.filter(pl.col("_id").is_in(set(hash_collection))).collect()
        res = {row['_id']: row for row in res.rows(named=True)}

        return res

    def get_annotation_kofam(self, hash_collection):
        res = self.ldf_annotation_kofam.filter(pl.col("_id").is_in(set(hash_collection))).collect()
        res = {row['_id']: row for row in res.rows(named=True)}

        return res

    @staticmethod
    def clean_bakta_value(ont, value):
        if ont == 'KEGG' and value.startswith('KEGG:'):
            return value[5:]
        if ont == 'COG' and value.startswith('COG:'):
            return value[4:]
        if ont in {'uniref_90', 'uniref_100', 'uniref_50'} and value.startswith('UniRef:'):
            return value[7:]
        return value

    def get_annotations(self, genome: MSGenome):
        ret = {}
        h_to_feature_id = {}
        feature_id_to_h = {}
        for f in genome.features:
            if f.seq:
                protein = ProteinSequence(f.seq)
                h = protein.hash_value
                if h not in h_to_feature_id:
                    h_to_feature_id[h] = set()
                h_to_feature_id[h].add(f.id)
                feature_id_to_h[f.id] = h

        h_to_kofam = self.get_annotation_kofam(h_to_feature_id)

        data = {
            'feature_id': [],
            'KEGG': [],
        }
        for f in genome.features:
            if f.seq:
                h = feature_id_to_h[f.id]
                _annotation = h_to_kofam.get(h)
                if _annotation:
                    kos = _annotation.get('kos', [])
                    if kos and len(kos) > 0:
                        data['feature_id'].append(f.id)
                        data['KEGG'].append('; '.join(kos))
        ret['kofam'] = pl.DataFrame(data)

        h_to_batka = self.get_annotation_bakta(h_to_feature_id)

        data = {
            'feature_id': [],
        }

        def list_tuple_to_dict(list_tuple, ontology_filster: set):
            _res = {_v: [] for _v in ontology_filster}
            for ont, value in list_tuple:
                _res[ont].append(self.clean_bakta_value(ont, value))
            return _res
        bakta_ontology_set = set()
        for doc in h_to_batka.values():
            o = collect_ontology(doc)
            bakta_ontology_set |= {t[0] for t in o}
        for v in bakta_ontology_set:
            data[v] = []

        for f in genome.features:
            if f.seq:
                h = feature_id_to_h[f.id]
                doc = h_to_batka.get(h)
                if doc:
                    o = collect_ontology(doc)
                    dict_list = list_tuple_to_dict(o, bakta_ontology_set)
                    data['feature_id'].append(f.id)
                    for k, l in dict_list.items():
                        data[k].append('; '.join(l) if len(l) > 0 else None)
        ret['bakta'] = pl.DataFrame(data)

        return ret

    def get_protein_ontology_bulk(self, features_or_seqs: list):
        _all_h = set()
        index_to_h = {}
        for i, f_or_s in enumerate(features_or_seqs):
            seq = f_or_s.seq if isinstance(f_or_s, MSFeature) else f_or_s
            protein = ProteinSequence(seq)
            h = protein.hash_value
            _all_h.add(h)
            index_to_h[i] = h

        collected_annotation = {i: {} for i in range(len(features_or_seqs))}
        # collect bakta ontology
        res = self.ldf_annotation_bakta.filter(pl.col("_id").is_in(_all_h)).collect()
        res = {row['_id']: row for row in res.rows(named=True)}

        for i, h in index_to_h.items():
            feature_annotation_bakta = res.get(h)
            if feature_annotation_bakta:
                feature_ontology = collect_ontology(feature_annotation_bakta)
                for on, v in feature_ontology:
                    if on not in collected_annotation[i]:
                        collected_annotation[i][on] = set()
                    collected_annotation[i][on].add(v)

        # collect kofam ontology
        res = self.ldf_annotation_kofam.filter(pl.col("_id").is_in(_all_h)).collect()
        res = {row['_id']: row['kos'] for row in res.rows(named=True)}

        for i, h in index_to_h.items():
            feature_annotation_kofam = res.get(h, [])
            if feature_annotation_kofam:
                if 'KEGG' not in collected_annotation[i]:
                    collected_annotation[i]['KEGG'] = set()
                for ko in feature_annotation_kofam:
                    collected_annotation[i]['KEGG'].add(ko)

        return collected_annotation
