from pathlib import Path
import polars as pl
import pandas as pd


class PhenotypePredict:

    # FIXME: not used
    ONTOLOGY_TERM_KEGG_KO = 'KEGG'
    ONTOLOGY_TERM_RAST = 'RAST'

    def __init__(self, model):
        self.model = model
        self.genome_features = {}

    def add_genome(self, path_genome, genome_id):
        features_sso_filter, features_ko_filter = self.build_X(path_genome, genome_id, self.model)
        self.genome_features[genome_id] = features_sso_filter | features_ko_filter

    def predict(self):
        if len(self.genome_features) == 0:
            raise ValueError('not genomes to predict')

        _data_X = {feature: [] for feature in self.model.feature_names_}
        labels = []
        for genome_id, features in self.genome_features.items():
            labels.append(genome_id)
            for k in self.model.feature_names_:
                _data_X[k].append(1 if k in features else 0)

        prediction = self.model.predict(pd.DataFrame(_data_X))
        return prediction

    @staticmethod
    def build_X(path_genome: Path, genome_id, model):
        path_genome_rast = path_genome / f'{genome_id}_rast.tsv'
        path_genome_kofam = path_genome / f'{genome_id}_annotation_kofam.tsv'

        if not path_genome_rast.exists():
            raise FileNotFoundError()
        if not path_genome_kofam.exists():
            raise FileNotFoundError()

        d_rast = {row['feature_id']: row['RAST'] for row in pl.read_csv(path_genome_rast, separator='\t').rows(named=True)}
        d_ko = {row['feature_id']: row['KEGG'] for row in pl.read_csv(path_genome_kofam, separator='\t').rows(named=True)}

        features_sso = set()
        for v in d_rast.values():
            _sso = set()
            if v:
                for s in v.split('; '):
                    _res = t.transform(s)
                    if _res:
                        _sso |= _res
            features_sso |= _sso
        features_ko = set()
        for v in d_ko.values():
            if v:
                features_ko |= set(v.split('; '))

        feature_names = set(model.feature_names_)

        features_sso_filter = features_sso & feature_names
        features_ko_filter = features_ko & feature_names

        return features_sso_filter, features_ko_filter
