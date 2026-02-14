import os
from pathlib import Path
from berdl.prediction.phenotype import PhenotypePredict
from catboost import CatBoostClassifier


def main():
    model = CatBoostClassifier()
    model.load_model(Path('./models/Glucose/Glucose_standard.cbm'))
    path_ml_models = Path('/data/reference_data/ml/phenotype')
    #for
    #PhenotypePredict(model=)


if __name__ == "__main__":
    pass
