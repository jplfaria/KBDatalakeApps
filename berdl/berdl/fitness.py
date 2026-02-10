import os
import json
import polars as pl
from pathlib import Path
from berdl.hash_seq import ProteinSequence


def map_protein_hash_to_fitness_records(fitness_path='/data/reference_data/phenotype_data/fitness_genomes/'):
    m_to_fitness_feature = {}
    for genome_id in os.listdir(fitness_path):
        if genome_id.endswith('.json'):
            with open(fitness_path + genome_id, 'r') as fh:
                fitness_genome_data = json.load(fh)
                for gene in fitness_genome_data['genes']:
                    _gene_data = fitness_genome_data['genes'][gene]
                    seq = _gene_data.get('protein_sequence')
                    _data_fitness = _gene_data.get('fitness')
                    if _data_fitness and len(_data_fitness) > 0 and seq:
                        protein = ProteinSequence(seq)
                        h = protein.hash_value
                        if h not in m_to_fitness_feature:
                            m_to_fitness_feature[h] = {}
                        m_to_fitness_feature[h][(genome_id[:-5], gene)] = {k: (v['fit'], v['t']) for k, v in
                                                                           _data_fitness.items() if 'fit' in v}
    return m_to_fitness_feature


def create_genome_fitness_table(input_genome, input_genome_id, m_to_r, r_to_m, m_to_fitness_feature):
    data = {
        'genome_id': [],
        'feature_id': [],
        'feature_h': [],
        'fitness_genome_id': [],
        'fitness_feature_id': [],
        'fitness_feature_h': [],
        'set_id': [],
        'value': [],
    }
    for feature in input_genome.features:
        if feature.seq:
            protein = ProteinSequence(feature.seq)
            h = protein.hash_value
            r = m_to_r[h]
            other_m = r_to_m[r]
            for other_h in other_m:
                genome_gene_fitness = m_to_fitness_feature.get(other_h, dict())
                for (fitness_genome_id, fitness_feature_id), _sets in genome_gene_fitness.items():
                    for set_id, (fit, t) in _sets.items():
                        data['genome_id'].append(input_genome_id)
                        data['feature_id'].append(feature.id)
                        data['feature_h'].append(h)
                        data['fitness_genome_id'].append(fitness_genome_id)
                        data['fitness_feature_id'].append(fitness_feature_id)
                        data['fitness_feature_h'].append(other_h)
                        data['set_id'].append(set_id)
                        data['value'].append(fit)

    df = pl.DataFrame(data)
    return df
