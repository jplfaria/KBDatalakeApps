import json
from modelseedpy.core.msgenome import normalize_role


class TransformOntologyRASTToSSO:

    def __init__(self, nmz_value_to_id):
        self.nmz_value_to_id = nmz_value_to_id

    @staticmethod
    def build_from_sso_dictionary(filename_sso_dict_json):
        with open(str(filename_sso_dict_json), 'r') as fh:
            dict_sso = json.load(fh)
        nmz_value_to_id = {}
        id_to_o = {}
        for _id, o in dict_sso['term_hash'].items():
            nmz_str = normalize_role(o['name'])
            if nmz_str not in nmz_value_to_id:
                nmz_value_to_id[nmz_str] = set()
            id_to_o[_id] = o
            nmz_value_to_id[nmz_str].add(_id)

        return TransformOntologyRASTToSSO(nmz_value_to_id)

    def transform(self, value):
        nmz_str = normalize_role(value)
        _values = self.nmz_value_to_id.get(nmz_str)

        return _values
