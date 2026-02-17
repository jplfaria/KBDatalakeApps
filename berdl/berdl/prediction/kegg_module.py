import json


class KeggModules:

    def __init__(self, kegg_pwy_def='/scratch/fliu/old_jupyter/modelseed2/data/modelseed2-data/data/kegg_pwy_def.json',
                 kegg_module_def='/scratch/fliu/old_jupyter/modelseed2/data/modelseed2-data/data/kegg_module_def.json',
                 kegg_to_seed='/scratch/fliu/data/KE/kegg_to_seed.json'):
        self.pathway_names = {}
        self.pathway_modules = {}
        self.pathway_module_data = {}
        with open(kegg_pwy_def, 'r') as fh:
            self.pathway_names = json.load(fh)
        with open(kegg_module_def, 'r') as fh:
            self.pathway_module_data = json.loads(fh.read())
        for module_id, module in self.pathway_module_data.items():
            for pwy in module['pathways']:
                if pwy not in self.patshway_modules:
                    self.pathway_modules[pwy] = set()
                self.pathway_modules[pwy].add(module_id)

        self.kegg_to_seed = None
        with open(kegg_to_seed, 'r') as fh:
            self.kegg_to_seed = {k: set(v) for k, v in json.load(fh).items()}
        self.seed_to_kegg = {}
        for kegg_id in self.kegg_to_seed:
            for seed_id in self.kegg_to_seed[kegg_id]:
                if seed_id not in self.seed_to_kegg:
                    self.seed_to_kegg[seed_id] = set()
                self.seed_to_kegg[seed_id].add(kegg_id)
        self.ignore_pwy = {
            'map01100',  # Metabolic pathways
            'map01110',  # Biosynthesis of secondary metabolites
            'map01120',  # Microbial metabolism in diverse environments
            'map01200',  # Carbon metabolism
            'map01230',  # Biosynthesis of amino acids
        }
        self.map_order = [
            "map00010", "map00030", "map00020", "map00620", "map00190", "map00630", "map00720", "map00710", "map00195",
            "map00230", "map00240", "map00040",
            "map00500", "map00051", "map00052", "map00966", "map00562", "map00053", "map00660", "map00260", "map00250",
            "map00270", "map00220", "map00340",
            "map00460", "map00330", "map00360", "map00400", "map00300", "map00350", "map00290", "map00380", "map00480",
            "map00310", "map00280", "map01220",
            "map00900", "map00100", "map00920", "map00910", "map00680", "map00650", "map00640", "map00072", "map00430",
            "map00627", "map00362", "map00790",
            "map00670", "map00410", "map00770", "map00860", "map00740", "map00730", "map00760", "map00750", "map00130",
            "map00780", "map01212", "map00061",
            "map00062", "map00071", "map01040", "map00592", "map00565", "map00604", "map00603", "map00601", "map00540",
            "map00600", "map00561", "map00564",
            "map00513", "map00510", "map00512", "map00515", "map00532", "map00534", "map00531", "map00563", "map01210",
            "map00520",
            # "map00941", "map00906", "map00940", "map00403", "map00522", "map00404", "map00140", "map01059", "map00333",
            "map00521",  # Streptomycin biosynthesis
            # "map00261", "map00311", "map00525", "map00332",
            # "map00905", "map00405", "map00331", "map00120", "map00253", "map00361", "map00621", "map00626", "map00623", "map00624", "map00622", "map00999",
            # "map00998", "map00997", "map01056", "map01057", "map00523", "map02020", "map01502", "map01503", "map01501", "map05111", "map05120", "map05130",
            # "map05133", "map05131", "map05110"
        ]
        self.module_order = []
        self.module_alias = {}
        for map_id in self.map_order:
            for module_id in self.pathway_modules[map_id]:
                if module_id not in self.module_order:
                    self.module_order.append(module_id)
                    # print(map_id, pathway_names[map_id])
                    self.module_alias[module_id] = self.pathway_names[map_id]

    def score(self, module, kegg_ids):
        c_total = 0
        c_match = 0
        for k in module['reactions']:
            c_total += 1
            rxn_set = set(k[1:-1].split(', '))
            if len(rxn_set & kegg_ids) > 0:
                c_match += 1
        if c_total == 0:
            return None
        score = c_match / c_total
        return score
