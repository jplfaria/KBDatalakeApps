from query_genome import QueryGenomeABC


class QueryGenomeLocal(QueryGenomeABC):

    def __init__(self):
        pass

    def get_genome_features(self, genome_id: str):
        """Return all gene clusters for a given clade."""
        raise NotImplementedError

    def get_genome_contigs(self, genome_id: str):
        """Return all gene clusters for a given clade."""
        raise NotImplementedError
