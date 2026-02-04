import os
from berdl.query.query_pangenome_local import QueryPangenomeLocal
from berdl.query.query_genome_local import QueryGenomeLocal
from berdl.pangenome import BERDLPangenome

def main():
    query_pg = QueryPangenomeLocal('/data/reference_data/berdl_db/ke-pangenomes')
    query_g = QueryGenomeLocal()
    berld_pangenome = BERDLPangenome(query_pg, query_g)


if __name__ == "__main__":
    pass