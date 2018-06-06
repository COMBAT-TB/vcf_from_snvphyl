import json
import os
import os.path
from py2neo import Graph

DBHOST = os.environ.get("DATABASE_URI", "combattb.sanbi.ac.za")

graph = Graph(host=DBHOST, bolt=True, password="")


def get_gene_data(cache_filename=None):
    if not cache_filename:
        cache_filename = os.path.join(os.environ.get('TMPDIR', '/tmp'), 'combattb_gene_list.cache')
    if os.path.exists(cache_filename):
        data = json.load(open(cache_filename))
    else:
        data = graph.run(
            "MATCH (g:Gene) -[:LOCATED_AT]->(l:Location)"
            "RETURN g.uniquename as uniquename, g.name as name,"
            "g.residues as residues, l.strand as strand, l.fmin as min,"
            "l.fmax as max").data()
        json.dump(data, open(cache_filename, 'w'))

    if data:
        return data
    else:
        return []
