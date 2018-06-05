from intervaltree import IntervalTree


def make_interval_trees(gene_info):
    """Take gene info from Neo4j COMBAT TB database and insert into Interval Trees for plus and minus strands"""
    plus_tree = IntervalTree()
    minus_tree = IntervalTree()
    for gene in gene_info:
        assert 'strand' in gene, "Gene with no strand encountered: {}".format(gene)
        assert gene['strand'] in (1, -1), "Gene with unknown strand type encountered: {}".format(gene)
        gene_info = dict(uniquename=gene['uniquename'], name=gene['name'], residues=gene['residues'])
        start = int(gene['min']) + 1  # NOTE NOTE NOTE: switch from 0 based to 1 base coordinates
        end = int(gene['max']) + 1
        if gene['strand'] == 1:
            plus_tree[start:end] = gene_info
        else:
            minus_tree[start:end] = gene_info
    return (plus_tree, minus_tree)
