import os.path
import pytest
import snptools.combattb_db
import snptools.itree
from snptools.vcf import compute_annotation_string


@pytest.fixture
def gene_data():
    stored_data = os.path.join(os.path.dirname(__file__), 'combattb_gene_list.cache')
    data = snptools.combattb_db.get_gene_data(cache_filename=stored_data)
    return data


@pytest.fixture
def itrees(gene_data):
    (plus_tree, minus_tree) = snptools.itree.make_interval_trees(gene_data)
    return (plus_tree, minus_tree)


def test_gene_data_keys(gene_data):
    for gene_info in gene_data:
        assert sorted(gene_info.keys()) == ['max', 'min', 'name', 'residues', 'strand', 'uniquename'], "Expected keys: {} but got {}".format(['max', 'min', 'name', 'residues', 'strand', 'uniquename'], sorted(gene_info.keys()))

    num_genes = 4018
    assert len(gene_data) == num_genes, "Expected to find {} genes, found {}".format(num_genes, len(gene_data))


def test_gene_data_item(gene_data):
    for gene_info in gene_data:
        test_locus = 'Rv0009'
        if gene_info['uniquename'] == test_locus:
            assert gene_info['name'] == 'ppiA', "Expected gene name ppiA for {}".format(test_locus)
            # these are zero-based coordinates
            assert gene_info['min'] == 12467, "Expected start position of ppiA to be 12467, got {}".format(gene_info['min'])
            assert gene_info['max'] == 13016, "Expected start position of ppiA to be 13016, got {}".format(gene_info['min'])
            assert gene_info['strand'] == 1, "Expected strand of ppiA to be 1, got {}".format(gene_info['strand'])
            assert gene_info['residues'].startswith('ATGGCA'), "Expected ppiA residues to start with 'ATGGCA', got {}".format(gene_info['residues'][:6])
            assert gene_info['residues'].endswith('TCCTGA'), "Expected ppiA residues to start with 'TTCTGA', got {}".format(gene_info['residues'][-6:])


def test_interval_trees(itrees):
    (plus_tree, minus_tree) = itrees

    # most 3' base of Rv2237A - minus strand
    minus_test_end = 2510351
    # most 5' base of Rv2237A - minus strand
    minus_test_start = 2510587
    for pos in (minus_test_end, minus_test_start):
        assert minus_tree.overlaps(pos), "Expected to find a gene at {} on the - strand.".format(pos)
    assert not minus_tree.overlaps(minus_test_end - 1), "Expected to find no gene at {} on the - strand".format(minus_test_end - 1)
    assert not minus_tree.overlaps(minus_test_start + 1), "Expected to find no gene at {} on the - strand".format(minus_test_start + 1)

    gene_info = list(minus_tree[minus_test_start])[0]
    assert gene_info.data['uniquename'] == 'Rv2237A', "Expected to find Rv2237A at {} on the - strand".format(minus_test_start)

    plus_test_start = 2509489
    plus_test_end = 2510256
    for pos in (plus_test_start, plus_test_end):
        assert plus_tree.overlaps(pos), "Expected to find a gene at {} on the + strand".format(pos)
    assert not plus_tree.overlaps(plus_test_start - 1), "Expected to find no gene at {} on the + strand".format(plus_test_start - 1)
    assert not plus_tree.overlaps(plus_test_end + 1), "Expected to find no gene at {} on the + strand".format(plus_test_end + 1)

    gene_info = list(plus_tree[plus_test_start])[0]
    assert gene_info.data['uniquename'] == 'Rv2237', "Expected to find Rv2237 at {} on the + strand".format(plus_test_start)


def test_annotation_creation(itrees):
    (plus_tree, minus_tree) = itrees
    position = 2509490
    gene_info = list(plus_tree[position])[0]
    # assert gene_info is None, "expected {}".format(gene_info)
    gene_length = gene_info.end - gene_info.begin
    assert gene_length % 3 == 0, "Expected a multiple of 3, got {} from length {}".format(gene_length % 3, gene_length)
    ref = 'T'
    alt = 'C'
    strand = 1
    ann_string = compute_annotation_string(gene_info,
                                           position,
                                           ref, alt,
                                           strand,
                                           True)

    for text in ('start_lost', 'HIGH', 'c.2T>C', 'p.(Met1Thr)', '2/768', '1/256'):
        assert text in ann_string, "Expected {} in annotation string: {}".format(text, ann_string)

    position = 2509489
    gene_info = list(plus_tree[position])[0]
    ref = 'A'
    alt = 'G'
    strand = 1
    ann_string = compute_annotation_string(gene_info,
                                           position,
                                           ref, alt,
                                           strand,
                                           True)

    for text in ('initiator_codon_variant', 'LOW', 'c.1A>G', 'p.(Met1Val)', '1/768', '1/256'):
        assert text in ann_string, "Expected {} in annotation string: {}".format(text, ann_string)

    position = 2509496
    gene_info = list(plus_tree[position])[0]
    ref = 'T'
    alt = 'A'
    strand = 1
    ann_string = compute_annotation_string(gene_info,
                                           position,
                                           ref, alt,
                                           strand,
                                           True)

    for text in ('stop_gained', 'HIGH', 'c.8T>A', 'p.(Leu3Ter)', '8/768', '3/256'):
        assert text in ann_string, "Expected {} in annotation string: {}".format(text, ann_string)

    ref = 'T'
    alt = 'G'
    strand = 1
    ann_string = compute_annotation_string(gene_info,
                                           position,
                                           ref, alt,
                                           strand,
                                           True)

    for text in ('missense_variant', 'MODERATE', 'c.8T>G', 'p.(Leu3Trp)', '8/768', '3/256'):
        assert text in ann_string, "Expected {} in annotation string: {}".format(text, ann_string)

    position = 2510255
    gene_info = list(plus_tree[position])[0]
    ref = 'G'
    alt = 'C'
    strand = 1
    ann_string = compute_annotation_string(gene_info,
                                           position,
                                           ref, alt,
                                           strand,
                                           True)

    for text in ('stop_lost', 'HIGH', 'c.767G>C', 'p.(Ter256Ser)', '767/768', '256/256'):
        assert text in ann_string, "Expected {} in annotation string: {}".format(text, ann_string)

    ref = 'G'
    alt = 'A'
    strand = 1
    ann_string = compute_annotation_string(gene_info,
                                           position,
                                           ref, alt,
                                           strand,
                                           True)

    for text in ('stop_retained_variant', 'LOW', 'c.767G>A', 'p.(Ter256Ter)', '767/768', '256/256'):
        assert text in ann_string, "Expected {} in annotation string: {}".format(text, ann_string)

    position = 2510351
    gene_info = list(minus_tree[position])[0]
    ref = 'G'
    alt = 'C'
    strand = -1
    ann_string = compute_annotation_string(gene_info,
                                           position,
                                           ref, alt,
                                           strand,
                                           True)

    for text in ('stop_lost', 'HIGH', 'c.237G>C', 'p.(Ter79Tyr)', '237/237', '79/79'):
        assert text in ann_string, "Expected {} in annotation string: {}".format(text, ann_string)

    position = 2510587
    gene_info = list(minus_tree[position])[0]
    ref = 'G'
    alt = 'A'
    strand = -1
    ann_string = compute_annotation_string(gene_info,
                                           position,
                                           ref, alt,
                                           strand,
                                           True)

    for text in ('initiator_codon_variant', 'LOW', 'c.1G>A', 'p.(Val1Met)', '1/237', '1/79'):
        assert text in ann_string, "Expected {} in annotation string: {}".format(text, ann_string)

    position = 2510588
    gene_info = sorted(minus_tree[position - 1000:position])[-1]
    ref = 'T'
    alt = 'A'
    strand = -1
    ann_string = compute_annotation_string(gene_info,
                                           position,
                                           ref, alt,
                                           strand,
                                           False)

    for text in ('upstream_gene_variant', 'MODIFIER', 'Rv2237A', 'c.-1T>A'):
        assert text in ann_string, "Expected {} in annotation string: {}".format(text, ann_string)

    position = 2498822
    gene_info = sorted(plus_tree[position:position + 1000])[0]
    ref = 'C'
    alt = 'A'
    strand = 1
    ann_string = compute_annotation_string(gene_info,
                                           position,
                                           ref, alt,
                                           strand,
                                           False)

    for text in ('upstream_gene_variant', 'MODIFIER', 'Rv2226', 'c.-10C>A'):
        assert text in ann_string, "Expected {} in annotation string: {}".format(text, ann_string)
