from Bio.Data import CodonTable
from Bio.SeqUtils import seq3


def compute_annotation_string(interval, var_position, ref, alt, strand, in_gene):
    effect_to_impact = dict(
        stop_lost='HIGH',
        stop_gained='HIGH',
        start_lost='HIGH',
        missense_variant='MODERATE',
        coding_sequence_variant='MODERATE',
        initiator_codon_variant='LOW',
        stop_retained_variant='LOW',
        upstream_gene_variant='MODIFIER',
    )
    codon_table = CodonTable.unambiguous_dna_by_name['Bacterial']

    def aa_name(c):
        return 'Ter' if c in codon_table.stop_codons else seq3(codon_table.forward_table[c])

    gene_info = interval.data
    gene_start = interval.begin
    gene_end = interval.end
    locus_name = gene_info['uniquename']
    gene_name = gene_info['name']
    gene_length = gene_end - gene_start
    # gene_start, gene_end and var_position are 1 based protein_positions
    # the resulting snp_gene_position is a zero-based offset into
    # the residues of the gene
    if strand == 1:
        # plus strand  relative position is past the feature start
        snp_gene_position = var_position - gene_start
    else:
        # minus strand - relative position is before feature end
        snp_gene_position = -(var_position - (gene_end - 1))
    if in_gene:
        dna_change = 'c.' + str(snp_gene_position + 1) + ref + '>' + alt
        # assertion below doesn't work for minus strand, thus removed.
        # assert gene_info['residues'][snp_gene_position] == ref, "Expected {} at {} ({} {}) but found {}".format(ref, snp_gene_position, var_position, strand, gene_info['residues'][snp_gene_position])
        # in gene, coding region variant
        distance_to_feature = ''
        # gene_position string is a 1-based position
        gene_position_string = '{}/{}'.format(snp_gene_position + 1, gene_length)
        snp_aa_position = snp_gene_position // 3
        protein_length = gene_length // 3
        # 1 based coordinates of amino acid changed
        protein_position_string = '{}/{}'.format(snp_aa_position + 1, protein_length)
        snp_codon_start = snp_aa_position * 3
        snp_ref_codon = gene_info['residues'][snp_codon_start:snp_codon_start + 3]
        residue_list = list(gene_info['residues'])
        residue_list[snp_gene_position] = alt  # insert the point mutation
        gene_length = gene_end - gene_start
        mutated_residues = ''.join(residue_list)
        snp_alt_codon = mutated_residues[snp_codon_start:snp_codon_start + 3]
        # print(strand, var_position, snp_ref_codon, snp_alt_codon, ref, alt, snp_codon_start, snp_gene_position, residue_list[snp_gene_position])
        # start codons start at position 0, other codon sequences look like start codons but aren't
        ref_start_codon = snp_codon_start == 0 and snp_ref_codon in codon_table.start_codons
        ref_stop_codon = snp_ref_codon in codon_table.stop_codons
        assert not (ref_stop_codon and ref_start_codon)
        alt_start_codon = snp_codon_start == 0 and snp_alt_codon in codon_table.start_codons
        alt_stop_codon = snp_alt_codon in codon_table.stop_codons
        assert not (alt_stop_codon and alt_start_codon)
        snp_ref_aa = aa_name(snp_ref_codon)
        snp_alt_aa = aa_name(snp_alt_codon)

        if ref_start_codon and not alt_start_codon:
            effect = 'start_lost'
        elif ref_stop_codon and not alt_stop_codon:
            effect = 'stop_lost'
        # "start" codon in middle of sequence is just a normal AA, so
        # 'start_gained' is not worth noting
        # elif not ref_start_codon and alt_start_codon:
        #     effect = 'start_gained'
        elif not ref_stop_codon and alt_stop_codon:
            effect = 'stop_gained'
        elif ref_start_codon and alt_start_codon:
            effect = 'initiator_codon_variant'
        elif ref_stop_codon and alt_stop_codon:
            effect = 'stop_retained_variant'
        elif snp_ref_aa != snp_alt_aa:
            effect = 'missense_variant'
        else:
            effect = 'coding_sequence_variant'
        assert effect in effect_to_impact, "Effect {} has unknown impact".format(effect)
        impact = effect_to_impact[effect]
        protein_change = 'p.(' + snp_ref_aa + str(snp_aa_position + 1) + snp_alt_aa + ')'
    else:
        # promoter region change
        dna_change = 'c.' + str(snp_gene_position) + ref + '>' + alt
        protein_change = ''
        gene_position_string = ''
        protein_position_string = ''
        effect = 'upstream_gene_variant'
        impact = effect_to_impact[effect]
        distance_to_feature = str(abs(snp_gene_position))
    # for SnpEff annotation format see:
    # http://snpeff.sourceforge.net/VCFannotationformat_v1.0.pdf
    annotation = 'ANN=' + '|'.join([alt, effect, impact, gene_name,
                                   locus_name, 'transcript', '',  # TODO: link to protein ID
                                   'Coding', '1/1', dna_change,
                                    protein_change, gene_position_string,
                                    gene_position_string, protein_position_string,
                                    distance_to_feature, ''])
    return annotation
