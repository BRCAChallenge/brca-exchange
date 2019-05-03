import logging
from collections import defaultdict

from variant_merging_constants import VCFVariant


def calculate_edited_seq(vcf_var, seq_provider):
    '''
    Applies sequence edit of a variant to a reference sequence fragment.

    :param vcf_var: VCFVariant record
    :param seq_provider: seq_provider instance returning a sequence fragment along with an offset
    :return: Tuple[int, int, str]: chromosome, offset (wrt chromosome), edited string fragment starting at offset
    '''
    seq, seq_start = seq_provider.get_seq_with_start(vcf_var.chr, vcf_var.pos)

    pos_seq = int(vcf_var.pos) - 1 - seq_start

    assert pos_seq >= 0, "position is below the reference for {}".format(vcf_var)
    assert pos_seq + len(vcf_var.ref) <= len(seq), "position is above the reference for {}".format(vcf_var)

    assert seq[pos_seq:].startswith(vcf_var.ref)

    edited = ''.join([seq[0:pos_seq], vcf_var.alt, seq[pos_seq + len(vcf_var.ref):]])

    return vcf_var.chr, seq_start, edited


def find_equivalent_variants_whole_seq(variants_dict, whole_seq_provider):
    '''
    Determines equivalent variants by editing the reference according to pos,
    ref, alt in the VCF ROW and comparing the resulting strings.

    Edit happens on the entire reference sequence of a gene.
    In order for not to have to keep the modified reference strings for all
    variants in memory, the modified string is hashed and the hashed values are
    used for comparisons. Since two distinct edited strings may in principle have
    the same hash (very unlikely though), some extra check is performed by
    recomputing and comparing the full strings.

    :param variants_dict: dictionary from variant (VCF String, e.g chr13:g.32326103:C>G) to its corresponding VCF row
    :param whole_seq_provider SeqProvider instance obtain sequence information
    :return: list of sets of equivalent variants represented as VCF string
    '''

    logging.info("Running find_equivalent_variants using whole seqs")

    variant_eq = [
        (v_name, hash(calculate_edited_seq(v_rec, whole_seq_provider))) for
        v_name, v_rec in variants_dict.items()]

    # dictionary from hashed edited references to a list of variant names
    hash_dict = defaultdict(list)
    for v_name, edited_hash in variant_eq:
        hash_dict[edited_hash].append(v_name)

    # list of sets
    equivalent_variants = []

    for var_lst in hash_dict.values():
        # var_lst should contain names of equivalent variants based on the hash
        # edited reference. Doing an extra check with the actual strings
        vd = defaultdict(list)
        for vn in var_lst:
            edited = calculate_edited_seq(variants_dict[vn],
                                          whole_seq_provider)
            vd[edited].append(vn)

        if len(vd) > 1:
            logging.debug(
                "Hash Collisions. Involved variants were {}".format(var_lst))

        for var_lst2 in vd.values():
            equivalent_variants.append(frozenset({vn for vn in var_lst2}))

    return equivalent_variants


def find_equivalent_variant(variants_dict, chunk_seq_provider):
    '''
    Determines equivalent variants by editing the reference according to pos,
    ref, alt in the VCF ROW and comparing the resulting strings.

    The edit doesn't happen on the entire reference sequence of a gene, but only
    on the corresponding 'chunk' of the variant. Chunk sequences are directly compared
    and are not hashed

    :param variants_dict: dictionary from variant (VCF String, e.g chr13:g.32326103:C>G) to its corresponding VCF row
    :param chunk_seq_provider:
    :return: list of sets of equivalent variants represented as VCF string
    '''

    logging.info("Running find_equivalent_variants using chunk strategy")

    variant_eq = [(v_name, calculate_edited_seq(v_rec, chunk_seq_provider)) for
                  v_name, v_rec in variants_dict.items()]

    # dictionary from edited references to a list of variant names
    edited_dict = defaultdict(list)
    for v_name, edited in variant_eq:
        edited_dict[edited].append(v_name)

    # list of sets
    equivalent_variants = []

    for var_lst in edited_dict.values():
        equivalent_variants.append(frozenset({vn for vn in var_lst}))

    return equivalent_variants


# for testing purposes
def variant_equal(v1, v2, ref_id, seq_provider):
    assert ref_id == 'hg38'
    v1_norm = calculate_edited_seq(VCFVariant(int(v1[0]), v1[1], v1[2], v1[3]), seq_provider)
    v2_norm = calculate_edited_seq(VCFVariant(int(v2[0]), v2[1], v2[2], v2[3]), seq_provider)

    return v1_norm == v2_norm
