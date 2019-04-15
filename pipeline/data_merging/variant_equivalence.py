import logging
from collections import defaultdict


def calculate_edited_seq(chr, pos_chr, ref, alt, seq_provider):
    seq, seq_start = seq_provider.get_seq_with_start(chr)

    pos_seq = int(pos_chr) - seq_start - 1

    assert pos_seq >= 0, "v1 positions is below the reference"
    assert pos_seq + len(ref) <= len(seq), "v1 position is above the reference"

    assert seq[pos_seq:].startswith(ref)

    edited = ''.join([seq[0:pos_seq], alt, seq[pos_seq + len(ref):]])

    return edited


def calculate_edited_seq_from_rec(v_rec, seq_provider):
    return calculate_edited_seq(int(v_rec[COLUMN_VCF_CHR]),
                                int(v_rec[COLUMN_VCF_POS]),
                                v_rec[COLUMN_VCF_REF],
                                v_rec[COLUMN_VCF_ALT], seq_provider)


from variant_merging_constants import COLUMN_VCF_CHR, COLUMN_VCF_POS, \
    COLUMN_VCF_REF, COLUMN_VCF_ALT


def find_equivalent_variant(variants_dict, seq_provider):
    '''
    Determines equivalent variants by editing the reference according to pos,
    ref, alt in the VCF ROW and comparing the resulting strings.

    In order for not to have to keep the modified reference strings for all
    variants in memory, the modified string is hashed and the hashed values are
    used for comparisons. Since two distinct edited strings may in principle have
    the same hash (very unlikely though), some extra check is performed by
    recomputing and comparing the full strings.

    :param variants_dict: dictionary from variant (VCF String, e.g chr13:g.32326103:C>G) to its corresponding VCF row
    :param seq_provider SeqProvider instance obtain sequence information
    :return: list of sets of equivalent variants represented as VCF string
    '''
    logging.info("Running find_equivalent_variants.")

    variant_eq = [(v_name, hash(calculate_edited_seq_from_rec(v_rec, seq_provider))) for
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
            edited = calculate_edited_seq_from_rec(variants_dict[vn],
                                                   seq_provider)
            vd[edited].append(vn)

        if len(vd) > 1:
            logging.debug(
                "Hash Collisions. Involved variants were {}".format(var_lst))

        for var_lst2 in vd.values():
            equivalent_variants.append(frozenset({vn for vn in var_lst2}))

    return equivalent_variants


# for testing purposes
def variant_equal(v1, v2, ref_id, seq_provider):
    assert ref_id == 'hg38'
    v1_norm = calculate_edited_seq(int(v1[0]), v1[1], v1[2], v1[3], seq_provider)
    v2_norm = calculate_edited_seq(int(v2[0]), v2[1], v2[2], v2[3], seq_provider)

    return v1_norm == v2_norm
