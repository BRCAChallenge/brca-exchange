import os
import logging
import hashlib

def brca_normalize_chunks(chr1, pos, ref1, alt1, seq_lookup):
    # chr1, pos1, ref1, alt1 = v1
    pos1 = int(pos)

    seq_dict = seq_lookup[chr1]

    seq = seq_dict["sequence"]
    pos1 = int(pos) - seq_dict["start"] - 1

    assert pos1 >= 0, "v1 positions is below the reference"
    reflen = len(seq)
    assert pos1 + len(ref1) <= reflen, "v1 position is above the reference"

    # assert seq[pos1:].startswith(ref1)
    edited = seq[0:pos1] + alt1 + seq[pos1 + len(ref1):]

    def compute_delta_joint(ref, alt):
        dref_left = 0
        dref_right = len(ref) - 1
        delta_left = 0
        dalt_left = 0
        dalt_right = len(alt) - 1
        delta_right = 0

        while dref_right >= dref_left and dalt_right >= dalt_left:
            # remove on the right first to obtain left aligned changes
            if ref[dref_right] == alt[dalt_right]:
                dref_right -= 1
                dalt_right -= 1
                delta_right += 1
            elif ref[dref_left] == alt[dalt_left]:
                dref_left += 1
                dalt_left += 1
                delta_left += 1
            else:
                break

        return delta_left, delta_right, ref[dref_left:dref_right + 1], alt[
                                                                       dalt_left:dalt_right + 1]

    delta_left, delta_right, ref_adj, alt_adj = compute_delta_joint(ref1,
                                                                    alt1)
    pos1 = pos1 + delta_left

    repeat_delta = 0

    # TODO: need to rule out trivial SNPs, like A>A
    # TODO: function
    if len(ref_adj) == 0:
        # insertion
        dt = len(alt_adj)
        if len(set(alt_adj)) == 1:
            dt = 1

        while (seq[(pos1 - repeat_delta - dt):(
                (pos1 - repeat_delta - dt + len(alt_adj)))] == alt_adj):
            repeat_delta += dt
    elif len(alt_adj) == 0:
        # deletion
        dt = len(ref_adj)
        if len(set(ref_adj)) == 1:
            dt = 1
        while (seq[(pos1 - repeat_delta - dt):(
                (pos1 - repeat_delta - dt + len(ref_adj)))] == ref_adj):
            repeat_delta += dt

    pos_repeat_adj = pos1 - repeat_delta

    return hashlib.sha1(edited).digest(), (
    pos_repeat_adj, alt_adj, pos_repeat_adj + len(ref_adj))


def get_seq_lookup():
    pwd = os.path.dirname(os.path.realpath(__file__))
    reference = os.path.join(pwd, '..', 'data') + '/'

    return {17: {"start": 43000000,
                 "sequence": open(reference + "brca1_hg38.txt",
                                  "r").read().upper()},
            13: {"start": 32300000,
                 "sequence": open(reference + "brca2_hg38.txt",
                                  "r").read().upper()}
            }


seq_lookup = get_seq_lookup()

from variant_merging_constants import COLUMN_VCF_CHR, COLUMN_VCF_POS, \
    COLUMN_VCF_REF, COLUMN_VCF_ALT


def find_equivalent_variant(variants):
    logging.info("Running find_equivalent_variants.")

    # return list of set of equivlanet variants

    variant_eq = [(vn, brca_normalize_chunks(int(v[COLUMN_VCF_CHR]),
                                             int(v[COLUMN_VCF_POS]),
                                             v[COLUMN_VCF_REF],
                                             v[COLUMN_VCF_ALT], seq_lookup)) for
                  vn, v in variants.items()]

    hash_dict = {}
    for vn, diff_req in variant_eq:
        h = diff_req[0]  # hash

        if h in hash_dict:
            l = hash_dict[h]
        else:
            l = list()
            hash_dict[h] = l
        l.append((vn, diff_req))

    # list of sets
    equivalent_variants = []
    for h, ldiff_req in hash_dict.items():
        print(ldiff_req)
        eqs = {d for _, (_, d) in ldiff_req}
        assert len(eqs) == 1, str(eqs)

        equivalent_variants.append({vn for vn, _ in ldiff_req})

    return equivalent_variants


# for testing purposes
def variant_equal(v1, v2, ref_id):
    assert ref_id == 'hg38'
    v1_norm = brca_normalize_chunks(int(v1[0]), v1[1], v1[2], v1[3], seq_lookup)
    v2_norm = brca_normalize_chunks(int(v2[0]), v2[1], v2[2], v2[3], seq_lookup)

    return v1_norm[0] == v2_norm[0]
