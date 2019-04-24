import logging
from collections import defaultdict

from variant_merging_constants import VCFVariant
import seq_utils
import pickle
import os
import pandas as pd

def calculate_edited_seq(vcf_var, seq_provider):
    seq, seq_start = seq_provider.get_seq_with_start(vcf_var.chr, vcf_var.pos)

    pos_seq = int(vcf_var.pos) - 1 - seq_start

    assert pos_seq >= 0, "position is below the reference for %s".format(vcf_var)
    assert pos_seq + len(vcf_var.ref) <= len(seq), "position is above the reference for %s".format(vcf_var)

    #assert seq[pos_seq:].startswith(vcf_var.ref)

    edited = ''.join([seq[0:pos_seq], vcf_var.alt, seq[pos_seq + len(vcf_var.ref):]])

    return vcf_var.chr, seq_start, edited


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

    #pickle.dump(variants_dict, open('/files/data/variant_dict.pkl', 'wb'))

    chunk_seq = seq_provider
    #chunk_seq = seq_utils.ChunkBasedSeqProvider(variants_dict.values(), margin, seq_provider)
    #return [{v} for v in variants_dict.keys()]

    variant_eq = [(v_name, calculate_edited_seq(v_rec, chunk_seq)) for
                  v_name, v_rec in variants_dict.items()]

    # dictionary from edited references to a list of variant names
    edited_dict = defaultdict(list)
    for v_name, edited in variant_eq:
        edited_dict[edited].append(v_name)

    # list of sets
    equivalent_variants = []

    for var_lst in edited_dict.values():
        equivalent_variants.append(frozenset({vn for vn in var_lst}))

    #print(len(equivalent_variants))
    print("non merged variants!")
    sing =[eq for eq in equivalent_variants if len(eq) != 2]
    print(len(sing), sing[0:10])

    #print("equiuvalent variants", len(equivalent_variants), "out of", len(variants_dict))
    return equivalent_variants


# for testing purposes
def variant_equal(v1, v2, ref_id, seq_provider):
    assert ref_id == 'hg38'
    v1_norm = calculate_edited_seq(VCFVariant(int(v1[0]), v1[1], v1[2], v1[3]), seq_provider)
    v2_norm = calculate_edited_seq(VCFVariant(int(v2[0]), v2[1], v2[2], v2[3]), seq_provider)

    return v1_norm == v2_norm


def trim_bases(chr, pos, ref, alt):
    #ref, alt = v.split(":")[2].split(">")
    if len(ref) <= 1 or len(alt) <= 1:
        return (chr, pos, ref, alt)
    else:
        (ref, alt) = trim_trailing(ref, alt)
        #v = ":".join(v.split(":")[0:2] + ["{0}>{1}".format(ref, alt)])
        (chr, pos, ref, alt) = trim_leading(chr, pos, ref, alt)
        return (chr, pos, ref, alt)


def trim_trailing(ref, alt):
    if len(ref) <= 1 or len(alt) <= 1:
        return ref, alt
    elif ref[-1] != alt[-1]:
        return ref, alt
    else:
        ref = ref[:-1]
        alt = alt[:-1]
        return trim_trailing(ref, alt)


def trim_leading(chr, pos, ref, alt):
    pos = int(pos)
    if len(ref) == 1 or len(alt) == 1 or ref[0] != alt[0]:
        return (chr, pos, ref, alt)
    else:
        ref = ref[1:]
        alt = alt[1:]
        return trim_leading(chr, str(pos+1), ref, alt)

def main():
    #var_dict = pickle.load(open('/Users/marc/brca/nobackup/enigma_wdir/brca_out/variant_dict.pkl', 'rb'))
    #find_equivalent_variant(var_dict)
    #seq_provider = seq_utils.ChunkBasedSeqProvider(var_dict.values(), 3,
    #                                               '/Users/marc/brca/nobackup/enigma_wdir/resources/seq_repo/latest')

    # 'chr13:g.32340603:T>C'
    # 'chr13:g.32376822:T>C'
    #print(calculate_edited_seq(VCFVariant(13, 32336222, 'T', 'C'), seq_provider))
    #calculate_edited_seq()

    df_eq = pd.read_csv('/Users/marc/brca/nobackup/enigma_wdir/equal_enigma_variants_20181005.txt', header=None, sep='\t')
    df_cases = (pd.concat(
        [df_eq.reset_index().drop(columns=1).rename(columns={0: "key"}),
         df_eq.reset_index().loc[:, ["index", 1]].rename(columns={1: "key"})],
        axis=0, ignore_index=True)
                .rename(columns={'index': 'pairid'}))

    #print(df_cases.head())

    #vd = {v : VCFVariant._make(
    #    (int(v.split(':')[0].lstrip('chr')),
    #    int(v.split(':')[1].lstrip('g.')),
    #    v.split(':')[2].split('>')[0],
    #    v.split(':')[2].split('>')[1])) for v in df_cases['key']}

    vd = {v : VCFVariant._make(trim_bases(
        int(v.split(':')[0].lstrip('chr')),
        int(v.split(':')[1].lstrip('g.')),
        v.split(':')[2].split('>')[0],
        v.split(':')[2].split('>')[1])) for v in df_cases['key']}

    #vd = { v: }

    df = pd.read_csv('/Users/marc/git/brca-exchange/pipeline/luigi/gene_config_brca_only.txt', sep=',', header=0)
    gene_config_df =  df.set_index('symbol', drop=False)

    gene_regions = [seq_utils.ChrInterval(a[0], a[1], a[2]) for a in
     gene_config_df.loc[:, ['chr', 'start_hg38', 'end_hg38']].values]

    seq_wrapper = seq_utils.SeqRepoWrapper('/Users/marc/brca/nobackup/enigma_wdir/resources/seq_repo/latest', regions_preload=gene_regions)

    chunk_seq_provider = seq_utils.ChunkBasedSeqProvider(vd.values(), 50, seq_wrapper)


    find_equivalent_variant(vd, chunk_seq_provider)

    print(calculate_edited_seq(VCFVariant(13, 32356477, 'TAAG', 'T'), chunk_seq_provider))

    print(calculate_edited_seq(VCFVariant(13, 32356482, 'AGAA', 'A'), chunk_seq_provider))

if __name__ == '__main__':
    main()