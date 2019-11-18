import itertools
import glob
import os
import unittest

import bioutils.seqfetcher
from mock import patch
import pytest
import vcf
from hypothesis import given, assume, settings
from hypothesis.strategies import integers, tuples, sampled_from, lists

from common import seq_utils
from common.config import load_config, extract_gene_regions_dict
from utilities import round_sigfigs
from variant_equivalence import variant_equal, find_equivalent_variant, find_equivalent_variants_whole_seq
from variant_merging import normalize_values, add_variant_to_dict, \
    COLUMN_SOURCE, append_exac_allele_frequencies, EXAC_SUBPOPULATIONS


from variant_merging_constants import VCFVariant

runtimes = 500000
settings.register_profile('ci', settings(max_examples=runtimes, max_iterations=runtimes, timeout=-1))

# Uncomment this for longer test runs.
#settings.load_profile('ci')

#
# initialize module
#

pwd = os.path.dirname(os.path.realpath(__file__))
data_dir = os.path.join(pwd, '..', 'data')
seq_provider = seq_utils.LegacyFileBasedSeqProvider(data_dir)

BRCA1_CHR = 17
reference_length = len(seq_provider.get_seq_with_start(BRCA1_CHR).sequence)

# variants are specfied in genome coords. BRCA1 and BRCA2 hold slices of the
# reference genomes, starting before the gene. Currently, all the reference
# slices have the same length. To test, we need to ensure that the variant
# falls within the region of the slice. If we pick a number in  [0, N] where N
# is the length of the reference, then add the start position of the reference,
# we have a genomic coordinate in the reference. Add one to have a 1-based
# coord.


#
# generators
#
default_reference_version = "hg38"
reference_id = sampled_from([default_reference_version])
chrom = sampled_from(('13', '17'))
pos = integers(min_value=1, max_value=reference_length) # pos relative to the reference start.
base = sampled_from(('A', 'C', 'T', 'G'))
subseq = lists(base).map(lambda x: ''.join(x))

# Note this doesn't catch all degenerate variants. A -> A would pass, for example.
def not_noop(v):
    "False if ref and alt are empty"
    (_, _, ref, alt) = v
    return ref != '' or alt != ''

# A variant with all random fields.
# (chrom, position, ref, alt)
variant = tuples(chrom, pos, subseq, subseq).filter(not_noop)


# Note this doesn't catch all degenerate variants. A -> A would pass, for example.
def not_on_ref_noop(v):
    "False if ref and alt are empty"
    (_, _, reflen, alt) = v
    return reflen != 0 or alt != ''

# A variant where 'ref' is a length. We will later copy in the reference base
# pairs from a reference sequence.
# (chrom, position, ref, alt)
variant_on_ref = tuples(chrom, pos, integers(min_value=0), subseq).filter(not_on_ref_noop)

chrom_ref = {
    '13': { default_reference_version : seq_provider.get_seq_with_start(13)._asdict() },
    '17': { default_reference_version :  seq_provider.get_seq_with_start(17)._asdict() }
}

# Add the start position of the slice of reference that we're holding, to create
# a valid chrom position from an integer in [0, len(reference)]. Also converts
# to 1-based coords.
def add_start(v, ref_id):
    "Add reference start position to a variant"
    (chrom, pos, ref, alt) = v

    return chrom, \
           pos + chrom_ref[chrom][ref_id]['start'] + 1, \
           chrom_ref[chrom][ref_id]['sequence'][pos: (pos + len(ref))], alt


#
# tests
#

def test_variant_equal_throws_below_reference():
    ref_id = chrom_ref['13'].keys()[0]
    v1 = ('13', chrom_ref['13'][ref_id]['start'], 'A', 'C')
    v2 = ('13', chrom_ref['13'][ref_id]['start'] + 10, 'A', 'C')
    with pytest.raises(AssertionError):
        variant_equal(v1, v2, ref_id, seq_provider)
    with pytest.raises(AssertionError):
        variant_equal(v2, v1, ref_id, seq_provider)

def test_variant_equal_throws_above_reference():
    ref_id = chrom_ref['13'].keys()[0]
    start = chrom_ref['13'][ref_id]['start']
    v1 = ('13', start + reference_length + 2, 'A', 'C')
    v2 = ('13', start + 10, 'A', 'C')
    with pytest.raises(AssertionError):
        variant_equal(v1, v2, ref_id, seq_provider)
    with pytest.raises(AssertionError):
        variant_equal(v2, v1, ref_id, seq_provider)

#@settings(max_examples=50000)
@given(variant, reference_id)
def test_variant_equal_identity(v, ref_id):
    "A variant should be equal to itself"
    (_, pos, ref, _) = v
    assume(pos + len(ref) <= reference_length)
    v = add_start(v, ref_id)
    assert variant_equal(v, v, ref_id, seq_provider)

#@settings(max_examples=50000)
@given(variant, variant, reference_id)
def test_variant_equal_commutative(v1, v2, ref_id):
    "Comparing x, y should be the same as comparing y, x"
    (_, pos1, ref1, _) = v1
    (_, pos2, ref2, _) = v2
    assume(pos1 + len(ref1) <= reference_length)
    assume(pos2 + len(ref2) <= reference_length)

    v1 = add_start(v1, ref_id)
    v2 = add_start(v2, ref_id)
    assert variant_equal(v1, v2, ref_id, seq_provider) == variant_equal(v2, v1, ref_id, seq_provider)


# The rules here describe how to generate equivalent variants.  I believe these
# cover all the equivalencies, however I haven't rigorously proved these
# assertions, so they may be incorrect.
#
# Equivalent deletes may be found by shifting the variant left or right by 1.
# The delete to the right is equivalent if the 1st base in the delete is equal
# to the 1st base after the delete. The delete to the left is equivalent if the
# 1st base before the delete is equal to the last base in the delete.
#
# Equivalent inserts may be found by computing the inverse operation (delete
# the insert), finding all the equivalent deletes, and inverting those deletes
# (to find the corresponding inserts).
#
# There are no equivalents for substitutions, or substitutions combined with an
# insert or delete.


# case insensitive string equality
def strieq(x, y):
    return x.upper() == y.upper()

def right_equiv_deletes(ref, pos, length):
    while pos < len(ref) - length and strieq(ref[pos], ref[pos + length]):
        yield pos + 1
        pos += 1

def left_equiv_deletes(ref, pos, length):
    while pos > 0 and strieq(ref[pos - 1], ref[pos + length - 1]):
        yield pos - 1
        pos -= 1

def equiv_deletes(ref, pos, length):
    return itertools.chain(left_equiv_deletes(ref, pos, length), right_equiv_deletes(ref, pos, length))

# normalize variant
def normalize_variant(variant):
    (chrom, pos, ref, alt) = variant
    while len(ref) > 0 and len(alt) > 0 and ref[0] == alt[0]:
        pos += 1
        ref = ref[1:]
        alt = alt[1:]
    while len(ref) > 0 and len(alt) > 0 and ref[-1] == alt[-1]:
        ref = ref[:-1]
        alt = alt[:-1]
    return (chrom, pos, ref, alt)

def equiv_variants_delete(reference, variant):
    (chrom, pos, ref, alt) = normalize_variant(variant)
    assert(alt == '') # must be a simple delete
    l = len(ref)
    equiv_pos = equiv_deletes(reference, pos, l)
    return [(chrom, epos, reference[epos : epos + l], alt) for epos in equiv_pos]

#
# We generate variants in zero-based coords, so an insert
# CT => CGGT is 1, '', 'GG'

def right_equiv_inserts(ref, pos, alt):
    while pos < len(ref) and strieq(ref[pos], alt[0]):
        alt = alt[1:] + alt[0]
        pos += 1
        yield (pos, alt)

def left_equiv_inserts(ref, pos, alt):
    while pos > 0 and strieq(ref[pos - 1], alt[-1]):
        alt = alt[-1] + alt[:-1]
        pos -= 1
        yield (pos, alt)

def equiv_inserts(ref, pos, alt):
    return itertools.chain(left_equiv_inserts(ref, pos, alt), right_equiv_inserts(ref, pos, alt))

def equiv_variants_insert(reference, variant):
    (chrom, pos, ref, alt) = normalize_variant(variant)
    assert(ref == '') # must be a simple insert
    equiv_pos = equiv_inserts(reference, pos, alt)
    return [(chrom, epos, ref, ealt) for (epos, ealt) in equiv_pos]

#
#

def is_deletion(v):
    (_, _, ref, alt) = v
    return len(ref) > 0 and len(alt) == 0

def is_insertion(v):
    (_, _, ref, alt) = v
    return len(ref) == 0 and len(alt) > 0

def is_in_bounds(v):
    (chrom, pos, ref, alt) = v
    return pos + len(ref) < reference_length

# Find all equivalent variant, neglecting changes
# in surrounding context.
def all_norm_equiv(refsequence, v):
    vnorm = normalize_variant(v)

    equivs = []
    if is_deletion(vnorm):
        equivs = equiv_variants_delete(refsequence, vnorm)
    elif is_insertion(vnorm):
        equivs = equiv_variants_insert(refsequence, vnorm)

    return equivs

# Take a variant with position and ref length, and copy
# in the ref bases from a reference sequence.
def inject_ref(refsequence, v_on_r):
    (chrom, pos, reflen, alt) = v_on_r
    return (chrom, pos, refsequence[pos : pos + reflen], alt)


# To test that equivalent variants test equal, generate a random variant, then
# compute variants equivalent to it, and assert that they all test equal to the
# first variant. Currently on works for deletes.


#@settings(max_examples=runtimes, max_iterations=runtimes, timeout=-1, database_file=None)
#@settings(max_examples=runtimes, max_iterations=runtimes, timeout=-1)
@given(variant_on_ref, reference_id)
def test_variant_equal_equiv(v, ref_id):
    (chrom, pos, reflen, alt) = v
    refsequence = chrom_ref[chrom][ref_id]["sequence"]
    assume(pos + reflen <= len(refsequence))
    v = inject_ref(refsequence, v)

    equivs = all_norm_equiv(refsequence, v)

    for veq in equivs:
        if is_in_bounds(veq):
            assert variant_equal(add_start(v, ref_id), add_start(veq, ref_id), ref_id, seq_provider)

def equiv_set(refsequence, v):
    return set(all_norm_equiv(refsequence, v) + [normalize_variant(v)])

#@settings(max_examples=runtimes, max_iterations=runtimes, timeout=-1, database_file=None)
#@settings(max_examples=runtimes, max_iterations=runtimes, timeout=-1)
@given(variant_on_ref, variant_on_ref, reference_id)
def test_variant_equal_not_equiv(v1, v2, ref_id):
    (chrom1, pos1, reflen1, alt1) = v1
    (chrom2, pos2, reflen2, alt2) = v2
    assume(pos1 + reflen1 <= reference_length)
    assume(pos2 + reflen2 <= reference_length)
    assume(not alt1 == '' and reflen1 == 0)
    assume(not alt2 == '' and reflen2 == 0)

    refsequence1 = chrom_ref[chrom1][ref_id]["sequence"]
    refsequence2 = chrom_ref[chrom2][ref_id]["sequence"]
    v1 = inject_ref(refsequence1, v1)
    v2 = inject_ref(refsequence2, v2)

    eq1 = equiv_set(refsequence1, v1)
    eq2 = equiv_set(refsequence2, v2)

    if len(eq1.intersection(eq2)) == 0: # should test not-equal
        assert not variant_equal(add_start(v1, ref_id), add_start(v2, ref_id), ref_id, seq_provider)
    else:
        assert variant_equal(add_start(v1, ref_id), add_start(v2, ref_id), ref_id, seq_provider)

# Do we need to explicitly test variations in surrounding reference length?
# The tests above only test random variants against normalized (minimum reference)
# variants.


class TestVariantMerging(unittest.TestCase):

    def setUp(self):
        self.variant_dict = {'chr13:g.32339228:GAA>G':
                             ['ENIGMA', 'BRCA2', 'chr13:32339228:GAA>G', '13', '32339228', 'GAA', 'G', 'NM_000059.3',
                              'c.4876_4877delAA ', '5104delAA', 'N1626Sfs*12', '', 'OMIM',
                              'BREAST-OVARIAN CANCER, FAMILIAL, SUSCEPTIBILITY TO, 2; BROVCA2 (612555)', 'Disease',
                              'Pathogenic', '22/4/2016', 'ENIGMA BRCA1/2 Classification Criteria (2015)',
                              'https://enigmaconsortium.org/wp-content/uploads/2016/06/ENIGMA_Rules_2015-03-26.pdf', '',
                              'Variant allele predicted to encode a truncated non-functional protein.', 'Curation',
                              'Germline', 'SCV000282396.1', 'p.(Asn1626SerfsTer12)', '46']
                             }
        self.genomic_coordinate = 'chr13:g.32339228:GAA>G'
        self.values_to_add = ['BIC', 'BRCA2', 'chr13:32339228:GAA>G', '13', '32339228', 'GAA', 'G', 'NM_000059.3',
                              'c.4876_4877delAA', '', 'N1626Sfs*12', '', 'OMIM',
                              'BREAST-OVARIAN CANCER, FAMILIAL, SUSCEPTIBILITY TO, 2; BROVCA2 (612555)', 'Disease',
                              'Pathogenic', '2016-09-08', 'ENIGMA BRCA1/2 Classification Criteria (2015)',
                              'https://enigmaconsortium.org/wp-content/uploads/2016/06/ENIGMA_Rules_2015-03-26.pdf', '',
                              'Variant allele predicted to encode a truncated non-functional protein.', 'Curation',
                              'Germline', 'SCV000282396.1', 'p.(Asn1626SerfsTer12)', '677']

    def test_normalize_values(self):
        empty_string = normalize_values('')
        none_value = normalize_values(None)
        whitespace = normalize_values(' value ')
        list_values = normalize_values(['dog ', ' fish ', '', None, 'dog'])
        list_values_default = normalize_values(['-'])
        list_values_empty_list = normalize_values([])
        self.assertEqual(empty_string, '-')
        self.assertEqual(none_value, '-')
        self.assertEqual(whitespace, 'value')
        self.assertEqual(list_values, ['dog', 'fish'])
        self.assertEqual(list_values_default, ['-'])
        self.assertEqual(list_values_empty_list, ['-'])

    def test_add_variant_to_dict(self):
        genomic_coordinate = 'chr13:g.32332705:GA>G'
        values = ['ENIGMA', 'BRCA2', 'chr13:32332705:GA>G', '13', '32332705', 'GA', 'G', 'NM_000059.3', 'c.1231delA',
                  '1459delA', 'I411Yfs*19', '', 'OMIM',
                  'BREAST-OVARIAN CANCER, FAMILIAL, SUSCEPTIBILITY TO, 2; BROVCA2 (612555)', 'Disease', 'Pathogenic',
                  '22/4/2016', 'ENIGMA BRCA1/2 Classification Criteria (2015)',
                  'https://enigmaconsortium.org/wp-content/uploads/2016/06/ENIGMA_Rules_2015-03-26.pdf', '',
                  'Variant allele predicted to encode a truncated non-functional protein.', 'Curation', 'Germline',
                  'SCV000282353.1', 'p.(Ile411TyrfsTer19)', '2']

        variant_dict = add_variant_to_dict(self.variant_dict, genomic_coordinate, values)
        self.assertIn('chr13:g.32332705:GA>G', variant_dict)
        self.assertEqual(variant_dict['chr13:g.32332705:GA>G'], values)

    def test_add_variant_to_dict_merge_different_values(self):
        variant_dict = add_variant_to_dict(self.variant_dict, self.genomic_coordinate, self.values_to_add)
        merged = variant_dict[self.genomic_coordinate]
        date_index = 16
        bx_id_index = -1
        self.assertIn('22/4/2016', merged[date_index])
        self.assertIn('2016-09-08', merged[date_index])
        self.assertIn('677', merged[bx_id_index])
        self.assertIn('46', merged[bx_id_index])
        self.assertIn('BIC', merged[COLUMN_SOURCE])
        self.assertIn('ENIGMA', merged[COLUMN_SOURCE])

    def test_add_variant_to_dict_merge_ignores_trailing_spaces(self):
        variant_dict = add_variant_to_dict(self.variant_dict, self.genomic_coordinate, self.values_to_add)
        merged = variant_dict[self.genomic_coordinate]
        self.assertEqual('c.4876_4877delAA', merged[8])

    def test_add_variant_to_dict_merge_adds_new_data_to_empty_fields(self):
        variant_dict = add_variant_to_dict(self.variant_dict, self.genomic_coordinate, self.values_to_add)
        merged = variant_dict[self.genomic_coordinate]
        self.assertEqual('5104delAA', merged[9])

    def test_append_exac_allele_frequencies_rounds_to_three_sig_figs(self):
        EXAC_VCF_FILENAME = os.path.join(os.path.dirname(__file__), 'test_files/ExAC_AF.vcf')
        for record in vcf.Reader(open(EXAC_VCF_FILENAME, 'r')):
            record = append_exac_allele_frequencies(record, new_record=None, i=None)
            for subpopulation in EXAC_SUBPOPULATIONS:
                val = record.INFO["AF_" + subpopulation]
                try:
                    float_val = float(val)
                    self.assertEqual(float_val, round_sigfigs(float(val), 3))
                except ValueError:
                    self.assertEqual(val, '-')

@pytest.fixture()
def fetch_seq_mock_data():
    mock_data = {}

    # mock data was generated by running tests using 'make test' and
    # saving the output of SeqRepoWrapper._fetch_seq into the corresponding files
    for path in glob.glob(os.path.join(data_dir, 'mock_fetch_seq*')):
        with open(path, 'r') as f:
            fn = os.path.basename(path).rstrip('.txt')
            key = tuple(fn.split('-')[-3:]) # ac, start, end
            mock_data[key] = f.readline().strip()

    return mock_data


def test_find_equivalent_variant(fetch_seq_mock_data):
    # mocking _fetch_seq method
    with patch.object(bioutils.seqfetcher, 'fetch_seq', side_effect=lambda ac, s, e: fetch_seq_mock_data[(str(ac), str(s), str(e))]):
        seq_wrapper = seq_utils.SeqRepoWrapper()
        # empty case
        assert [] == find_equivalent_variant({}, seq_wrapper)

        # a bunch of variants. If they appear in the same set, they are considered equivalent
        example_variants = [
            frozenset({'chr13:g.32355030:A>AA'}),
            frozenset({'chr13:g.32339774:GAT>G', 'chr13:g.32339776:TAT>T'}),
            frozenset({'chr17:g.43090921:G>GCA', 'chr17:g.43090921:GCA>GCACA'})
        ]

        # construct variant dict (flattening example_variants!)
        variant_dict = {v: VCFVariant(
            int(v.split(':')[0].lstrip('chr')),
            int(v.split(':')[1].lstrip('g.')),
            v.split(':')[2].split('>')[0],
            v.split(':')[2].split('>')[1]) for eq_variants in example_variants for v
        in
            eq_variants
        }

        margin = 20
        chunk_provider = seq_utils.ChunkBasedSeqProvider(variant_dict.values(), margin, seq_wrapper)

        assert frozenset(example_variants) == frozenset(
            find_equivalent_variant(variant_dict, chunk_provider))


def test_find_equivalent_variant_whole_seq(fetch_seq_mock_data):
    with patch.object(bioutils.seqfetcher, 'fetch_seq', side_effect=lambda ac, s, e: fetch_seq_mock_data[(str(ac), str(s), str(e))]):
        gene_config_path = os.path.join(pwd, 'test_files', 'gene_config_test.txt')

        cfg = load_config(gene_config_path)
        regions = extract_gene_regions_dict(cfg, 'start_hg38_legacy_variants', 'end_hg38_legacy_variants').keys()
        seq_wrapper = seq_utils.SeqRepoWrapper(regions_preload=regions)

        # empty case
        assert [] == find_equivalent_variants_whole_seq({}, seq_wrapper)

        # a bunch of variants. If they appear in the same set, they are considered equivalent
        example_variants = [
            frozenset({'chr13:g.32355030:A>AA'}),
            frozenset({'chr13:g.32339774:GAT>G', 'chr13:g.32339776:TAT>T'}),
            frozenset({'chr17:g.43090921:G>GCA', 'chr17:g.43090921:GCA>GCACA'})
        ]

        # construct variant dict (flattening example_variants!)
        variant_dict = {v: VCFVariant(
            int(v.split(':')[0].lstrip('chr')),
            int(v.split(':')[1].lstrip('g.')),
            v.split(':')[2].split('>')[0],
            v.split(':')[2].split('>')[1]) for eq_variants in example_variants for v
        in
            eq_variants
        }

        whole_seq_provider = seq_utils.WholeSeqSeqProvider(seq_wrapper)

        assert frozenset(example_variants) == frozenset(
            find_equivalent_variants_whole_seq(variant_dict, whole_seq_provider))


def test_chunking():
    def chunker(vars, margin):
        return seq_utils.ChunkBasedSeqProvider.generate_chunks(vars, margin)

    margin = 2
    assert chunker([], margin) == []

    some_chr = 7
    v1 = VCFVariant(some_chr, 10, 'AGAGT', 'G')
    assert chunker([v1], margin) == [(some_chr, 10 - margin, 10 + 5 + margin)]

    # merging intervals of v1 and v2 together
    v2 = VCFVariant(some_chr, 11, 'GAGT', 'G')
    assert chunker([v1, v2], margin) == [(some_chr, min(10 - margin, 11 - margin), max(10 + 5 + margin, 11 + 4 + margin))]


if __name__ == "__main__":
    # To reproduce failure conditions, paste them in here and run as
    # python ./test_variant_merging.py
    # print variant_equal(v1 = ('17', 41100001, 'gcttccca', ''), v2 = ('17', 41100002, 'cttcccag', ''), version = 'hg38')
    print variant_equal(('13', 32800003, '', 'A'), ('13', 32800005, '', 'A'), 'hg19')
    pass
