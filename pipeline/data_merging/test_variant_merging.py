from hypothesis import given, assume, settings
from hypothesis.strategies import integers, tuples, text, sampled_from, lists
from variant_merging import variant_equal, init
import itertools
import os
import pytest

#
# initialize module
#

pwd=os.path.dirname(os.path.realpath(__file__))

# XXX instead of adding '/' we should fix variant_merging to use
# os.join.
class Args:
    reference = os.path.join(pwd, '..', 'data') + '/'

init(Args())
# This is really hacky, but we have to re-import from variant_merging to
# get the BRCA1 and BRCA2 after init. I'm sure there's a better way.
from variant_merging import BRCA1, BRCA2

# If this isn't true, things are bad.
assert(set(BRCA1.keys()) == set(BRCA2.keys()))
lengths = [len(ref["sequence"]) for ref in BRCA1.values()] + [len(ref["sequence"]) for ref in BRCA1.values()]
assert(len(lengths) > 0)
assert(len(set(lengths)) == 1)

reference_length = lengths[0]

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
reference_id = sampled_from(BRCA1.keys())
chrom = sampled_from(('13', '17'))
pos = integers(min_value=0, max_value=reference_length) # pos relative to the reference start.
base = sampled_from(('A', 'C', 'T', 'G'))
subseq = lists(base).map(lambda x: ''.join(x))

# A variant with all random fields.
# (chrom, position, ref, alt, src)
variant = tuples(chrom, pos, subseq, subseq, text())

# A variant where 'ref' is a length. We will later copy in the reference base
# pairs from a reference sequence.
# (chrom, position, ref, alt, src)
variant_on_ref = tuples(chrom, pos, integers(min_value=0), subseq, text())

chrom_ref = {
    '13': BRCA2,
    '17': BRCA1
}

# Add the start position of the slice of reference that we're holding, to create
# a valid chrom position from an integer in [0, len(reference)]. Also converts
# to 1-based coords.
def add_start(v, ref_id):
    "Add reference start position to a variant"
    (chrom, pos, ref, alt, src) = v
    return (chrom, pos + chrom_ref[chrom][ref_id]['start'] + 1, ref, alt, src)


#
# tests
#

import sys
def spy(msg, x):
    print >> sys.stderr, msg, x
    return x

def test_variant_equal_throws_below_reference():
    ref_id = chrom_ref['13'].keys()[0]
    v1 = ('13', chrom_ref['13'][ref_id]['start'], 'A', 'C', 'foo')
    v2 = ('13', chrom_ref['13'][ref_id]['start'] + 10, 'A', 'C', 'foo')
    with pytest.raises(AssertionError):
        variant_equal(v1, v2, ref_id)
    with pytest.raises(AssertionError):
        variant_equal(v2, v1, ref_id)

def test_variant_equal_throws_above_reference():
    ref_id = chrom_ref['13'].keys()[0]
    start = chrom_ref['13'][ref_id]['start']
    v1 = ('13', start + reference_length + 2, 'A', 'C', 'foo')
    v2 = ('13', start + 10, 'A', 'C', 'foo')
    with pytest.raises(AssertionError):
        variant_equal(v1, v2, ref_id)
    with pytest.raises(AssertionError):
        variant_equal(v2, v1, ref_id)

#@settings(max_examples=50000)
@given(variant, reference_id)
def test_variant_equal_identity(v, ref_id):
    "A variant should be equal to itself"
    v = add_start(v, ref_id)
    assert variant_equal(v, v, ref_id)

#@settings(max_examples=50000)
@given(variant, variant, reference_id)
def test_variant_equal_commutative(v1, v2, ref_id):
    "Comparing x, y should be the same as comparing y, x"
    v1 = add_start(v1, ref_id)
    v2 = add_start(v2, ref_id)
    assert variant_equal(v1, v2, ref_id) == variant_equal(v2, v1, ref_id)


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
    (chrom, pos, ref, alt, src) = variant
    while len(ref) > 0 and len(alt) > 0 and ref[0] == alt[0]:
        pos += 1
        ref = ref[1:]
        alt = alt[1:]
    while len(ref) > 0 and len(alt) > 0 and ref[-1] == alt[-1]:
        ref = ref[:-1]
        alt = alt[:-1]
    return (chrom, pos, ref, alt, src)

def equiv_variants_delete(reference, variant):
    (chrom, pos, ref, alt, src) = normalize_variant(variant)
    assert(alt == '') # must be a simple delete
    l = len(ref)
    equiv_pos = equiv_deletes(reference, pos, l)
    return [(chrom, epos, reference[epos : epos + l], alt, src) for epos in equiv_pos]


# To test that equivalent variants test equal, generate a random variant, then
# compute variants equivalent to it, and assert that they all test equal to the
# first variant. Currently on works for deletes.

#runtimes = 5000
#@settings(max_examples=runtimes, max_iterations=runtimes, timeout=-1, database_file=None)
@given(variant_on_ref, reference_id)
def test_variant_equal_equiv(v, ref_id):
    (chrom, pos, reflen, alt, src) = v
    refsequence = chrom_ref[chrom][ref_id]["sequence"]
    assume(pos + reflen <= len(refsequence))
    v1 = (chrom, pos, refsequence[pos : pos + reflen], alt, src)
    if alt == '' and reflen != 0:
        equivs = equiv_variants_delete(refsequence, v1)
        for veq in equivs:
            assert variant_equal(add_start(v1, ref_id), add_start(veq, ref_id), ref_id)
    else:
        pass

if __name__ == "__main__":
    # To reproduce failure conditions, paste them in here and run as
    # python ./test_variant_merging.py
    #print variant_equal(v1 = ('17', 41100001, 'gcttccca', '', ''), v2 = ('17', 41100002, 'cttcccag', '', ''), version = 'hg38')
    pass
