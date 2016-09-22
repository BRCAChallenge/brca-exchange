from hypothesis import given, settings
from hypothesis.strategies import integers, tuples, text, sampled_from, lists
from variant_merging import variant_equal, init
import os

#
# generators
#
chrom = sampled_from(('13', '17'))
pos = integers(min_value=1)
base = sampled_from(('A', 'C', 'T', 'G'))
subseq = lists(base).map(lambda x: ''.join(x))

# (chrom, position, ref, alt, src)
variant = tuples(chrom, pos, subseq, subseq, text())

#
# initialize module
#

pwd=os.path.dirname(os.path.realpath(__file__))

# XXX instead of adding '/' we should fix variant_merging to use
# os.join.
class Args:
    reference = os.path.join(pwd, '..', 'data') + '/'

init(Args())

#
# tests
#

#@settings(max_examples=50000)
@given(variant)
def test_variant_equal_identity(v):
    assert variant_equal(v, v)

#@settings(max_examples=50000)
@given(variant, variant)
def test_variant_equal_commutative(v1, v2):
	assert variant_equal(v1, v2) == variant_equal(v2, v1)

#
#
#

if __name__ == "__main__":
#    variant_equal(('13', 1, 'A', 'T', None), ('13', 2, 'A', 'T', None))
    pass
