#!/usr/bin/env python
"""spliceai_runner: run SpliceAI, writing the delta scores at a precision we
choose rather than SpliceAI's hard-coded two decimal places.

SpliceAI 1.3.1 formats its four delta scores with ``{:.2f}`` on the way into
the VCF INFO field, so every score below 0.005 is recorded as ``0.00``.  That
floor is coarse next to the thresholds the ENIGMA BRCA1/2 VCEP applies, and it
makes our scores impossible to reconcile with services that keep the raw
floats -- a variant the ARIANE service scores at 0.017 reaches our database as
a flat zero.

The scoring itself is fine, so this does not fork it.  ``get_delta_scores`` is
rebuilt with its one format-string constant replaced and SpliceAI's own
``main()`` does the rest: same models, same annotation handling, same CLI.  If
a future SpliceAI release changes that constant the patch raises instead of
silently reverting to two decimals.

Usage is SpliceAI's own, plus --precision:

    python spliceai_runner.py -I in.vcf -O out.vcf -R hg38.fa \\
        -A insilico/spliceai_annotations/brca_mane_grch38.txt -D 4999
"""

import sys
import types

DEFAULT_PRECISION = 3

# The single constant inside spliceai.utils.get_delta_scores that fixes the
# precision of DS_AG|DS_AL|DS_DG|DS_DL.  The delta *positions* around it are
# integers and are left alone.
TWO_DP_FORMAT = '{}|{}|{:.2f}|{:.2f}|{:.2f}|{:.2f}|{}|{}|{}|{}'


def reformat(function, precision):
    """Return a copy of `function` whose delta scores carry `precision` digits."""
    wanted = TWO_DP_FORMAT.replace('{:.2f}', '{:.%df}' % precision)
    constants = function.__code__.co_consts
    if TWO_DP_FORMAT not in constants:
        raise RuntimeError(
            'spliceai.utils.get_delta_scores no longer builds its output with '
            '{!r}. Refusing to run rather than write scores at an unknown '
            'precision -- check what the installed SpliceAI now emits and '
            'update TWO_DP_FORMAT.'.format(TWO_DP_FORMAT))
    patched = tuple(wanted if constant == TWO_DP_FORMAT else constant
                    for constant in constants)
    return types.FunctionType(function.__code__.replace(co_consts=patched),
                              function.__globals__, function.__name__,
                              function.__defaults__, function.__closure__)


def take_precision(argv):
    """Pull --precision out of `argv`, leaving SpliceAI's own arguments."""
    if '--precision' not in argv:
        return argv, DEFAULT_PRECISION
    index = argv.index('--precision')
    if index + 1 >= len(argv):
        raise SystemExit('--precision needs a number of decimal places')
    precision = int(argv[index + 1])
    if not 1 <= precision <= 6:
        raise SystemExit('--precision must be between 1 and 6')
    return argv[:index] + argv[index + 2:], precision


def main(argv=None):
    argv, precision = take_precision(
        list(sys.argv[1:] if argv is None else argv))

    import spliceai.utils
    import spliceai.__main__ as spliceai_main

    patched = reformat(spliceai.utils.get_delta_scores, precision)
    spliceai.utils.get_delta_scores = patched
    # __main__ did `from spliceai.utils import get_delta_scores`, so it holds
    # its own reference to the original and has to be patched separately.
    spliceai_main.get_delta_scores = patched

    sys.argv = ['spliceai'] + argv
    spliceai_main.main()


if __name__ == '__main__':
    main()
