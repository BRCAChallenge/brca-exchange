import pytest

from insilico import spliceai_runner


def sample_formatter():
    """A stand-in for get_delta_scores.

    The format string is written out as a literal rather than referenced
    through the module, because reformat() rewrites constants in the code
    object and a global lookup would not be one.
    """
    return '{}|{}|{:.2f}|{:.2f}|{:.2f}|{:.2f}|{}|{}|{}|{}'.format(
        'A', 'BRCA1', 0.0174, 0.0021, 0.0, 0.0067, 946, -3, 881, 886)


def test_two_decimal_places_floor_small_scores():
    # The behaviour we are correcting: 0.0174 reaches the VCF as 0.02, and
    # anything under 0.005 as a flat 0.00.
    assert sample_formatter() == 'A|BRCA1|0.02|0.00|0.00|0.01|946|-3|881|886'


def test_reformat_raises_precision_without_touching_positions():
    patched = spliceai_runner.reformat(sample_formatter, 3)
    assert patched() == 'A|BRCA1|0.017|0.002|0.000|0.007|946|-3|881|886'


def test_reformat_honours_other_precisions():
    assert spliceai_runner.reformat(sample_formatter, 4)().split('|')[2] == '0.0174'
    assert spliceai_runner.reformat(sample_formatter, 2)() == sample_formatter()


def test_reformat_leaves_the_original_alone():
    spliceai_runner.reformat(sample_formatter, 3)
    assert sample_formatter().split('|')[2] == '0.02'


def test_reformat_refuses_a_function_without_the_known_format():
    def unrecognised():
        return '{}|{}|{:.5f}'.format('A', 'BRCA1', 0.0174)

    with pytest.raises(RuntimeError, match='unknown precision'):
        spliceai_runner.reformat(unrecognised, 3)


def test_take_precision_defaults_and_strips_the_flag():
    argv = ['-I', 'in.vcf', '-O', 'out.vcf', '-D', '4999']
    assert spliceai_runner.take_precision(list(argv)) == (argv, 3)
    assert spliceai_runner.take_precision(argv + ['--precision', '4']) == (argv, 4)
    assert spliceai_runner.take_precision(
        ['--precision', '5'] + argv) == (argv, 5)


def test_take_precision_rejects_nonsense():
    with pytest.raises(SystemExit):
        spliceai_runner.take_precision(['--precision'])
    with pytest.raises(SystemExit):
        spliceai_runner.take_precision(['--precision', '9'])


def test_installed_spliceai_still_uses_the_expected_format():
    """Canary: a SpliceAI upgrade that changes this must not pass silently."""
    utils = pytest.importorskip('spliceai.utils')
    assert spliceai_runner.TWO_DP_FORMAT in utils.get_delta_scores.__code__.co_consts
