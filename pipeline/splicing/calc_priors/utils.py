# coding=utf-8
import time
import csv


class Benchmark(object):
    def __init__(self,name):
        self.name = name

    def __enter__(self):
        self.start = time.time()

    def __exit__(self,ty,val,tb):
        end = time.time()
        print("%s : %0.3f seconds" % (self.name, end-self.start))
        return False


def isclose(a, b, rel_tol=1e-09, abs_tol=0.0):
    """
    Return True if the values a and b are close to each other and False otherwise. Whether or not two values are
    considered close is determined according to given absolute and relative tolerances.

    (definition from https://docs.python.org/3/library/math.html#math.isclose)

    rel_tol is the relative tolerance – it is the maximum allowed difference between a and b, relative to the larger
    absolute value of a or b. For example, to set a tolerance of 5%, pass rel_tol=0.05. The default tolerance is 1e-09,
    which assures that the two values are the same within about 9 decimal digits. rel_tol must be greater than zero.

    abs_tol is the minimum absolute tolerance – useful for comparisons near zero. abs_tol must be at least zero.

    If no errors occur, the result will be: abs(a-b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol).

    The IEEE 754 special values of NaN, inf, and -inf will be handled according to IEEE rules. Specifically, NaN is not
    considered close to any other value, including NaN. inf and -inf are only considered close to themselves.
    :param a:
    :param b:
    :param rel_tol: the relative tolerance
    :param abs_tol: the minimum absolute tolerance
    :return: True if the values a and b are close to each other and False otherwise.
    """
    return abs(a - b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)


class MismatchException(Exception):
    pass


def handle_mismatch(x, warn):
    if warn:
        print x
    else:
        raise MismatchException(x)


def approximate_compare_tsv(aPath, bPath, warnOnMismatch=False):
    """
    For two TSV files with the same columns and rows, ensures that the values of each cell
    betweeen the two files are at least approximately equal if the cell contains a float-parseable value,
    and exactly equal if not.
    :param aPath: path to TSV file 'a'
    :param bPath: path to TSV file 'b'
    :return: True if they are approximately identical, False otherwise
    """
    with open(aPath, "r") as originalFP, open(bPath, "r") as outputFP:
        orig_variants = csv.DictReader(originalFP, delimiter="\t")
        new_variants = csv.DictReader(outputFP, delimiter="\t")

        for line_no, (old, new) in enumerate(zip(orig_variants, new_variants)):
            for col in old:
                try:
                    if not isclose(float(new[col]), float(old[col])):
                        handle_mismatch(
                            "Numeric mismatch on line %d (%s, %s), cell %s: %s vs. %s" % (
                                line_no+1,
                                old["Gene_Symbol"], old["HGVS_cDNA"],
                                col, old[col], new[col]
                            ),
                            warn=warnOnMismatch
                        )
                except ValueError:
                    if new[col] != old[col]:
                        handle_mismatch(
                            "Non-numeric mismatch on line %d (%s, %s), cell %s: '%s' vs. '%s'" % (
                                line_no+1,
                                old["Gene_Symbol"], old["HGVS_cDNA"],
                                col, old[col], new[col]
                            ),
                            warn=warnOnMismatch
                        )

        return True
