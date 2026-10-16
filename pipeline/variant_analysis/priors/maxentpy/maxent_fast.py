"""
Pure-Python stand-in for the maxent_fast C-extension module.

The original maxent_fast used a compiled _hashseq extension that only works
with Python 2.7. This wrapper delegates to maxent.py (pure Python) and
accepts the same load_matrix/score5/score3 interface.
"""

from maxentpy import maxent as _maxent


def load_matrix(donor_or_acceptor):
    """Return a preloaded matrix object (here, just a sentinel — ignored by score functions)."""
    return donor_or_acceptor  # passed back as the matrix kwarg; ignored below


def score5(sequence, matrix=None):
    return _maxent.score5(sequence)


def score3(sequence, matrix=None):
    return _maxent.score3(sequence)
