"""
Shared fixtures mainly around mocking requests via bioutils.seqfetcher to an
external webservice during tests. This requests to the webservice can be problematic
especially on the CI (intermittent failures).

When new tests are added, GENERATE_MOCK_DATA can be set to True, s.t. the webservice
is called and the result of this call (a sequence) added to a file in the 'data' directory.
After this happened, the GENERATE_MOCK_DATA can be set to False and the call to that
webservice is mocked using the data from the file just generated.
"""
import glob
import os

import bioutils
import pytest
from bioutils import seqfetcher
from mock import patch

from common import seq_utils
from common.hgvs_utils import HgvsWrapper


pwd = os.path.dirname(os.path.realpath(__file__))
data_dir = os.path.join(pwd, 'data')


@pytest.fixture(scope="session")
def fetch_seq_mock_data():
    mock_data = {}

    for path in glob.glob(os.path.join(data_dir, 'mock_fetch_seq*')):
        with open(path, 'r') as f:
            fn = os.path.basename(path).rstrip('.txt')
            key = tuple(fn.split('-')[-3:]) # ac, start, end
            mock_data[key] = f.readline().strip()

    return mock_data

# see comment in file header
GENERATE_MOCK_DATA = False

def generate_mock_data(ac, s, e):
    """
    Calls webservice and stores data into a file for test usage e.g. in CI
    """
    fetched = bioutils.seqfetcher.fetch_seq(ac, s, e)

    with open(os.path.join(data_dir, "mock_fetch_seq-" + str(ac) + '-' + str(s) + '-' + str(e) + ".txt"), 'w') as ff:
        ff.write(fetched)

    return fetched


@pytest.fixture(scope="session")
def hgvs_wrapper(fetch_seq_mock_data):
    if not GENERATE_MOCK_DATA:
        with patch.object(bioutils.seqfetcher, 'fetch_seq',
                              side_effect=lambda ac, s, e: fetch_seq_mock_data[(str(ac), str(s), str(e))]):
            return HgvsWrapper()
    else:
        with patch.object(bioutils.seqfetcher, 'fetch_seq',
                              side_effect=lambda ac, s, e: generate_mock_data(ac, s, e)):
            return HgvsWrapper()


@pytest.fixture(scope="module")
def seq_fetcher(fetch_seq_mock_data):
    if not GENERATE_MOCK_DATA:
        with patch.object(bioutils.seqfetcher, 'fetch_seq',
                              side_effect=lambda ac, s, e: fetch_seq_mock_data[(str(ac), str(s), str(e))]):
            return seq_utils.SeqRepoWrapper()
    else:
        with patch.object(bioutils.seqfetcher, 'fetch_seq',
                              side_effect=lambda ac, s, e: generate_mock_data(ac, s, e)):
            return seq_utils.SeqRepoWrapper()


