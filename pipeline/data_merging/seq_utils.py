import os
from collections import namedtuple

SeqWithStart = namedtuple("SeqWithStart", "sequence, start")

class SeqProvider:
    '''
    Provides reference sequences.

    Very prototypy, will evolve by calling an appropriate external API.

    '''

    def __init__(self, reference_path):
        self._ref_dict = {17: {"start": 43000000,
                               "sequence": open(os.path.join(reference_path,
                                                             "brca1_hg38.txt"),
                                                "r").read().upper()},
                          13: {"start": 32300000,
                               "sequence": open(os.path.join(
                                   reference_path, "brca2_hg38.txt"),
                                   "r").read().upper()}
                          }

    def get_seq_with_start(self, chr):
        if chr not in {13, 17}:
            raise ValueError(
                "Only chromosomes 13 and 17 are supported but got {}".format(
                    chr))
        return SeqWithStart(self._ref_dict[chr]["sequence"], self._ref_dict[chr]["start"])
