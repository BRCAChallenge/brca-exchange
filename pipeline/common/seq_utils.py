import logging
import os
from collections import namedtuple

from biocommons.seqrepo import SeqRepo
from bioutils import assemblies, seqfetcher

from .utils import build_interval_trees_by_chr, ChrInterval

SeqWithStart = namedtuple("SeqWithStart", "sequence, start")


class SeqRepoWrapper:
    '''
    Wrap access to biocommons seqrepo.

    Has a mechanism to preload certain genomic regions. Queries falling into these
    regions are then served from memory.

    '''
    DEFAULT_ASSY_NAME = 'GRCh38.p11'

    def __init__(self, seq_repo_path=None, regions_preload=None, preload_pos_margin=500):
        '''
        :param seq_repo_path: Path to local seqrepo directory. If None, read HGVS_SEQREPO_DIR environment variable
        :param regions_preload: Iterable[ChrInterval], optionally preload these genomic regions
        :param preload_pos_margin: adding margin at the end of a preloaded genome
          in order to have data to verify structural variants across the end of a gene
        '''

        if not seq_repo_path:
            seq_repo_path = os.environ.get("HGVS_SEQREPO_DIR")

        if seq_repo_path:
            seq_repo = SeqRepo(seq_repo_path)
            self.seq_repo_fetcher = seq_repo.fetch
        else:
            logging.warn("Using remote sequence provider.")
            self.seq_repo_fetcher = seqfetcher.fetch_seq

        self.assy_map = assemblies.make_name_ac_map(self.DEFAULT_ASSY_NAME)

        self.preloaded_regions = {}
        if regions_preload:
            self.preloaded_regions = build_interval_trees_by_chr(regions_preload,
                                                                 lambda c, s, e: self._fetch_seq(c, s, e + preload_pos_margin))

    def get_seq_at(self, chr, pos, length):
        return self.get_seq(chr, pos, pos + length)

    def get_seq(self, chr, start_pos, end_pos):
        preloaded = self.get_preloaded_seq_at(chr, start_pos)

        if preloaded:
            first_pos, _, seq = preloaded[0]
            return seq[start_pos - first_pos: end_pos - first_pos]
        else:
            return self._fetch_seq(chr, start_pos, end_pos)

    def get_preloaded_seq_at(self, chr, pos):
        preloaded = []
        if chr in self.preloaded_regions:
            preloaded = list(self.preloaded_regions[chr].at(pos))

        return preloaded

    def _fetch_seq(self, chr, start_pos, end_pos):
        accession = self.assy_map[str(chr)]

        # SeqRepo used 0-based coordinates
        return self.seq_repo_fetcher(accession, start_pos - 1, end_pos - 1)

    __instance = None

    @staticmethod
    def get_instance():
        """ Static access method. """
        if not SeqRepoWrapper.__instance:
            SeqRepoWrapper.__instance = SeqRepoWrapper()
        return SeqRepoWrapper.__instance

class WholeSeqSeqProvider:
    '''
    Sequence provider returning sequence of an entire gene to apply edits to verify
    if variants are equivalent
    '''
    def __init__(self, seq_wrapper):
        '''

        :param seq_wrapper: SeqRepoWrapper instance with preloaded regions
        '''
        if not seq_wrapper.preloaded_regions:
            raise ValueError("need to have preloaded regions in seq wrapper")

        self.seq_wrapper = seq_wrapper

    def get_seq_with_start(self, chr, pos):
        preloaded = self.seq_wrapper.get_preloaded_seq_at(chr, pos)

        if not preloaded:
            raise ValueError("Expected to have a sequence preloaded at chr {} pos {}".format(chr, pos))

        first_pos, _, seq = preloaded[0]
        return SeqWithStart(seq, first_pos)


class ChunkBasedSeqProvider:
    '''
    Sequence provider not returning sequence of an entire gene, but only of the 'chunk'
    a variant belongs in.

    Needs to have all variants which are going to be processed up front.
    For every variant the 'interval of influence' is determined.
    This is a range of positions on the genome where differences to the reference may possibly occur.
    It is basically the position on the chromosome plus the length of the ref
    or alt string whatever is longer. Then, to deal with repeat regions and
    differing alignment some constant margin (say 100 nucleotides) is added at each side.

    Then the 'interval of influence' are together into chunks, s.t. all overlapping
    intervals end up in the same chunk. One chunk would then have a start and end
    position associated. After this preprocessing, we determine for every variant
    the chunk it belongs and calculate the edited sequence wrt to the sequence
    corresponding to the chunk.
    '''
    def __init__(self, variant_records, margin, seq_wrapper):
        '''

        :param variant_records: Iterable[VCFVariant] complete list of variants to be processed
        :param margin: margin to add to the left and right of the initial interval of influence of
          a variant to account for repeat regions in the genome
        :param seq_wrapper: instance of SeqRepoWrapper
        '''
        self.sr_wrapper = seq_wrapper

        chunks = self.generate_chunks(variant_records, margin)

        self.itree_dict = build_interval_trees_by_chr([ChrInterval._make(c) for c in chunks], lambda c, s, e: self.sr_wrapper.get_seq(c, s, e))

        logging.debug("Number of chunks: %d", len(chunks))
        logging.debug("Total bytes in memory from chunks: %d", sum(
            [len(inter[2]) for tree in self.itree_dict.values() for inter in
             tree]))


    def get_seq_with_start(self, chr, pos):
        itree = self.itree_dict[chr]

        intervals = list(itree.at(pos))
        assert len(
            intervals) == 1, "expect exactly one chunk per chr, got {}".format(
            len(intervals))
        start, _, seq = intervals[0]

        return (seq, start)

    @staticmethod
    def generate_chunks(variant_records, margin):
        if not variant_records:
            return []

        l = sorted([(v.chr,
                     v.pos - margin,
                     v.pos + max(len(v.ref), len(v.alt)) + margin)
                    for v in variant_records])

        chunks = []  # (start, end) of chunk

        # initial interval
        cur_chr, cur_start, cur_end = l[0]

        for c, start, end in l[1:]:
            if cur_chr != c or start > cur_end:
                chunks.append((cur_chr, cur_start, cur_end))

                # start new chunk
                cur_chr = c
                cur_start = start
                cur_end = end

            cur_end = max(end, cur_end)

        chunks.append((cur_chr, cur_start, cur_end))

        return chunks


class LegacyFileBasedSeqProvider:
    '''
    Sequence provider returning sequence of an entire gene read from a file
    Deprecated, only for maintaining compatibility with testing code.
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

    def get_seq_with_start(self, chr, pos=None):
        if chr not in {13, 17}:
            raise ValueError(
                "Only chromosomes 13 and 17 are supported but got {}".format(
                    chr))
        return SeqWithStart(self._ref_dict[chr]["sequence"],
                            self._ref_dict[chr]["start"])
