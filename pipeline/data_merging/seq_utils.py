import logging
import os
from collections import namedtuple

from biocommons.seqrepo import SeqRepo
from bioutils import assemblies
from intervaltree import IntervalTree
from toolz import groupby

SeqWithStart = namedtuple("SeqWithStart", "sequence, start")

ChrInterval = namedtuple("ChrInterval", "chr, start, end")


def build_interval_trees_by_chr(chr_intervals, interval_tuple_builder):
    d = {}

    for c, regs in groupby(lambda r: r.chr,
                           chr_intervals).iteritems():
        interval_tuples = [
            (r.start, r.end, interval_tuple_builder(c, r.start, r.end)) for
            r in regs]

        d[c] = IntervalTree.from_tuples(
            interval_tuples)

    return d

class SeqRepoWrapper:
    DEFAULT_ASSY_NAME = 'GRCh38.p11'  # TODO: move somewhere else?

    def __init__(self, seq_repo_path=None, regions_preload=None, preload_pos_margin=100):
        if not seq_repo_path:
            seq_repo_path = os.environ.get("HGVS_SEQREPO_DIR")

        self.seq_repo = SeqRepo(seq_repo_path)
        self.assy_map = assemblies.make_name_ac_map(self.DEFAULT_ASSY_NAME)

        self.preloaded_regions = {}
        if regions_preload:
            self.preloaded_regions = build_interval_trees_by_chr(regions_preload, lambda c, s, e: self._fetch_seq(c, s, e + preload_pos_margin))

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
        ac = self.assy_map[str(chr)]
        return self.seq_repo.fetch(ac, start_pos, end_pos)


class WholeSeqSeqProvider:
    def __init__(self, seq_wrapper):
        if not seq_wrapper.preloaded_regions:
            raise ValueError("need to have preloaded regions in seq wrapper")

        self.seq_wrapper = seq_wrapper

    def get_seq_with_start(self, chr, pos):
        preloaded = self.seq_wrapper.get_preloaded_seq_at(chr, pos)

        if not preloaded:
            raise ValueError("Expected to have a sequence preloaded at chr {} pos {}".format(chr, pos))

        first_pos, _, seq = preloaded[0]
        return SeqWithStart(seq, first_pos)


# TODO rename?
class ChunkBasedSeqProvider:
    def __init__(self, variant_records, margin, seq_wrapper):
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


# TODO mark as deprecated
class LegacyFileBasedSeqProvider:
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
