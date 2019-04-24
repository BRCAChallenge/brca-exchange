import os
from collections import namedtuple
from toolz import groupby, memoize
from intervaltree import IntervalTree
from bioutils import assemblies
from biocommons.seqrepo import SeqRepo

SeqWithStart = namedtuple("SeqWithStart", "sequence, start")

ChrInterval = namedtuple("ChrInterval", "chr, start, end")


class SeqRepoWrapper:
    DEFAULT_ASSY_NAME = 'GRCh38.p11'  # TODO: move somewhere else?

    def __init__(self, seq_repo_path=None, regions_preload=None):
        if not seq_repo_path:
            seq_repo_path = os.environ.get("HGVS_SEQREPO_DIR")

        self.seq_repo = SeqRepo(seq_repo_path)
        self.assy_map = assemblies.make_name_ac_map(self.DEFAULT_ASSY_NAME)

        interval_tuples = []
        if regions_preload:
            interval_tuples = [(r.start, r.end, self._fetch_seq(r.chr, r.start, r.end)) for r in regions_preload]

        self.preloaded_regions = IntervalTree.from_tuples(interval_tuples)


    def get_seq_at(self, chr, pos, length):
        return self.get_seq(chr, pos, pos + length)

    def get_seq(self, chr, start_pos, end_pos):
        preloaded = list(self.preloaded_regions.at(start_pos + 1))
        if preloaded:
            first_pos, _, seq = preloaded[0]
            return seq[start_pos - 1 - first_pos : end_pos - 1 - first_pos]
        else:
            return self._fetch_seq(chr, start_pos, end_pos)

    def _fetch_seq(self, chr, start_pos, end_pos):
        ac = self.assy_map[str(chr)]
        return self.seq_repo.fetch(ac, start_pos, end_pos)


# TODO rename?
class ChunkBasedSeqProvider:
    def __init__(self, variant_records, margin, seq_wrapper):
        self.sr_wrapper = seq_wrapper

        chunks = self.generate_chunks(variant_records, margin)

        self.itree_dict = {}

        for c, tuples in groupby(lambda x: x[0], chunks).iteritems():
            interval_tuples = [self._get_interval_tree_el(c, t[1], t[2]) for t in tuples]
            self.itree_dict[c] = IntervalTree.from_tuples(interval_tuples)

        # TODO: check how many chars kept in memory?
        print("number of chunks", len(chunks))
        print(sum([len(inter[2]) for tree in self.itree_dict.values() for inter in tree]))


    def _get_interval_tree_el(self, chr, start, end):
        seq = self.sr_wrapper.get_seq(chr, start - 1, end)
        return start, end, seq

    def get_seq_with_start(self, chr, pos):
        itree = self.itree_dict[chr]

        intervals = list(itree.at(pos))

        assert len(intervals) == 1, "expect only one chunk per chr"
        chunk = intervals[
            0]  # TODO make more robust, should only have one element!

        ret = (chunk[2], chunk[0] - 1)

        # print(ret, pos, chunk[0])
        return ret

    def _repeat_seq(self, chr, start, delta):
        # go over in chunks of 10, find repeat seq
        pass

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
