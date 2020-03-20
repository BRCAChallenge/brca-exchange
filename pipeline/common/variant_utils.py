from collections import namedtuple

import hgvs

from common import seq_utils
from bioutils.sequences import reverse_complement

class VCFVariant(namedtuple("VCFVariant", "chr,pos,ref,alt")):
    __slots__ = ()

    def to_hgvs_obj(self, contig_ac_map):
        ref = self.ref
        alt = self.alt

        start = int(self.pos)
        end = start + len(ref) - 1

        ac = contig_ac_map[str(self.chr)]

        posedit = hgvs.posedit.PosEdit(
            hgvs.location.Interval(
                start=hgvs.location.SimplePosition(start),
                end=hgvs.location.SimplePosition(end),
                uncertain=False),
            hgvs.edit.NARefAlt(
                ref=ref if ref != '' else None,
                alt=alt if alt != '' else None,
                uncertain=False))

        return hgvs.sequencevariant.SequenceVariant(ac=ac,
                                                    type='g',
                                                    posedit=posedit)

    def __str__(self):
        return "chr{}:g.{}:{}>{}".format(self.chr, self.pos, self.ref, self.alt)

    @staticmethod
    def from_str(s):
        return VCFVariant(
            int(s.split(':')[0].lstrip('chr')),
            int(s.split(':')[1].lstrip('g.')),
            s.split(':')[2].split('>')[0],
            s.split(':')[2].split('>')[1])

    @staticmethod
    def from_hgvs_obj(hgvs_var, seq_fetcher=seq_utils.SeqRepoWrapper.get_instance()):
        chr = int(hgvs_var.ac.split("_")[1].split('.')[0])

        alt = hgvs_var.posedit.edit.alt if hasattr(hgvs_var.posedit.edit, 'alt') else ''

        if not alt:
            alt = ''

        edit_type = str(hgvs_var.posedit.edit)

        pos = hgvs_var.posedit.pos.start.base
        ref = hgvs_var.posedit.edit.ref
        if not ref:
            if edit_type.startswith('ins'):
                ref = str(seq_fetcher.get_seq(str(chr), hgvs_var.posedit.pos.start.base, hgvs_var.posedit.pos.start.base + 1))
            else:
                ref = str(seq_fetcher.get_seq(str(chr), hgvs_var.posedit.pos.start.base, hgvs_var.posedit.pos.end.base + 1))

        if len(ref) >= 1 and len(alt) >= 1 and not edit_type.startswith('ins'):
            return VCFVariant(int(chr), int(pos), ref, alt)

        # require padding, i.e. inserting previous base to avoid empty alt
        # e.g. instead of 'C'>'' do 'AC'>'A'
        if edit_type.startswith('del') or edit_type.startswith('ins') or edit_type.startswith('dup') or edit_type.startswith('inv'):
            if not edit_type.startswith('ins') and not edit_type.startswith('inv'):
                pos -= 1

            # transforming 'del' to a delins
            padding = str(seq_fetcher.get_seq_at(str(chr), pos, 1))

            if edit_type.startswith('ins'):
                alt = padding + alt
            elif edit_type.startswith('dup'):
                alt = padding + ref
                ref = padding
            elif edit_type.startswith('del'):
                ref = padding + ref
                alt = padding + alt
            elif edit_type.startswith('inv'):
                alt = reverse_complement(ref)

        return VCFVariant(int(chr), int(pos), ref, alt)
