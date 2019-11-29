from collections import namedtuple

import hgvs

from common import seq_utils


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

        alt = hgvs_var.posedit.edit.alt if hasattr(hgvs_var.posedit.edit,
                                                   'alt') else ''

        if not alt:
            alt = ''

        edit_type = str(hgvs_var.posedit.edit)

        pos = hgvs_var.posedit.pos.start.base
        ref = hgvs_var.posedit.edit.ref
        if not ref:
            ref = str(seq_fetcher.get_seq(str(chr), hgvs_var.posedit.pos.start.base, hgvs_var.posedit.pos.end.base + 1))

        if len(ref) >= 1 and len(alt) >= 1 and not edit_type.startswith('ins'):
            return VCFVariant(chr, pos, ref, alt)

        # require padding

        if edit_type.startswith('del') or edit_type.startswith('ins') or edit_type.startswith('dup'):
            # transforming 'del' to a delins
            padding = str(seq_fetcher.get_seq_at(str(chr), pos-1, 1))

            if not edit_type.startswith('ins'):
                pos -= 1

            if edit_type.startswith('ins'):
                ref = padding
                alt = padding + alt
            elif not edit_type.startswith('dup'):
                ref = padding + ref
                alt = padding + alt
            else:
                alt = padding + ref
                ref = padding
        return VCFVariant(chr, pos, ref, alt)
