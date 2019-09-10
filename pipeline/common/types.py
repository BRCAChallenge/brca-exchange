from collections import namedtuple
import hgvs


class VCFVariant(namedtuple("VCFVariant", "chr,pos,ref,alt")):
    __slots__ = ()

    def to_hgvs_obj(self, contig_ac_map):
        ref = self.ref
        alt = self.alt

        # TODO check updated library
        start = self.pos
        end = start + len(
            ref) - 1  # TODO: more than SNP, TODO: verify how to calculated

        ac = contig_ac_map[str(self.chr)]

        var_g = hgvs.sequencevariant.SequenceVariant(ac=ac,
                                                     type='g',
                                                     posedit=hgvs.posedit.PosEdit(
                                                         hgvs.location.Interval(
                                                             start=hgvs.location.SimplePosition(
                                                                 start),
                                                             end=hgvs.location.SimplePosition(
                                                                 end),
                                                             uncertain=False),
                                                         hgvs.edit.NARefAlt(
                                                             ref=ref if ref != '' else None,
                                                             alt=alt if alt != '' else None,
                                                             uncertain=False)))

        return var_g

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
    def from_hgvs_obj(hgvs_var):
        chr = int(hgvs_var.ac.split("_")[1].split('.')[0])

        alt = hgvs_var.posedit.edit.alt if hasattr(hgvs_var.posedit.edit, 'alt') else 'noalt'

        return VCFVariant(chr, hgvs_var.posedit.pos.start, hgvs_var.posedit.edit.ref, alt)
