import hgvs
import hgvs.parser
import hgvs.dataproviders.uta
import hgvs.assemblymapper

import logging


# TODO: does not belong here. move to utilties or similar. but then need to change how luigi stuff is ran.
class HGVSWrapper:
    def __init__(self):
        self.hgvs_dp = hgvs.dataproviders.uta.connect()
        logging.info("Connecting to %s", self.hgvs_dp.url)
        self.hgvs_hp = hgvs.parser.Parser()
        self.hgvs_am = hgvs.assemblymapper.AssemblyMapper(self.hgvs_dp)


    # TODO: think about interface
    def compute_protein_change(self, hgvs_cnda, include_braces=True):
        try:
            v = self.hgvs_hp.parse_hgvs_variant(hgvs_cnda)
            vp = self.hgvs_am.c_to_p(v)

            if include_braces and vp and vp.posedit:
                vp.posedit.uncertain = True
        except Exception as e:
            logging.warning('Issues converting %s. %s', hgvs_cnda, e.message)
            return None
        return vp