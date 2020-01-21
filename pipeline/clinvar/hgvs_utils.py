import hgvs
import hgvs.parser
import hgvs.dataproviders.uta
import hgvs.assemblymapper

import logging
import time

# TODO: does not belong here. move to utilties or similar. but then need to change how luigi stuff is ran.
class HGVSWrapper:
    def __init__(self, hgvs_dp=hgvs.dataproviders.uta.connect(pooling=True)):
        self.hgvs_dp = hgvs_dp
        logging.info("Connecting to %s", self.hgvs_dp.url)
        self.hgvs_hp = hgvs.parser.Parser()
        self.hgvs_am = hgvs.assemblymapper.AssemblyMapper(self.hgvs_dp)


    # TODO: think about interface

    def compute_protein_change(self, hgvs_cdna, include_braces=True, retry_cnts = 5):

        try:
            v = self.hgvs_hp.parse_hgvs_variant(hgvs_cdna)
            vp = self.hgvs_am.c_to_p(v)

            if include_braces and vp and vp.posedit:
                vp.posedit.uncertain = True
        except hgvs.exceptions.HGVSDataNotAvailableError as e:
            time.sleep(1)

            if(retry_cnts > 0):
                logging.warn("Retrying another time for %s", hgvs_cdna)

                return self.compute_protein_change(hgvs_cdna, include_braces, retry_cnts - 1)
            else:
                # fail the pipeline. We actually should get data, but don't due to some intermittent issue
                # see also https://github.com/biocommons/hgvs/issues/519
                raise RuntimeError("HGVS data not available for %s . Cannot continue. Error was %s", hgvs_cdna, e.message)
        except Exception as e:
            logging.warning('Issues converting %s. %s type %s', hgvs_cdna, e.message, type(e))
            return None

        return vp
