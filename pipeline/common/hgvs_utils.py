import logging
import os

import bioutils
import hgvs
import hgvs.assemblymapper
import hgvs.dataproviders.uta
import hgvs.normalizer
import hgvs.parser
import hgvs.projector
import hgvs.validator
from hgvs.exceptions import HGVSError


class HgvsWrapper:
    GRCh38_Assem = 'GRCh38'
    GRCh37_Assem = 'GRCh37'

    def __init__(self):
        logging.info("HGVS_SEQREPO_DIR: {}".format(os.environ.get('HGVS_SEQREPO_DIR', "Not set. Using public instance")))

        self.hgvs_dp = hgvs.dataproviders.uta.connect()
        logging.info("Using UTA instance at {}".format(self.hgvs_dp.url))

        self.hgvs_parser = hgvs.parser.Parser()
        self.hgvs_norm = hgvs.normalizer.Normalizer(self.hgvs_dp)

        assemblies = [self.GRCh37_Assem, self.GRCh38_Assem]

        self.hgvs_ams = { a : (hgvs.assemblymapper.AssemblyMapper(self.hgvs_dp, assembly_name=a, normalize=False, prevalidation_level=None)) for a in assemblies }

        all_assemblies = bioutils.assemblies.get_assemblies()

        self.contig_maps = {}
        for a in assemblies:
            m = {}
            for s in all_assemblies[a]['sequences']:
                if s['refseq_ac'] is not None:
                    m[s['name']] = s['refseq_ac']
            self.contig_maps[a] = m

    def to_cdna(self, hgvs_obj, assembly=GRCh38_Assem):
        am = self.hgvs_ams[assembly]

        try:
            tr = am.relevant_transcripts(hgvs_obj)

            if tr:
                return am.g_to_c(hgvs_obj, tr[0])

        except HGVSError as e:
            logging.info("CDNA conversion issues " + str(e))

        return None

    def to_protein(self, hgvs_cdna):
        if not hgvs_cdna:
            return None

        try:
            return str(self.hgvs_ams[HgvsWrapper.GRCh38_Assem].c_to_p(
                hgvs_cdna))
        except (hgvs.exceptions.HGVSError, IndexError) as e:
            logging.info(
                "Protein conversion issues with " + str(hgvs_cdna) + ": " + str(
                    e))

        return None

    def normalizing(self, v):
        if v:
            try:
                return self.hgvs_norm.normalize(v)
            except (hgvs.exceptions.HGVSError, IndexError) as e:
                logging.info(
                    "Issues with normalizing " + str(v) + ": " + str(e))
        return None
