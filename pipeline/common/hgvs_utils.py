import logging
import os

import bioutils
import hgvs
import hgvs.assemblymapper
import hgvs.dataproviders.uta
import hgvs.normalizer
import hgvs.parser
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

    def parse_hgvs_string(self, hgvs_string):
        hp = self.hgvs_parser
        hgvs_obj = None
        try:
            hgvs_obj = hp.parse_hgvs_variant(hgvs_string)
        except hgvs.exceptions.HGVSParseError as e:
            logging.info("HGVS Parsing error " + str(e))
        return(hgvs_obj)

    def to_cdna(self, hgvs_obj, target_transcript=None):
        cdna_hgvs_obj = None
        if hgvs_obj.type == 'c':
            genomic_hgvs_obj = self.nm_to_genomic(hgvs_obj)
            if genomic_hgvs_obj:
                cdna_hgvs_obj = self.genomic_to_cdna(genomic_hgvs_obj,
                                                     target_transcript=target_transcript)
        elif hgvs_obj.type == 'g':
            cdna_hgvs_obj = self.genomic_to_cdna(genomic_hgvs_obj,
                                                 target_transcript=target_transcript)
        return(cdna_hgvs_obj)
    
        
    def genomic_to_cdna(self, hgvs_obj, assembly=GRCh38_Assem,
                        target_transcript=None):
        am = self.hgvs_ams[assembly]

        try:
            tr = am.relevant_transcripts(hgvs_obj)

            if tr:
                if target_transcript:
                    return am.g_to_c(hgvs_obj, target_transcript)
                else:
                    return am.g_to_c(hgvs_obj, tr[0])

        except HGVSError as e:
            logging.info("CDNA conversion issues " + str(e))

        return None

    def cdna_to_protein(self, hgvs_cdna, return_str=True):
        if not hgvs_cdna:
            return None

        try:
            if return_str:
                return str(self.hgvs_ams[HgvsWrapper.GRCh38_Assem].c_to_p(
                    hgvs_cdna))
            else:
                return self.hgvs_ams[HgvsWrapper.GRCh38_Assem].c_to_p(
                    hgvs_cdna)
        except (hgvs.exceptions.HGVSError, IndexError, ValueError) as e:
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

    __instance = None

    def hg19_to_hg38(self, v):
        """
        Conversion from hg19 (GRCh37) to hg38 via NM_ transcripts (doesn't work for variants outside transcripts)
        :param v:
        :return:
        """
        am37 = self.hgvs_ams[self.GRCh37_Assem]

        transcripts = [t for t in am37.relevant_transcripts(v) if t.startswith('NM_')]

        if not transcripts:
            raise ValueError("Didn't find transcripts for " + str(v))

        v_c = am37.g_to_c(v, transcripts[0])
        return self.hgvs_ams[self.GRCh38_Assem].c_to_g(v_c)

    def u_to_genomic(self, v, target_assembly=GRCh38_Assem):
        v37 = hgvs.assemblymapper.AssemblyMapper(self.hgvs_dp,
                                                 assembly_name=self.GRCh37_Assem,
                                                 alt_aln_method='BLAST').n_to_g(v)

        if target_assembly == self.GRCh38_Assem:
            return self.hg19_to_hg38(v37)
        elif target_assembly == self.GRCh37_Assem:
            return v37
        else:
            raise ValueError("Unknown assembly " + target_assembly)

    def ng_to_genomic(self, v, target_assembly=GRCh38_Assem):
        am = self.hgvs_ams[target_assembly]

        rel = [t for t in am.relevant_transcripts(v) if t.startswith('NM_')]

        if not rel:
            logging.warning("No transcripts could be found for " + str(v) + " in " + str(am.relevant_transcripts(v)) +
                         " and target assembly " + str(target_assembly))
            return None

        v_c = am.g_to_c(v, rel[0])

        return self.hgvs_ams[target_assembly].c_to_g(v_c)

    def nm_to_genomic(self, v, target_assembly=GRCh38_Assem):
        try:
            genomic_v = self.hgvs_ams[target_assembly].c_to_g(v)
            return(genomic_v)
        except hgvs.exceptions.HGVSDataNotAvailableError as e:
            logging.info("Assembly mapping data not available " + str(e))
        except hgvs.exceptions.HGVSInvalidIntervalError as e:
            logging.info("Out of bounds error " + str(e))
        return(None)

    @staticmethod
    def get_instance():
        if not HgvsWrapper.__instance:
            HgvsWrapper.__instance = HgvsWrapper()
        return HgvsWrapper.__instance
