#!/usr/bin/env python
"""
Wrapper class for interpreting GA4GH VR objects.
"""
from ga4gh.core import ga4gh_digest
from ga4gh.vr import __version__, ga4gh_identify, ga4gh_serialize, models, normalize
from ga4gh.vr.extras.dataproxy import SeqRepoRESTDataProxy
from ga4gh.vr.extras.translator import Translator
import hgvs.assemblymapper
import hgvs.dataproviders.uta
import hgvs.parser
import json
import sys

class Mapper:
    """
    Utility class for mapping variants to hg38
    """

    parser = None
    hg38_mapper = None

    def __init__(self):
        self.parser = hgvs.parser.Parser()
        hdp = hgvs.dataproviders.uta.connect()
        self.hg38_mapper = hgvs.assemblymapper.AssemblyMapper(hdp,
                                                    assembly_name='GRCh38',
                                                    alt_aln_method='splign',
                                                    replace_reference=True)

    def _parser(self, variant_string):
        """
        Parses a variant from string representation to HGVS representation
        """
        try:
            parsed_variant = self.parser.parse_hgvs_variant(variant_string)
        except hgvs.exceptions.HGVSParseError:
            sys.stderr.write("Parse error: invalid variant %s\n" % (variant_string))
            parsed_variant = None
        except hgvs.exceptions.HGVSInvalidVariantError:
            sys.stderr.write("Error: invalid variant %s\n"  % (variant_string));
            parsed_variant = None
        return(parsed_variant)

    def map_variant(self, variant_string):
        """
        Given a variant represented as a HGVS string, return the variant 
        mapped to the target assembly in hgvs representation
        """
        mapped_variant = None
        parsed_variant = self._parser(variant_string)
        if parsed_variant is not None:
            try:
                mapped_variant = self.hg38_mapper.c_to_g(parsed_variant)
            except hgvs.exceptions.HGVSInvalidIntervalError:
                sys.stderr.write("Error: invalid interval (out of bounds?) %s\n"
                                 % (variant_string))
                parsed_variant = None
            except hgvs.exceptions.HGVSInvalidVariantError:
                sys.stderr.write("Error: invalid variant %s\n" % (variant_string))
                parsed_variant = None
        return(mapped_variant)
                        


class VRUtils:
    """
    Utility class for converting to VR variant representation
    """

    hg38_mapper = None
    translator = None
    
    def __init__(self, seqrepo_url, hg38_mapper=None):
        if hg38_mapper is None:
            self.hg38_mapper = Mapper()
        else:
            self.hg38_mapper = hg38_mapper
        dp = SeqRepoRESTDataProxy(base_url=seqrepo_url)
        self.translator = Translator(data_proxy=dp)

    def to_vr(self, variant_string, verbose=False):
        hg38_variant = self.hg38_mapper.map_variant(variant_string)
        if hg38_variant is None:
            vr_variant = None
        else:
            if verbose:
                print("Variant", variant_string, "translated to",
                      str(hg38_variant))
            try:
                vr_variant = self.translator.from_hgvs(str(hg38_variant))
            except ValueError:
                vr_variant = None
        return(vr_variant)
                                         


class VR:
    """Class for interpreting VR variant representation"""
    
    def __init__(self, ga4gh_vr):
        self.vr = ga4gh_vr.as_dict()

    def start(self):
        """Return the start coordinate"""
        return self.vr["location"]["interval"]["start"]

    def end(self):
        """Return the start coordinate"""
        return self.vr["location"]["interval"]["end"]

    def alt(self):
        """Return the alt base(s)"""
        return self.vr['state']['sequence']
