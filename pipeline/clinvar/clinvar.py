"""
ClinVarUtils: basic
"""

import re
from common import hgvs_utils, variant_utils
from hgvs.exceptions import HGVSError
import hgvs
import logging

def isCurrent(element):
    """Determine if the indicated clinvar set is current"""
    rr = element.find("RecordStatus")
    if rr == None:
        return False
    else:
        return(rr.text == "current")

def textIfPresent(element, field):
    """Return the text associated with a field under the element, or
    None if the field is not present"""
    ff = element.find(field)
    if ff == None or ff.text == None:
        return None
    else:
        return(ff.text.encode('utf-8'))


def processClinicalSignificanceElement(el, obj):
    if el != None:
        obj.reviewStatus = textIfPresent(el, "ReviewStatus")
        obj.clinicalSignificance = textIfPresent(el, "Description")
        obj.summaryEvidence = textIfPresent(el, "Comment")
        obj.dateSignificanceLastEvaluated = el.get('DateLastEvaluated', None)
    else:
        obj.reviewStatus = None
        obj.clinicalSignificance = None
        obj.summaryEvidence = None
        obj.dateSignificanceLastEvaluated = None


def extractSynonyms(el):
    include_types = {'ProteinChange3LetterCode', 'ProteinChange1LetterCode',
                     'nucleotide change', 'protein change, historical'}
    include_types_norm = {s.lower() for s in include_types}

    exclude_hgvs = {'HGVS, protein, RefSeq', 'HGVS, coding, RefSeq'}
    exclude_hgvs_norm = {s.lower() for s in exclude_hgvs}

    sy_alt = [a.text for a in el.findall(
        'MeasureSet/Measure/Name/ElementValue') if a.get('Type').lower() ==
              'alternate']

    sy = []
    for a in el.findall('MeasureSet/Measure/AttributeSet/Attribute'):
        type = a.get('Type').lower()
        if type in include_types_norm or ('hgvs' in type and type
                                               not in exclude_hgvs_norm):
            sy.append(a.text)

    return sy + sy_alt


# TODO tune interface (which element)
# TODO: add assembly
def extract_genomic_coordinates_from_measure(meas_el):
    sequence_locations = meas_el.findall('SequenceLocation')

    coords = {}
    for el in sequence_locations:
        assembly = el.attrib['Assembly'] # GRCh38

        if el.get('referenceAlleleVCF'):
            coords[assembly] = variant_utils.VCFVariant(
                el.get('Chr'),
                el.get('positionVCF'),
                el.get('referenceAlleleVCF'),
                el.get('alternateAlleleVCF')
            )

    # if no reference/alternate allele found, compute (assuming the missingness
    # is consistent across the different assemblies)
    if not coords:
        coords = _extract_genomic_coordinates_from_non_genomic_fields(meas_el)

    return coords


def _preprocess_variant(var_str):
    var_str = re.sub(r'\s*\(p[^)]+\)', '', var_str)
    #var_str = re.sub(r'(del[TCGA]+)ins[0-9]+$', r'\1', var_str)
    # test re.sub('ins[0-9]+$', '', 'NM_007294.3(BRCA1):c.1387_1390delAAAAins4')
    return var_str


def _extract_genomic_coordinates_from_non_genomic_fields(meas_el, assemblies = [hgvs_utils.HgvsWrapper.GRCh38_Assem], hgvs_wrapper = hgvs_utils.HgvsWrapper.get_instance()):
    pref_el_lst = meas_el.findall('Name/ElementValue[@Type="Preferred"]')

    coords = {}
    if not pref_el_lst:
        return coords

    pref = pref_el_lst[0].text

    try:
        preprocessed_var = _preprocess_variant(pref)
        v = hgvs_wrapper.hgvs_parser.parse(preprocessed_var)

        for assembly in assemblies:
            if v.ac.startswith('U'):
                v37 = hgvs.assemblymapper.AssemblyMapper(hgvs_wrapper.hgvs_dp,
                                                         assembly_name=hgvs_utils.HgvsWrapper.GRCh37_Assem,
                                                         alt_aln_method='BLAST').n_to_g(v)

                if assembly == hgvs_utils.HgvsWrapper.GRCh38_Assem:
                    v_g = hgvs_wrapper.hg19_to_hg38(v37)
                else:
                    v_g = v37

            elif v.ac.startswith('NG_'):
                am38 = hgvs_wrapper.hgvs_ams[hgvs_utils.HgvsWrapper.GRCh38_Assem]

                rel = [t for t in am38.relevant_transcripts(v) if t.startswith('NM_')]

                if not rel:
                    logging.warn("No transcripts could be found for " + preprocessed_var + " in " + str(am38.relevant_transcripts(v)))
                    continue

                v_c = am38.g_to_c(v, rel[0])

                v_g = hgvs_wrapper.hgvs_ams[
                    hgvs_utils.HgvsWrapper.GRCh38_Assem].c_to_g(v_c)
            elif v.ac.startswith('NM_'):
                v_g = hgvs_wrapper.hgvs_ams[assembly].c_to_g(v)
            else:
                logging.warn("Skipping genomic coordinate extraction for " + preprocessed_var)
                continue

            vcf = variant_utils.VCFVariant.from_hgvs_obj(v_g)
            coords[assembly] = vcf
    except HGVSError as e:
        logging.warn("HGVS Error while attempting to process " + preprocessed_var + " : " + str(e))

    return coords


class genomicCoordinates:
    """Contains the genomic information on the variant"""

    def __init__(self, element, useNone=False, debug=False):
        if debug:
            print("Parsing genomic coordinates")
        if useNone:
            self.element = None
            self.chrom = None
            self.start = None
            self.stop = None
            self.length = None
            self.referenceAllele = None
            self.alternateAllele = None
        else:
            self.element = element
            self.chrom = element.get("Chr")
            self.stop = element.get("stop")
            self.length = element.get("variantLength")
            self.start = element.get("positionVCF")
            self.referenceAllele = element.get("referenceAlleleVCF")
            self.alternateAllele = element.get("alternateAlleleVCF")


class variant:
    """The Measure set.  We are interested in the variants specifically,
    but measure sets can be other things as well, such as haplotypes"""

    def __init__(self, element, name, id, debug=False):
        self.element = element
        self.id = id
        if debug:
            print("Parsing variant", self.id)
        self.name = name

        self.attribute = dict()
        for attrs in element.findall("AttributeSet"):
            for attrib in attrs.findall("Attribute"):
                self.attribute[attrib.get("Type")] = attrib.text

        self.coordinates = extract_genomic_coordinates_from_measure(element)

        self.geneSymbol = None
        symbols = element.findall("MeasureRelationship/Symbol")
        for symbol in symbols:
            symbol_val = textIfPresent(symbol, "ElementValue")
            if symbol_val.startswith('BRCA'):
                self.geneSymbol = symbol_val



class referenceAssertion:
    """For gathering the reference assertion"""

    def __init__(self, element, debug=False):
        self.element = element
        self.id = element.get("ID")
        if debug:
            print("Parsing ReferenceClinVarAssertion", self.id)

        processClinicalSignificanceElement(element.find(
            "ClinicalSignificance"), self)


        obs = element.find("ObservedIn")
        if obs == None:
            self.origin = None
            self.ethnicity = None
            self.geographicOrigin = None
            self.age = None
            self.gender = None
            self.familyData = None
            self.method = None
        else:
            sample = obs.find("Sample")
            if sample != None:
                self.origin = textIfPresent(sample, "Origin")
                self.ethnicity = textIfPresent(sample, "Ethnicity")
                self.geographicOrigin = textIfPresent(sample, "GeographicOrigin")
                self.age = textIfPresent(sample, "Age")
                self.gender = textIfPresent(sample, "Gender")
                self.familyData = textIfPresent(sample, "FamilyData")
            method = obs.find("Method")
            if method != None:
                self.method = textIfPresent(method, "MethodType")
        self.variant = None
        self.synonyms = []

        measureSet = element.find("MeasureSet")
        #if measureSet.get("Type") == "Variant":

        if len(measureSet.findall("Measure")) > 1:
            logging.warn("Assertion with ID " + str(self.id) + " has multiple measures. Taking first one.")
        if len(measureSet.findall("Measure")) >= 1:
            name = measureSet.find("Name")
            if name == None:
                variantName = None
            else:
                variantName = name.find("ElementValue").text
            self.variant = variant(measureSet.find("Measure"), variantName,
                                   measureSet.get("ID"),
                                   debug=debug)

        self.synonyms = extractSynonyms(element)

        # extract condition
        self.condition_type = None
        self.condition_value = None
        self.condition_db_id = None
        traitSet = element.find("TraitSet")
        if traitSet != None:
            self.condition_type = traitSet.attrib["Type"]
            trait = traitSet.find("Trait")
            if trait != None:
                names = trait.findall("Name")
                if names != None and len(names) > 0:
                    for name in names:
                        ev = name.find("ElementValue")
                        if ev != None and ev.attrib["Type"] == "Preferred":
                            self.condition_value = textIfPresent(name, "ElementValue")
                        break
                xrefs = trait.findall("XRef")
                if xrefs != None and len(xrefs) > 0:
                    self.condition_db_id = []
                    for xref in xrefs:
                        self.condition_db_id.append(xref.attrib["DB"] + "_" + xref.attrib["ID"])



class clinVarAssertion:
    """Class for representing one submission (i.e. one annotation of a
    submitted variant"""

    def __init__(self, element, debug=False):
        self.element = element
        self.id = element.get("ID")
        if debug:
            print("Parsing ClinVarAssertion", self.id)
        cvsd = element.find("ClinVarSubmissionID")
        if cvsd == None:
            self.submitter = None
            self.dateSubmitted = None
        else:
            self.submitter = cvsd.get("submitter", default=None)
            self.dateSubmitted = cvsd.get("submitterDate")
        cva = element.find("ClinVarAccession")
        if cva == None:
            self.accession = None
        else:
            self.accession = cva.get("Acc", default=None)
            self.accession_version = cva.get("Version", default=None)

        self.origin = None
        self.method = None
        self.description = None
        oi = element.find("ObservedIn")
        if oi != None:
            sample = oi.find("Sample")
            if sample != None:
                self.origin = textIfPresent(sample, "Origin")
            method = oi.find("Method")
            if method != None:
                self.method = textIfPresent(method, "MethodType")
            description = oi.find("ObservedData")
            if description != None:
                for attr in description.findall("Attribute"):
                    if attr.attrib["Type"] == 'Description':
                        self.description = textIfPresent(description, "Attribute")

        processClinicalSignificanceElement(element.find(
            "ClinicalSignificance"), self)

        self.dateLastUpdated = cva.get("DateUpdated")

        self.synonyms = extractSynonyms(element)


class clinVarSet:
    """Container class for a ClinVarSet record, which is a set of submissions
    that were submitted to ClinVar together.  In the ClinVar terminology,
    each ClinVarSet is one aggregate record ("RCV Accession"), which contains
    one or more submissions ("SCV Accessions").
    """

    def __init__(self, element, debug=False):
        self.element = element
        self.id = element.get("ID")
        if debug:
            print("Parsing ClinVarSet ID", self.id)
        rcva = element.find("ReferenceClinVarAssertion")
        if isCurrent(rcva):
            self.referenceAssertion = referenceAssertion(rcva, debug=debug)
        self.otherAssertions = dict()

        for item in element.findall("ClinVarAssertion"):
            if isCurrent(item):
                cva = clinVarAssertion(item)
                accession = cva.accession
                self.otherAssertions[accession] = cva

        if self.referenceAssertion.variant:
            self.referenceAssertion.hgvs_cdna = self.extract_hgvs_cdna(self.referenceAssertion.variant.name, element)

    def extract_hgvs_cdna(self, variant_name, clinvar_set_el):
        """
        Finds a HGVS CDNA representation of a variant within a ClinVarSet.
        If possible, avoid repeat representations using the "[]" synatax, since
        we are currently not able to handle it further downstream (https://github.com/biocommons/hgvs/issues/113)

        :param variant_name: variant name from title
        :param clinvar_set_el: clinvar set element
        :return: HGVS CDNA representation as string
        """
        hgvs_cand = re.sub("\(" + "(BRCA[1|2])" + "\)",
                      "", variant_name.split()[0])

        # only take variants starting with NM_ and not containing []
        hgvs_cdna_re = 'NM_.*:[^\[]*$'

        if not re.match(hgvs_cdna_re, hgvs_cand):
            # taking Attribute of 'HGVS', 'HGVS, coding' or 'HGVS, coding, RefSeq' in
            # both ReferenceClinVarAssertion and ClinVarAssertion's
            hgvs_candidates = [ev.text for ev in clinvar_set_el.findall('.//Measure/AttributeSet/Attribute')
                               if ev.attrib['Type'].startswith('HGVS, coding') or ev.attrib['Type'] == 'HGVS']

            filtered = [s for s in hgvs_candidates if re.match(hgvs_cdna_re, s)]
            if filtered:
                return filtered[0]

        return hgvs_cand
