"""
ClinVarUtils: basic
"""

import re
from common import hgvs_utils, variant_utils
from hgvs.exceptions import HGVSError
import hgvs
import logging
import requests
import xmltodict

clinvar_eutils = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

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
        return ff.text


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

def build_xpath_filter_for_cv_assertions(gene_symbols):
    symbols_str = [ f'text()="{s}"' for s in gene_symbols]
    symbols_pred = ' or '.join(symbols_str)

    # filter assertion if it contains a Symbol we are interested in
    return f"ReferenceClinVarAssertion/MeasureSet/Measure/MeasureRelationship/Symbol/ElementValue[({symbols_pred}) and @Type=\"Preferred\"]"


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
            if a.text is not None:
                sy.append(a.text)

    return sy + sy_alt


def extractSynonymsFromApiRecord(vcv_data):
    synonyms = []
    assert('VariationArchive' in vcv_data['ClinVarResult-Set'])
    assert('InterpretedRecord' in vcv_data['ClinVarResult-Set']['VariationArchive'])
    simple_allele = vcv_data['ClinVarResult-Set']['VariationArchive']['InterpretedRecord']['SimpleAllele']
    hgvs_list = simple_allele['HGVSlist']
    hgvs_type = {'NucleotideExpression', 'ProteinExpression'}
    for hgvs_record in hgvs_list['HGVS']:
        for this_type in hgvs_type:
            if (this_type in hgvs_record):
                synonyms.append(hgvs_record[this_type]['Expression'])
    if 'ProteinChange' in simple_allele:
        if isinstance(simple_allele['ProteinChange'], str):
            synonyms.append(simple_allele['ProteinChange'])
        elif isinstance(simple_allele['ProteinChange'], list):
            for item in simple_allele['ProteinChange']:
                synonyms.append(item)
    if 'OtherNameList' in simple_allele:
        if isinstance(simple_allele['OtherNameList']['Name'], str):
            synonyms.append(simple_allele['OtherNameList']['Name'])
    return(synonyms)



def extract_genomic_coordinates_from_measure(meas_el):
    """
    meas_el: `xml` module object of a ClinVar `Measure` element

    returns: dictionary of assembly (str) to genomic coordinates (VCFVariant object)
    """

    sequence_locations = meas_el.findall('SequenceLocation')

    coords = {}
    for el in sequence_locations:
        assembly = el.attrib['Assembly'] # GRCh38

        if el.get('referenceAlleleVCF'):
            coords[assembly] = variant_utils.VCFVariant(
                int(el.get('Chr')),
                int(el.get('positionVCF')),
                el.get('referenceAlleleVCF'),
                el.get('alternateAlleleVCF')
            )

    # if no reference/alternate allele found, compute (assuming genomic coordinates
    # are either present for all assemblies or for none)
    if not coords:
        coords = _extract_genomic_coordinates_from_non_genomic_fields(meas_el)

    return coords

def extract_genomic_coordinates_from_location(location):
    coords = {}
    for assembly in location['SequenceLocation']:
        if assembly['@AssemblyStatus'] == 'current':
            if '@referenceAlleleVCF' in assembly:
                assembly_name = assembly["@Assembly"]
                coords[assembly_name] = variant_utils.VCFVariant(int(assembly['@Chr']),
                                                                 int(assembly['@positionVCF']),
                                                                 assembly['@referenceAlleleVCF'],
                                                                 assembly['@alternateAlleleVCF'])

def _preprocess_element_value(var_str):
     # removing dangling protein level changes like '(p.Glu2198fs)'
    return re.sub(r'\s*\(p[^)]+\)', '', var_str)


def _extract_genomic_coordinates_from_non_genomic_fields(meas_el, assemblies = [hgvs_utils.HgvsWrapper.GRCh38_Assem], hgvs_wrapper = hgvs_utils.HgvsWrapper.get_instance()):
    pref_el_lst = meas_el.findall('Name/ElementValue[@Type="Preferred"]')

    coords = {}
    if not pref_el_lst:
        return coords

    pref = pref_el_lst[0].text
    preprocessed_var = _preprocess_element_value(pref)

    hutils = hgvs_wrapper.get_instance()

    try:
        v = hgvs_wrapper.hgvs_parser.parse(preprocessed_var)

        for assembly in assemblies:
            if v.ac.startswith('U'):
                v_g = hutils.u_to_genomic(v, assembly)
            elif v.ac.startswith('NG_'):
                v_g = hutils.ng_to_genomic(v, assembly)
            elif v.ac.startswith('NM_'):
                v_g = hutils.nm_to_genomic(v, assembly)
            else:
                logging.warning("Skipping genomic coordinate extraction for " + preprocessed_var)
                continue

            if v_g:
                vcf = variant_utils.VCFVariant.from_hgvs_obj(v_g)
                coords[assembly] = vcf
    except HGVSError as e:
        logging.warning("HGVS Error while attempting to process " + preprocessed_var + " : " + str(e))

    return coords

def variant_ids_from_gene(gene, max_expected_variants = 100000):
    """
    Given a gene, return a list of variants (by ClinVar ID) for that gene
    """
    api_url = "%s/esearch.fcgi?db=clinvar&term=%s[gene]&retmax=%d&format=json" \
                                   % (clinvar_eutils, gene, max_expected_variants)
    response = requests.get(api_url)
    rj = response.json()
    number_variants = rj["esearchresult"]["count"]
    logging.info("Retrieved %s IDs for gene %s" % (number_variants, gene))
    assert(int(number_variants) == len(rj["esearchresult"]["idlist"]))
    return(rj["esearchresult"]["idlist"])

def variant_summary_from_id(id):
    """
    Given a ClinVar variant ID, return the summary record for the variant
    """
    variant_summary_url = "%s/esummary.fcgi?db=clinvar&id=%s&retmode=json" \
                                                % (clinvar_eutils, id)
    variant_summary_response = requests.get(variant_summary_url)
    variant_summary = variant_summary_response.json()
    return(variant_summary)


def vcv_from_accession(accession):
    """ 
    Given the accession of a ClinVar VCV record, retrieve the record
    The VCV is the aggregate record representing the set of all submissions
    to ClinVar for the indicated ID, with any phenotype.  The RCV record aggregates
    information from all submissions for a specific phenotype.
    """
    vcv_fetch_url = "%s/efetch.fcgi?db=clinvar&rettype=vcv&id=%s" % (clinvar_eutils, accession)
    vcv_fetch_response = requests.get(vcv_fetch_url)
    vcv_data = xmltodict.parse(vcv_fetch_response.text)
    logging.debug("retrieving data for variant ID %s, VCV %s, url %s" % (id, accession, vcv_fetch_url))
    return(vcv_data)
    

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

    def __init__(self, variant_obj, hgvs_protein_refseq, gene):
        self.id = variant_obj['measure_id']
        if debug:
            print("Parsing variant", self.id)
        self.name = variant_obj['variation_name']
        self.attribute = dict()
        self.attribute['HGVS, protein, RefSeq'] = hgvs_protein_refseq
        self.coordinates = extract_genomic_coordinates() # here
        self.geneSymbol = gene
        
        


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
            logging.warning("Assertion with ID " + str(self.id) + " has multiple measures. Taking first one.")
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

    def __init__(self, this_id, variant_summary, vcv_data, gene, debug=False):
        self.origin = None
        self.ethnicity = None
        self.geographicOrigin = None
        self.age = None
        self.gender = None
        self.familyData = None
        self.method = None
        self.variant = None
        self.synonyms = []
        self.condition_type = None
        self.condition_value = None
        self.condition_db_id = None
        self.id = this_id
        if 'variation_set' in variant_summary:
            self.variant = variant(variant_summary['result'][self.id]['variation_set'][0],
                                   self._find_hgvs_protein(vcv_data), gene)
        if 'clinical_significance' in variant_summary:
            clinical_significance = variant_summary['result'][self.id]['clinical_significance']
            self.reviewStatus = clincial_significance['review_status']
            self.clinicalSignificance = clinical_significance['description']
            self.dateSignificanceLastEvaluated = clinical_significance['last_evaluated']
        else:
            self.reviewStatus = None
            self.clinicalSignificance = None
            self.dateSignificanceLastEvaluated = None
        traitSet = variant_summary['result'][self.id]['trait_set']
        self.condition_value = traitSet[0]['trait_name']
        self.condition_db_id = []
        for xref in traitSet[0]['trait_xrefs']:
            self.condition_db_id.append(xref['db_source'] + "_" + xref['db_id'])
        self.synonyms = extractSynonymsFromApiRecord(vcv_data)

        def _find_hgvs_protein(vcv_data):
            hgvs_protein = None
            simple_allele = vcv_data['ClinVarResult-Set']['VariationArchive']['InterpretedRecord']['SimpleAllele']              
            hgvs_list = simple_allele['HGVSlist']['HGVS']
            for hgvs_record in hgvs_list:
                if 'NucleotideExpression' in hgvs_record:
                    if 'ProteinExpression' in hgvs_record:
                        if '@MANESelect' in hgvs_record['NucleotideExpression']:
                            if hgvs_record['NucleotideExpression']['@MANESelect'] == "true":
                                hgvs_protein = hgvs_record['ProteinExpression']['Expression']
            return(hgvs_protein)

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

        
    def __init__(self, clinical_assertion, debug=False):
        self.id = clinical_assertion['@ID']
        if debug:
            print("Parsing ClinVarAssertion", self.id)
        if not 'ClinVarSubmissionID' in clinical_assertion:
            self.submitter = None
            self.dateSubmitted = None
        else:
            self.dateSubmitted = clinical_assertion['@SubmissionDate']
        if not 'ClinVarAccession' in clinical_assertion:
            self.accession = None
        else:
            cva = clinical_assertion['ClinVarAccession']
            self.submitter = cva['@SubmitterName']
            self.accession = cva['@Accession']
            self.accession_version = cva['@Version']
        self.origin = None
        self.method = None
        self.description = None
        if 'ObservedInList' in clinical_assertion:
            oi = clinical_assertion['ObservedInList']['ObservedIn']
            if 'Sample' in oi:
                if 'Origin' in oi['Sample']:
                    self.origin = oi['Sample']['Origin']
            if 'Method' in oi:
                if 'MethodType' in oi['Method']:
                    self.method = oi['Method']['MethodType']
            if 'ObservedData' in oi:
                if 'Attribute' in oi['ObservedData']:
                    attr = oi['ObservedData']['Attribute']
                    if attr['@Type'] == "Description":
                        self.description = attr['#text']
        if 'ReviewStatus' in clinical_assertion:
            self.reviewStatus = clinical_assertion['ReviewStatus']
        else:
            self.reviewStatus = None
        if 'Interpretation' in clinical_assertion:
            self.clinicalSignificance = clinical_assertion['Interpretation']['Description']
            self.dateSignificanceLastUpdated = clinical_assertion['Interpretation']['@DateLastEvaluated']
            self.summaryEvidence = None  # Here
        else:
            self.clinicalSignificance = None
            self.summaryEvidence = None
            self.dateSignificanceLastUpdated = None
        self.synonyms = None


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

    def __init__(self, this_id, gene, debug=False):
        """Given the ID of a ClinVar variant, build and represent the set of submissions"""
        self.this_id = this_id
        variant_summary = variant_summary_from_id(this_id)
        assert('result' in variant_summary)
        assert(this_id in variant_summary['result'])
        vcv_data = vcv_from_accession(variant_summary['result'][this_id]['accession'])
        self.referenceAssertion = referenceAssertion(this_id, variant_summary, vcv_data, gene, debug=debug)
        self.otherAssertions = dict()
        clinical_assertion_list = vcv_data['ClinVarResult-Set']['VariationArchive']['InterpretedRecord']['ClinicalAssertionList']
        if len(clinical_assertion_list) == 1:
            cva = clinVarAssertion(clinical_assertion_list['ClinicalAssertion'])
            accession = cva.accession
            self.otherAssertions[accession] = cva
        self.referenceAssertion.hgvs_cdna = extract_hgvs_cdna_from_api_record(variant_summary['result'][this_id]['title'],
                                                                              None) ### Here
        print("hgvs cdna", self.referenceAssertion.hgvs_cdna)
        
        
    def extract_hgvs_cdna(self, variant_name, clinvar_set_el):
        """
        Finds a HGVS CDNA representation of a variant within a ClinVarSet.
        If possible, avoid repeat representations using the "[]" synatax, since
        we are currently not able to handle it further downstream (https://github.com/biocommons/hgvs/issues/113)

        :param variant_name: variant name from title
        :param clinvar_set_el: clinvar set element
        :return: HGVS CDNA representation as string
        """
        hgvs_cand = re.sub(r"\(" + "(BRCA[1|2])" + r"\)",
                      "", variant_name.split()[0])

        # only take variants starting with NM_ and not containing []
        hgvs_cdna_re = r'NM_.*:[^\[]*$'

        if not re.match(hgvs_cdna_re, hgvs_cand):
            # taking Attribute of 'HGVS', 'HGVS, coding' or 'HGVS, coding, RefSeq' in
            # both ReferenceClinVarAssertion and ClinVarAssertion's
            hgvs_candidates = [ev.text for ev in clinvar_set_el.findall('.//Measure/AttributeSet/Attribute')
                               if ev.attrib['Type'].startswith('HGVS, coding') or ev.attrib['Type'] == 'HGVS']

            filtered = [s for s in hgvs_candidates if re.match(hgvs_cdna_re, s)]
            if filtered:
                return filtered[0]

        return hgvs_cand


def extract_hgvs_cdna_from_api_record(variant_name, clinvar_set_el):
    """
    Finds a HGVS CDNA representation of a variant within a ClinVarSet.
    If possible, avoid repeat representations using the "[]" synatax, since
    we are currently not able to handle it further downstream (https://github.com/biocommons/hgvs/issues/113)

    :param variant_name: variant name from title
    :param clinvar_set_el: clinvar set element
    :return: HGVS CDNA representation as string
    """
    hgvs_cand = re.sub(r"\(" + "[A-Za-z0-9]*" + r"\)",
                       "", variant_name.split()[0])

    # only take variants starting with NM_ and not containing []
    hgvs_cdna_re = r'NM_.*:[^\[]*$'
    if re.match(hgvs_cdna_re, hgvs_cand):
        return(hgvs_cand)
    else:
        return(None)

    #if not re.match(hgvs_cdna_re, hgvs_cand):
    #    # taking Attribute of 'HGVS', 'HGVS, coding' or 'HGVS, coding, RefSeq' in
    #    # both ReferenceClinVarAssertion and ClinVarAssertion's
    #    hgvs_candidates = [ev.text for ev in clinvar_set_el.findall('.//Measure/AttributeSet/Attribute')
    #                       if ev.attrib['Type'].startswith('HGVS, coding') or ev.attrib['Type'] == 'HGVS']
#
#            filtered = [s for s in hgvs_candidates if re.match(hgvs_cdna_re, s)]
#            if filtered:

