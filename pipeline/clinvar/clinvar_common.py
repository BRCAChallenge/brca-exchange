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
    if ff is not None:
        return ff.text
    elif field in element.attrib:
        return element.attrib[field]
    else:
        return None

def findUniqueElement(name, parent):
    """Find a child element directly or indirectly underneath this parent
       element which should occur only once (i.e. there should be no other
       elements of the same name).  Test this assumption and return
       the child element"""
    child_found = False
    for next_child in parent.iter(name):
        assert(child_found is False)
        child = next_child
        child_found = True
    return(child)



#def build_xpath_filter_for_cv_assertions(gene_symbols):
##    symbols_str = [ f'text()="{s}"' for s in gene_symbols] ### HERE
#    symbols_pred = ' or '.join(symbols_str)

    # filter assertion if it contains a Symbol we are interested in
#    return f"ReferenceClinVarAssertion/MeasureSet/Measure/MeasureRelationship/Symbol/ElementValue[({symbols_pred}) and @Type=\"Preferred\"]"

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


def extract_genomic_coordinates_from_location(loc):
    """
    loc: `xml` module object of a ClinVar `Location` element
    returns: dictionary of assembly (str) to genomic coordinates (VCFVariant object)
    """
    coords = {}
    for sl in loc.iter('SequenceLocation'):
        assembly = sl.attrib['Assembly'] # GRCh38
        if sl.get('referenceAlleleVCF'):
            coords[assembly] = variant_utils.VCFVariant(
                int(sl.attrib['Chr']),
                int(sl.attrib['positionVCF']),
                sl.attrib['referenceAlleleVCF'],
                sl.attrib['alternateAlleleVCF']
            )
    # if no reference/alternate allele found, compute (assuming genomic coordinates
    # are either present for all assemblies or for none)
    #if not coords:
    #    coords = _extract_genomic_coordinates_from_non_genomic_fields(meas_el)
    return coords


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
    """The variant (SimpleAllele) set.  """

    def __init__(self, element, debug=True):
        self.element = element
        self.valid = True
        self.variantType = textIfPresent(element, "VariantType")
        if self.variantType in ["Deletion", "Duplication", "Insertion"]:
            self.valid = False
            return
        if "AlleleID" in element.attrib:
            self.id = element.get("AlleleID")
        else:
            self.id = None
        if debug:
            print("Parsing variant", self.id)
        self.geneSymbol = None
        for gene in element.iter("GeneList/Gene"):
            geneSymbol = gene.get("Symbol")
            if self.geneSymbol is None:
                self.geneSymbol = geneSymbol
            else:
                self.geneSymbol = self.geneSymbol + "," + geneSymbol
        location = element.find("Location")
        if location is None:
            self.valid = False
            return
        self.coordinates = extract_genomic_coordinates_from_location(location)
        self.hgvs_cdna = None
        self.proteinChange = None
        self.synonyms = list()
        spdi = textIfPresent(element, "CanonicalSPDI")
        pc = textIfPresent(element, "ProteinChange")
        if pc is not None:
            self.synonyms.append(pc)
        if element.find("HGVSlist"):
            for hgvs in element.find("HGVSlist").iter("HGVS"):
                pe = hgvs.find("ProteinExpression")
                if pe:
                    proteinHgvs = pe.find("Expression").text
                    if debug:
                        print("protein HGVS", proteinHgvs)
                    self.synonyms.append(proteinHgvs)
                ne = hgvs.find("NucleotideExpression")
                if ne:
                    nucleotideHgvs = ne.find("Expression").text
                    if debug:
                        print("Nucleotide hgvs", nucleotideHgvs)
                    self.synonyms.append(nucleotideHgvs)
                    if "MANESelect" in ne.attrib:
                        if ne.attrib["MANESelect"] == "true":
                            self.hgvs_cdna = nucleotideHgvs
                            if pe:
                                self.proteinChange = proteinHgvs
                                if debug:
                                    print("MANE Select HGVS cDNA ",
                                          nucleotideHgvs,
                                          "protein", proteinHgvs)
                
                


class referenceAssertion:
    """For gathering the reference assertion"""

    def __init__(self, element, debug=False):
        self.element = element
        self.valid = True
        self.title = element.get("Title")
        if debug:
            print("Parsing ReferenceClinVarAssertion", self.title)
        gc = element.find(".//GermlineClassification")
        if gc is None:
            self.valid = False
            return
        self.reviewStatus = textIfPresent(gc, "ReviewStatus")
        self.clinicalSignificance = textIfPresent(gc, "Description")
        if self.clinicalSignificance is not None:
            self.dateSignificanceLastEvaluated = gc.find("Description").get("DateLastEvaluated")
        else:
            self.dateSignificanceLastEvaluated = None
        self.summaryDescription = None


class classification:
    """For gathering data on the trait.  This code expects only one trait"""
    def __init__(self, element, debug=False):
        self.element = element        
        self.condition_type = None
        self.condition_value = None
        self.condition_db_id = list()
        for trait in element.iter("Trait"):
            self.condition_type = trait.get("Type")
            for name in trait.iter("Name"):
                ev = name.find("ElementValue")
                if ev.get("Type") == "Preferred":
                    self.condition_value = ev.text
                for xref in trait.iter("XRef"):
                    xref_string = xref.get("DB") + "_" + xref.get("ID")
                    if not xref_string in self.condition_db_id:
                        self.condition_db_id.append(xref_string)



class clinicalAssertion:
    """Class for representing one submission (i.e. one annotation of a
dir(referen    submitted variant"""

    
    
    def __init__(self, element, debug=False):
        self.element = element
        self.id = element.get("ID")
        if debug:
            print("Parsing ClinicalAssertion", self.id)
        self.dateSubmitted = element.get("SubmissionDate")
        self.dateLastUpdated = element.get("DateLastUpdated")
        cva = element.find("ClinVarAccession")
        if cva == None:
            self.accession = None
        else:
            self.accession = cva.get("Accession")
            self.accession_version = cva.get("Version")
            self.submitter = cva.get("SubmitterName")
        self.reviewStatus = textIfPresent(element, "ReviewStatus")
        classif = element.find("Classification")
        self.dateSignificanceLastEvaluated = classif.get("DateLastEvaluated")
        self.clinicalSignificance = textIfPresent(classif,
                                                  "GermlineClassification")
        self.summaryEvidence = textIfPresent(classif, "Comment")
        oil = element.find("ObservedInList")
        self.origin = list()
        self.method = list()
        self.description = list()
        for oi in oil.iter("ObservedIn"):
            sample = oi.find("Sample")
            if sample != None:
                origin = textIfPresent(sample, "Origin")
                if not origin in self.origin:
                    self.origin.append(origin)
            method = oi.find("Method")
            if method != None:
                newMethod = textIfPresent(method, "MethodType")
                if not newMethod in self.method:
                    self.method.append(newMethod)
            od = oi.find("ObservedData")
            if od != None:
                for attr in od.iter("Attribute"):
                    if attr.get("Type") == 'Description':
                        newDescription = textIfPresent(od, "Attribute")
                    else:
                        newDescription = "None"
                    if not newDescription in self.description:
                        self.description.append(newDescription)



class variationArchive:
    """Container class for a variationArchive record, which is a set of 
    submissions
    that were submitted to ClinVar together.  In the ClinVar terminology,
    each ClinVarSet is one "SimpleAllele" with one aggregate record 
   ("RCV Accession"), which contains
    one or more submissions ("SCV Accessions").
    """

    def __init__(self, element, debug=True):
        self.element = element
        self.id = element.get("VariationID")
        if debug:
            print("Parsing ClinVarSet ID", self.id)
        self.valid = True
        sa = element.find("./ClassifiedRecord/SimpleAllele")
        if sa is None:
            self.valid = False
            return
        self.variant = variant(sa, debug=False)
        if not self.variant.valid:
            self.valid = False
            return

        #
        # Look for the RCVAccession object.  There should be exactly one.
        rcva = element.find("./ClassifiedRecord/RCVList/RCVAccession")
        if rcva is None:
            self.valid = False
            return
        self.referenceAssertion = referenceAssertion(rcva, debug=False)
        if not self.referenceAssertion.valid:
            self.valid = False
            return
        
        #
        # Look for the GermlineClassification object.
        cl = element.find("./ClassifiedRecord/Classifications")
        self.classification = classification(cl, debug=debug)
        
        self.otherAssertions = dict()

        for item in element.iter("ClinicalAssertion"):
            if isCurrent(item):
                ca = clinicalAssertion(item)
                accession = ca.accession
                self.otherAssertions[accession] = ca
        
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
