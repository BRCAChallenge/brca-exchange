#!/usr/bin/env python
"""
clinVarParse: parse the ClinVar XML file and output the data of interest
"""
import argparse
import logging
import xml.etree.ElementTree as ET
import re

from clinvar import clinvar_common as clinvar
from common import utils


def printHeader():
    print("\t".join(("HGVS", "Submitter", "ClinicalSignificance",
                     "DateLastUpdated", "DateSignificanceLastEvaluated", "SCV",
                     "SCV_Version", "ID", "Origin", "Method", "Genomic_Coordinate",
                     "Symbol", "Protein", "Description", "SummaryEvidence",
                     "ReviewStatus", "ConditionType", "ConditionValue",
                     "ConditionDB_ID", "Synonyms")))

MULTI_VALUE_SEP = ','

def processSubmission(submissionSet, assembly):
    ra = submissionSet.referenceAssertion
    classification = submissionSet.classification
    variant = submissionSet.variant

    if variant is None:
        logging.warning("No variant information could be extracted for ReferenceClinVarAssertion ID %s %s",
                     submissionSet.referenceAssertion.id, [c.accession for c in submissionSet.otherAssertions.values()])
        return None

    debug = False
    for oa in list(submissionSet.otherAssertions.values()):
        if not ("somatic" in oa.origin and len(oa.origin) == 1):
            hgvs = variant.hgvs_cdna
            if assembly in variant.coordinates:
                synonyms = MULTI_VALUE_SEP.join(variant.synonyms)
                vcf_var = variant.coordinates[assembly]

                # Omit the variants that don't have any genomic start coordinate indicated.
                if vcf_var and _bases_only(vcf_var.ref) and _bases_only(vcf_var.alt):
                    print("\t".join((str(hgvs),
                                     oa.submitter,
                                     str(oa.clinicalSignificance),
                                     str(oa.dateLastUpdated),
                                     str(oa.dateSignificanceLastEvaluated),
                                     str(oa.accession),
                                     str(oa.accession_version),
                                     str(oa.id),
                                     ",".join(oa.origin),
                                     ",".join(oa.method),
                                     str(vcf_var).replace('g.', ''), #change
                                     str(variant.geneSymbol),
                                     str(variant.proteinChange),
                                     ",".join(oa.description),
                                     str(oa.summaryEvidence),
                                     str(oa.reviewStatus),
                                     str(classification.condition_type),
                                     str(classification.condition_value),
                                     ",".join(classification.condition_db_id),
                                     str(synonyms))))
                    #print("\t".join((str(hgvs),
                    #                 oa.submitter,
                    #                 str(oa.clinicalSignificance),
                    #                 str(oa.dateLastUpdated),
                    #                 str(oa.dateSignificanceLastEvaluated),
                    #                 str(oa.accession),
                    #                 str(oa.accession_version),
                    #                 str(oa.id),
                    #                 str(oa.origin),
                    #                 str(oa.method),
                    #                 str(vcf_var).replace('g.', ''),
                    #                 str(variant.geneSymbol),
                    #                 str(proteinChange),
                    #                 str(oa.description),
                    #                 str(oa.summaryEvidence),
                    #                 str(oa.reviewStatus),
                    #                 str(ra.condition_type),
                    #                 str(ra.condition_value),
                    #                 ",".join(ra.condition_db_id) if isinstance(ra.condition_db_id, list) else str(ra.condition_db_id),
                    #                 str(synonyms))))


def _bases_only(seq):
    # only allow bases, a not other IUPAC codes such as N, B, S etc
    return all(s in set(['-', 'A', 'C', 'T', 'G']) for s in seq)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("clinVarXmlFilename")
    parser.add_argument('-a', "--assembly", default="GRCh38")
    parser.add_argument('-l', "--logs")
    args = parser.parse_args()

    utils.setup_logfile(args.logs)

    printHeader()


    with open(args.clinVarXmlFilename) as inputFile:
        for event, elem in ET.iterparse(inputFile, events=('start', 'end')):
            if event == 'end' and elem.tag == 'VariationArchive':
                if clinvar.isCurrent(elem):
                    submissionSet = clinvar.variationArchive(elem, debug=False)
                    if submissionSet.valid:
                        processSubmission(submissionSet, args.assembly)
                elem.clear()

if __name__ == "__main__":
    # execute only if run as a script
    main()
