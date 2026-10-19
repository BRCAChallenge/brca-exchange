#!/usr/bin/env python
"""
clinVarParse: parse the ClinVar XML file and output the data of interest
"""
import argparse
import logging
import xml.etree.ElementTree as ET
import re

from clinvar import clinvar_common as clinvar
from common import config, utils


# Columns describing the VCV-level aggregate classification.  These are the
# same on every row (submission) of a given variant.
VCV_COLUMNS = ("VCV_VariationID", "VCV_Accession", "VCV_Version",
               "VCV_ClinicalSignificance", "VCV_ReviewStatus",
               "VCV_DateLastEvaluated", "VCV_NumberOfSubmissions",
               "VCV_NumberOfSubmitters", "VCV_MostRecentSubmission",
               "VCV_Explanation", "VCV_Conditions", "VCV_ConditionDB_IDs",
               "VCV_DateLastUpdated")


def printHeader():
    print("\t".join(("HGVS", "Submitter", "ClinicalSignificance",
                     "DateLastUpdated", "DateSignificanceLastEvaluated", "SCV",
                     "SCV_Version", "ID", "Origin", "Method", "Genomic_Coordinate",
                     "Symbol", "Protein", "Description", "SummaryEvidence",
                     "ReviewStatus", "ConditionType", "ConditionValue",
                     "ConditionDB_ID", "Synonyms", "BIC_Nomenclature")
                    + VCV_COLUMNS))

MULTI_VALUE_SEP = ','
CONDITION_SEP = '|'


def _vcv_value(value):
    """Render an aggregate classification value, '-' if absent.  ';' becomes
    ',' because convert_tsv_to_vcf.py rewrites ';' (its INFO delimiter) to '.',
    which would mangle an Explanation like 'Pathogenic(3); Benign(1)'."""
    if value is None or value == '':
        return '-'
    return str(value).replace(';', ',')


def aggregateFields(agg):
    """The values for VCV_COLUMNS, in order"""
    conditions = CONDITION_SEP.join(agg.conditions)
    # Keep empty per-condition slots so IDs stay aligned with VCV_Conditions
    condition_db_ids = (CONDITION_SEP.join(MULTI_VALUE_SEP.join(ids)
                                           for ids in agg.conditionDbIds)
                        if any(agg.conditionDbIds) else None)
    return tuple(_vcv_value(v) for v in (
        agg.variationID, agg.accession, agg.version,
        agg.clinicalSignificance, agg.reviewStatus, agg.dateLastEvaluated,
        agg.numberOfSubmissions, agg.numberOfSubmitters,
        agg.mostRecentSubmission, agg.explanation,
        conditions, condition_db_ids,
        agg.dateLastUpdated))


def processSubmission(submissionSet, assembly):
    classification = submissionSet.classification
    variant = submissionSet.variant

    if variant is None:
        logging.warning("No variant information could be extracted for VariationArchive %s %s",
                        submissionSet.id, [c.accession for c in submissionSet.otherAssertions.values()])
        return None

    hgvs = submissionSet.name
    aggregate = aggregateFields(submissionSet.aggregateClassification)
    for oa in list(submissionSet.otherAssertions.values()):
        if ("somatic" in oa.origin and len(oa.origin) == 1):
            logging.warning("HGVS %s because submissions are only somatic in origin", hgvs)
        else:
            if not assembly in variant.coordinates:
                logging.warning("HGVS %s rejected for poor variant coordinate data", hgvs)
            else:
                synonyms = MULTI_VALUE_SEP.join(variant.synonyms)
                vcf_var = variant.coordinates[assembly]

                # Omit the variants that don't have any genomic start coordinate indicated.
                if not (vcf_var and _bases_only(vcf_var.ref) and _bases_only(vcf_var.alt)):
                    logging.warning("HGVS %s rejected for poor VCF data", hgvs)
                else:
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
                                     str(synonyms),
                                     variant.bic_nomenclature or '-')
                                    + aggregate))


def _bases_only(seq):
    # only allow bases, a not other IUPAC codes such as N, B, S etc
    return all(s in set(['-', 'A', 'C', 'T', 'G']) for s in seq)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("clinVarXmlFilename")
    parser.add_argument('-a', "--assembly", default="GRCh38")
    parser.add_argument('-l', "--logs")
    parser.add_argument('--configfile',
                        help="path to gene configuration file; when given, "
                             "GeneList entries whose gene's known chromosome "
                             "doesn't match the variant's own coordinates are rejected")
    args = parser.parse_args()

    utils.setup_logfile(args.logs)

    gene_chromosomes = None
    if args.configfile:
        gene_config_df = config.load_config(args.configfile)
        gene_chromosomes = dict(zip(gene_config_df[config.SYMBOL_COL], gene_config_df[config.CHROM_COL]))

    printHeader()


    with open(args.clinVarXmlFilename) as inputFile:
        for event, elem in ET.iterparse(inputFile, events=('start', 'end')):
            if event == 'end' and elem.tag == 'VariationArchive':
                if clinvar.isCurrent(elem):
                    submissionSet = clinvar.variationArchive(elem, gene_chromosomes=gene_chromosomes, debug=False)
                    if submissionSet.valid:
                        processSubmission(submissionSet, args.assembly)
                    else:
                        if hasattr(submissionSet, 'name'):
                            logging.warning("Submission %s not valid", submissionSet.name)
                        else:
                            logging.warning("Submission not valid, name not available")
                elem.clear()

if __name__ == "__main__":
    # execute only if run as a script
    main()
