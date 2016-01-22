#!/usr/bin/env python
"""
clinVarParse: parse the ClinVar XML file and output the data of interest
"""
import argparse
import dipper.utils.ClinVar as clinvar 
import codecs
import sys
import xml.etree.ElementTree as ET

def printHeader():
    print("\t".join(("Chrom", "Pos", "Ref", "Alt", "Symbol", "Name", "ID", 
                     "RefClinicalSignificance", "Submitter", 
                     "ClinicalSignificance")))

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("clinVarXmlFilename")
    parser.add_argument('-a', "--assembly", default="GRCh37")
    args = parser.parse_args()

    printHeader()


    tree = ET.parse(args.clinVarXmlFilename)
    root = tree.getroot()
    for cvs in root.findall("ClinVarSet"):
        if clinvar.isCurrent(cvs):
            submissionSet = clinvar.clinVarSet(cvs)
            ra = submissionSet.referenceAssertion
            for oa in submissionSet.otherAssertions.values():
                submitter = oa.submitter
                if oa.method != "literature only" or oa.submitter == "Counsyl":
                    if oa.origin != "somatic" or oa.clinicalSignificance != "none provided":
                        variant = ra.variant
                        if variant != None:
                            chrom = None
                            start = None
                            referenceAllele = None
                            alternateAllele = None
                            if args.assembly in variant.coordinates:
                                genomicData = variant.coordinates[args.assembly]
                                chrom = genomicData.chrom
                                start = genomicData.start
                                referenceAllele = genomicData.referenceAllele
                                alternateAllele = genomicData.alternateAllele
                            print("\t".join((str(chrom), str(start),
                                             str(referenceAllele),
                                             str(alternateAllele),
                                             str(variant.geneSymbol), 
                                             str(variant.name), 
                                             str(variant.id),
                                             str(ra.clinicalSignificance), 
                                             str(oa.submitter), 
                                             str(oa.clinicalSignificance))))
                        

if __name__ == "__main__":
    # execute only if run as a script
    main()
