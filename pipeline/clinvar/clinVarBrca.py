#!/usr/bin/env python
"""clinVarBrca: parse the ClinVar XML file, and output a subset containing the BRCA records

The ClinVar XML file contains ClinVarSet records and XML header/footer
information.  The purpose of this script is to output the
header/footer information plus those ClinVarSet records containing
BRCA variants.  It reads through the file line by line.  When it's
reading a ClinVarSet record, it concatenates the contents of the
record into a buffer.  When it comes to the end of a ClinVarSet
record, it checks whether the record contains a BRCA variant, and if
so, it prints it.  If it's not reading a ClinVarSet record, then it
echoes each line to stdout.
"""
import argparse
import clinvar
import gzip
import xml.etree.ElementTree as ET
import sys
import logging

# Below is a magic formula that keeps the code from choking on characters
# beyond the standard ASCII set.
reload(sys)
sys.setdefaultencoding('utf8')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("clinVarXmlFilename")
    parser.add_argument('-a', "--artifacts_dir", help='Artifacts directory with pipeline artifact files.')
    args = parser.parse_args()

    logging_level = logging.DEBUG
    log_file_path = args.artifacts_dir + "clinvarbrcapy.log"
    logging.basicConfig(filename=log_file_path, filemode="w", level=logging_level)

    inputBuffer = ""

    with gzip.open(args.clinVarXmlFilename) as inputFile:
        inClinVarSet = False
        for line in inputFile:
            if "<ClinVarSet" in line:
                inHeader = False
                inputBuffer = line
                inClinVarSet = True
            elif "</ClinVarSet>" in line:
                inputBuffer += line
                inClinVarSet = False
                cvs = ET.fromstring(inputBuffer)
                # Ignore non-brca variants
                if clinvar.isCurrent(cvs) and "brca" in inputBuffer.lower():
                    try:
                        submissionSet = clinvar.clinVarSet(cvs)
                    except AttributeError:
                        logging.debug("AttributeError running clinvar.clinVarSet(cvs), inputBuffer: %s, cvs: %s", inputBuffer, cvs)
                        continue
                    variant = submissionSet.referenceAssertion.variant
                    if variant != None:
                        if variant.geneSymbol == "BRCA1" or variant.geneSymbol == "BRCA2":
                            print inputBuffer
                inputBuffer = None
            elif inClinVarSet:
                inputBuffer += line
            else:
                if len(line) > 1:
                    print line.rstrip()

if __name__ == "__main__":
    # execute only if run as a script
    main()
