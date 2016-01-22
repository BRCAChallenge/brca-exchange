#!/usr/bin/env python
"""
clinVarParse: parse the ClinVar XML file and output the data of interest
"""
import argparse
import dipper.utils.ClinVar as clinvar 
import importlib
import xml.etree.ElementTree as ET
from xml.dom import minidom
import sys

# Below is a magic formula that keeps the code from choking on characters
# beyond the standard ASCII set.
#importlib.reload(sys)
#sys.setdefaultencoding('utf8')

def prettify(elem):
    """Return a pretty-printed XML string for the Element.
    """
    rough_string = ET.tostring(elem, 'utf-8')
    reparsed = minidom.parseString(rough_string)
    return reparsed.toprettyxml(indent="  ")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("clinVarXmlFilename")
    args = parser.parse_args()
    tree = ET.parse(args.clinVarXmlFilename)
    root = tree.getroot()
    newRoot = ET.Element(root.tag, attrib=root.attrib)
    for cvs in root.findall("ClinVarSet"):
        if clinvar.isCurrent(cvs):
            submissionSet = clinvar.clinVarSet(cvs)
            variant = submissionSet.referenceAssertion.variant
            if variant != None:
                if variant.geneSymbol == "BRCA1" or variant.geneSymbol == "BRCA2":
                    newRoot.append(cvs)
    print(prettify(newRoot))

if __name__ == "__main__":
    # execute only if run as a script
    main()
