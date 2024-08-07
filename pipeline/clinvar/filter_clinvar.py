#!/usr/bin/env python
"""

filter_clinvar.py: Filter the ClinVar VCV Release XML file, and build a subset
XML file containing the VariationArchive records for the genes of interest.

This file filters the XML by reading the XML file as text, assembling text 
objects of the VariationArchive records, and then parsing those records as XML
to see if they specify the target gene symbols.  It's done this way because:
- reading the whole file into one big ElementTree XML object doesn't work 
  because it's too large and we run out of memory
- reading the file iteratively, with ElementTree iterparse hasn't worked, 
  because the ClinVar archive contains some badly-formatted XML records.
  By not parsing the XML as XML until we know that we're looking at an object
  that we want, we can wager that the broken XML records are for other genes.
"""
import argparse
import gzip
import io
import xml.etree.ElementTree as ET

def contains_target_gene_symbol(variationArchiveElement, target_symbol_values):
    root = ET.fromstring(variationArchiveElement)
    geneList = root.find(".//GeneList")
    if geneList:
        for geneElement in geneList.iter("Gene"):
            if "Symbol" in geneElement.attrib:
                if geneElement.attrib["Symbol"] in target_symbol_values:
                    return(True)
    return(False)



    
def filter_xml_for_gene_symbol(input_fp, output_fp, target_symbol_values,
                               chunk_size=1024):
    in_variation_archive = False
    in_header = True
    for line in input_fp:
        if in_header:
            output_fp.write(line)
        if '<VariationArchive' in line:
            in_variation_archive = True
            in_header = False
            variation_archive = io.StringIO()
        if in_variation_archive:
            variation_archive.write(line)
        if '</VariationArchive>' in line:
            in_variation_archive = False
            xml_string = variation_archive.getvalue()
            variation_archive.close()
            if contains_target_gene_symbol(xml_string,
                                           target_symbol_values):
                output_fp.write(xml_string)


def old(input_fp, output_fp, target_symbol_values,  chunk_size=1024):
    for event, elem in ET.iterparse(input_fp, events=('start', 'end')):
        if event == 'start' and elem.tag == "ClnVarVariationRelease":
            output_FP.write(ET.tostring(elem, encoding='UTF-8'))
        if event == 'end' and elem.tag == 'VariationArchive':
            gene_element = elem.find('.//Gene')
            if gene_element is not None:
                if "Symbol" in gene_element.attrib:
                    if gene_element.attrib["Symbol"] in target_symbol_values:
                        output_fp.write(ET.tostring(elem, encoding='UTF-8'))
            elem.clear()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', "--input")
    parser.add_argument('-o', "--output")
    parser.add_argument('-g', "--gene", action='append')  
    args = parser.parse_args()
    if args.input.endswith('gz'):
        input_fp = io.TextIOWrapper(gzip.open(args.input, 'rb'),
                                    encoding='utf-8')
    else:
        input_fp = open(args.input, 'r')
    output_fp = open(args.output, 'w')
    filter_xml_for_gene_symbol(input_fp, output_fp, args.gene)
    output_fp.write("</ClinVarVariationRelease>".encode('UTF-8'))

if __name__ == "__main__":
    main()


