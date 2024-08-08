#!/usr/bin/env python

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
        if '<VariationArchive' in line:
            in_variation_archive = True
            in_header = False
            variation_archive = io.StringIO()
        if in_header:
            output_fp.write(line)
        if in_variation_archive:
            variation_archive.write(line)
        if '</VariationArchive>' in line:
            in_variation_archive = False
            xml_string = variation_archive.getvalue()
            variation_archive.close()
            if contains_target_gene_symbol(xml_string,
                                           target_symbol_values):
                output_fp.write(xml_string)




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
    output_fp.write("</ClinVarVariationRelease>")

if __name__ == "__main__":
    main()


