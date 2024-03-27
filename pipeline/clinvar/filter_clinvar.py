#!/usr/bin/env python

import argparse
import gzip
import xml.etree.ElementTree as ET

def filter_xml_for_gene_symbol(input_fp, output_fp, target_symbol_values,
                               chunk_size=1024):
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
    output_fp.write("</ClinVarVariationRelease>".encode('UTF-8'))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', "--input")
    parser.add_argument('-o', "--output")
    parser.add_argument('-g', "--gene", action='append')  
    args = parser.parse_args()
    if args.input.endswith('gz'):
        input_fp = gzip.open(args.input)
    else:
        input_fp = open(args.input, 'r')
    output_fp = open(args.output, 'wb')
    filter_xml_for_gene_symbol(input_fp, output_fp, args.gene)

if __name__ == "__main__":
    main()


