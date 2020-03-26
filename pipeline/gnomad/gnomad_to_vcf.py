#!/usr/bin/env python

"""
Description:
    Takes in a gnomad table and converts it to vcf format.
"""

from __future__ import print_function, division
import argparse
from collections import defaultdict
import logging

from common import vcf_files_helper

def parse_args():
    parser = argparse.ArgumentParser(description='Convert database table to VCF format.')
    parser.add_argument('-i', '--input', type=argparse.FileType('r'),
                        help='Input file for conversion.')
    parser.add_argument('-a', '--inAnnot', default='/hive/groups/cgl/brca/phase1/data/resources/gnomadAnnotation',
                        help='Input annotation file for conversion. Tab-delimited with 1st column representing field name and 2nd column representing the field description. Default(/hive/groups/cgl/brca/phase1/data/resources/exLOVDAnnotation)')
    parser.add_argument('-o', '--out', type=argparse.FileType('w'),
                        help='Ouput VCF file result.')
    parser.add_argument('-l', '--logfile', default='/tmp/gnomad_to_vcf.log')
    parser.add_argument('-s', '--source', default='gnomAD')
    parser.add_argument('-v', '--verbose', action='count', default=False, help='determines logging')

    options = parser.parse_args()
    return options


def main():
    options = parse_args()
    inputFile = options.input
    annotFile_path = options.inAnnot
    vcfFile = options.out
    source = options.source
    logfile = options.logfile

    if options.verbose:
        logging_level = logging.DEBUG
    else:
        logging_level = logging.CRITICAL

    logging.basicConfig(filename=logfile, filemode="w", level=logging_level)

    # open and store annotation fields in a dictionary
    annotDict = defaultdict()
    with open(annotFile_path) as inAnnotFile:
        for line in inAnnotFile:
            line = line.strip().split('\t')
            annotDict[line[0]] = line[1]

    # print header lines to vcf file
    print('##fileformat=VCFv4.0', file=vcfFile)
    print('##source={0}'.format(source), file=vcfFile)
    print('##reference=GRCh37', file=vcfFile)
    for annotation, description in annotDict.items():
        print('##INFO=<ID={0},Number=.,Type=String,Description="{1}">'.format(annotation.replace(' ', '_'), description), file=vcfFile)
    print('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO', file=vcfFile)

    # extract INFO field column indicies for annotation terms
    headerline = inputFile.readline().strip().replace(' ', '_').replace('"', '').split('\t')

    fieldIdxDict = defaultdict()
    for index, field in enumerate(headerline):
        fieldIdxDict[field] = index

    # extract info from each line of the flat file
    for line in inputFile:
        line = line.replace('"', '')
        INFO_field = list()
        parsedLine = line.strip().split('\t')
        for field in headerline:
            field_index = fieldIdxDict[field]
            field_value = parsedLine[field_index]

            field_value = vcf_files_helper.normalize_field_value(field_value)
            if any(s in field for s in ['_ac', '_an']):
                field_value = field_value.replace('.0', '')

            INFO_field.append('{0}={1}'.format(field, field_value))

        # extract hgvs cDNA term for variant and cleanup formatting
        hgvsName = parsedLine[fieldIdxDict['hgvs']]
        if hgvsName == '-':
            logging.debug("hgvs name == '-' for line: %s", parsedLine)
            continue
        chrom = parsedLine[fieldIdxDict['chrom']].lower()
        # TODO improve this heuristics using config?
        if chrom == '17':
            transcript = 'NM_007294.3'
        elif chrom == '13':
            transcript = 'NM_000059.3'
        else:
            logging.debug("improper gene symbol: %s", chrom)
            continue
        queryHgvsName = transcript + ':' + hgvsName.rstrip().split(';')[0]

        INFO_field_string = ';'.join(INFO_field)

        print('{0}\t{1}\t{2}\t{3}\t{4}\t.\t.\t{5}'.format(parsedLine[fieldIdxDict['chrom']],
                                                          parsedLine[fieldIdxDict['pos']],
                                                          queryHgvsName,
                                                          parsedLine[fieldIdxDict['ref']],
                                                          parsedLine[fieldIdxDict['alt']],
                                                          INFO_field_string),
              file=vcfFile)



if __name__ == "__main__":
    main()
