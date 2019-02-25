#!/usr/bin/env python

"""
Description:
    Takes in a functional assay table and converts it to vcf format.
"""

from __future__ import print_function, division
import argparse
import sys
import os
from collections import defaultdict
import pyhgvs as hgvs
import pyhgvs.utils as hgvs_utils
from pygr.seqdb import SequenceFileDB
import logging


def parse_args():
    parser = argparse.ArgumentParser(description='Convert database table to VCF format.')
    parser.add_argument('-i', '--input', type=argparse.FileType('r'),
                        help='Input file for conversion.')
    parser.add_argument('-a', '--inAnnot', default='/hive/groups/cgl/brca/phase1/data/resources/exLOVDAnnotation',
                        help='Input annotation file for conversion. Tab-delimited with 1st column representing field name and 2nd column representing the field description. Default(/hive/groups/cgl/brca/phase1/data/resources/exLOVDAnnotation)')
    parser.add_argument('-o', '--out', type=argparse.FileType('w'),
                        help='Ouput VCF file result.')
    parser.add_argument('-l', '--logfile', default='/tmp/functional_assays_to_vcf.log')
    parser.add_argument('-g', '--gpath', default='/hive/groups/cgl/brca/phase1/data/resources/hg19.fa',
                        help='Whole path to genome file. Default: (/hive/groups/cgl/brca/phase1/data/resources/hg19.fa)')
    parser.add_argument('-r', '--rpath', default='/hive/groups/cgl/brca/phase1/data/resources/refseq_annotation.hg19.gp',
                        help='Whole path to refSeq file. Default: (/hive/groups/cgl/brca/phase1/data/resources/refseq_annotation.hg19.gp)')
    parser.add_argument('-s', '--source', default='FunctionalAssay')
    parser.add_argument('-v', '--verbose', action='count', default=False, help='determines logging')

    options = parser.parse_args()
    return options


def main():
    options = parse_args()
    inputFile = options.input
    annotFile_path = options.inAnnot
    vcfFile = options.out
    genome_path = options.gpath
    refseq_path = options.rpath
    source = options.source
    logfile = options.logfile

    if options.verbose:
        logging_level = logging.DEBUG
    else:
        logging_level = logging.CRITICAL

    logging.basicConfig(filename=logfile, filemode="w", level=logging_level)

    with open(refseq_path) as infile:
        transcripts = hgvs_utils.read_transcripts(infile)

    genome = SequenceFileDB(genome_path)

    def get_transcript(name):
        return transcripts.get(name)

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
            field_value = normalize(field, field_value)
            INFO_field.append('{0}={1}'.format(field, field_value))

        # extract hgvs cDNA term for variant and cleanup formatting
        hgvsName = parsedLine[fieldIdxDict['hgvs_nucleotide']]
        if hgvsName == '-':
            logging.debug("hgvs name == '-' for line: %s", parsedLine)
            continue
        gene_symbol = parsedLine[fieldIdxDict['gene_symbol']].lower()
        if gene_symbol == 'brca1':
            transcript = 'NM_007294.3'
        elif gene_symbol == 'brca2':
            transcript = 'NM_000059.3'
        else:
            logging.debug("improper gene symbol: %s", gene_symbol)
            continue
        queryHgvsName = transcript + ':' + hgvsName.rstrip().split(';')[0]
        INFO_field_string = ';'.join(INFO_field)
        try:
            chrom, offset, ref, alt = hgvs.parse_hgvs_name(queryHgvsName, genome, get_transcript=get_transcript)
            chrom = chrom.replace('chr', '')
            print('{0}\t{1}\t{2}\t{3}\t{4}\t.\t.\t{5}'.format(chrom, offset, queryHgvsName, ref, alt, INFO_field_string), file=vcfFile)
        except Exception as e:
            logging.debug("could not parse hgvs field: %s", queryHgvsName)


def normalize(field, field_value):
    if not is_empty(field_value):
        if field_value[0] == ';':
            field_value = field_value[1:]
        if field_value[-1] == ';':
            field_value = field_value[:-1]
        if ';' in field_value:
            field_value = field_value.replace(';', '')
    return field_value


def is_empty(field_value):
    return field_value == '' or field_value is None


if __name__ == "__main__":
    main()
