#!/usr/bin/env python

"""
Charles Markello
6/23/2015

Description:
    Python script 'exLOVD2vcf' takes in a ex-LOVD table flat file and converts it to
    vcf format. Used primarily for the purposes of data extraction for integration
    into the ga4gh reference server.

    Version 1, the basic flat representation of sample data
"""


from __future__ import print_function, division
import argparse
import sys
import os
from collections import defaultdict
import pyhgvs as hgvs
import pyhgvs.utils as hgvs_utils
from pygr.seqdb import SequenceFileDB
import urllib


LOVD_LIST_FIELDS = ["genetic_origin", "RNA", "variant_effect", "individuals", "submitters", "Protein", "cDNA"]


def parse_args():
    """
    Description:
        function 'parse_args' parses arguments from command-line and returns an argparse
        object containing the arguments and their values. Default values are 'False' if option
        is not listed in the command, else the option value is set to True.
    """
    parser = argparse.ArgumentParser(description='Convert exLOVD database table to VCF format.')
    parser.add_argument('-i', '--inEXLOVD', type=argparse.FileType('r'),
                        help='Input exLOVD file for conversion.')
    parser.add_argument('-a', '--inAnnot', default='/hive/groups/cgl/brca/phase1/data/resources/exLOVDAnnotation',
                        help='Input exLOVD annotation file for conversion. Tab-delimited with 1st column representing field name and 2nd column representing the field description. Default(/hive/groups/cgl/brca/phase1/data/resources/exLOVDAnnotation)')
    parser.add_argument('-o', '--out', type=argparse.FileType('w'),
                        help='Ouput VCF file result.')
    parser.add_argument('-e', '--errors', type=argparse.FileType('w'),
                        help='File containing all LOVD variants that could not be parsed.')
    parser.add_argument('-g', '--gpath', default='/hive/groups/cgl/brca/phase1/data/resources/hg19.fa',
                        help='Whole path to genome file. Default: (/hive/groups/cgl/brca/phase1/data/resources/hg19.fa)')
    parser.add_argument('-r', '--rpath', default='/hive/groups/cgl/brca/phase1/data/resources/refseq_annotation.hg19.gp',
                        help='Whole path to refSeq file. Default: (/hive/groups/cgl/brca/phase1/data/resources/refseq_annotation.hg19.gp)')

    options = parser.parse_args()
    return options


def main(args):
    options = parse_args()
    exLOVDFile = options.inEXLOVD
    annotFile_path = options.inAnnot
    vcfFile = options.out
    genome_path = options.gpath
    refseq_path = options.rpath
    errorsFile = options.errors

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
    print('##source=exLOVD', file=vcfFile)
    print('##reference=GRCh37', file=vcfFile)
    for annotation, description in annotDict.items():
        print('##INFO=<ID={0},Number=.,Type=String,Description="{1}">'.format(annotation.replace(' ', '_'), description), file=vcfFile)
    print('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO', file=vcfFile)

    # extract INFO field column indicies for annotation terms
    headerline = exLOVDFile.readline().strip().replace(' ', '_').replace('"', '').split('\t')

    fieldIdxDict = defaultdict()
    for index, field in enumerate(headerline):
        fieldIdxDict[field] = index

    # extract info from each line of the bic flat file
    for line in exLOVDFile:
        line = line.replace('"', '')
        INFO_field = list()
        parsedLine = line.strip().split('\t')
        for field in headerline:
            field_index = fieldIdxDict[field]
            field_value = parsedLine[field_index]
            field_value = normalize(field, field_value)
            INFO_field.append('{0}={1}'.format(field, field_value))

        # extract hgvs cDNA term for variant and cleanup formatting
        # Sometimes dna_change is in the field cDNA, sometimes it's labeled dna_change.
        if 'cDNA' in fieldIdxDict:
            hgvsName = parsedLine[fieldIdxDict['cDNA']]
        elif 'dna_change' in fieldIdxDict:
            hgvsName = parsedLine[fieldIdxDict['dna_change']]
        else:
            sys.exit("ERROR: could not parse hgvs name.")
        if hgvsName == '-':
            print(parsedLine)
            continue
        queryHgvsName = hgvsName.rstrip().split(';')[0]
        INFO_field_string = ';'.join(INFO_field)
        try:
            chrom, offset, ref, alt = hgvs.parse_hgvs_name(queryHgvsName, genome, get_transcript=get_transcript)
            chrom = chrom.replace('chr', '')
            print('{0}\t{1}\t{2}\t{3}\t{4}\t.\t.\t{5}'.format(chrom, offset, queryHgvsName, ref, alt, INFO_field_string), file=vcfFile)
        except Exception as e:
            print(str(e)+': could not parse hgvs field '+queryHgvsName, file=errorsFile)


def normalize(field, field_value):
    if not is_empty(field_value):
        if field_value[0] == ';':
            field_value = field_value[1:]
        if field_value[-1] == ';':
            field_value = field_value[:-1]
        if field not in LOVD_LIST_FIELDS and ';' in field_value:
            field_value = field_value.replace(';', '')
        if field in LOVD_LIST_FIELDS and ';' in field_value:
            # Semicolons are sometimes used as a list delimiter,
            # this changes them to commas for consistency with other fields.
            field_value = field_value.replace(';', ', ')
        if field in LOVD_LIST_FIELDS:
            # Use url encoding to prevent issues in VCF file format.
            # Decoded during merging and reports aggregation.
            field_value = urllib.quote_plus(field_value)
    return field_value


def is_empty(field_value):
    return field_value == '' or field_value is None


if __name__ == "__main__":
    sys.exit(main(sys.argv))
