#!/usr/bin/env python

"""
Description:
    Python script 'lovd2VCF' takes in a LOVD/exLOVD table flat file and converts it to
    vcf format. Used primarily for the purposes of data extraction for integration
    into the ga4gh reference server.
"""


from __future__ import print_function, division

import argparse
import re
import sys
from collections import defaultdict

from common import vcf_files_helper
from common.hgvs_utils import HgvsWrapper
from common.variant_utils import VCFVariant
from hgvs.exceptions import HGVSError
from common.seq_utils import SeqRepoWrapper
import hgvs

def parse_args():
    """
    Description:
        function 'parse_args' parses arguments from command-line and returns an argparse
        object containing the arguments and their values. Default values are 'False' if option
        is not listed in the command, else the option value is set to True.
    """
    parser = argparse.ArgumentParser(description='Convert database table to VCF format.')
    parser.add_argument('-i', '--input', type=argparse.FileType('r'),
                        help='Input file for conversion.')
    parser.add_argument('-a', '--inAnnot', default='/hive/groups/cgl/brca/phase1/data/resources/exLOVDAnnotation',
                        help='Input annotation file for conversion. Tab-delimited with 1st column representing field name and 2nd column representing the field description. Default(/hive/groups/cgl/brca/phase1/data/resources/exLOVDAnnotation)')
    parser.add_argument('-o', '--out', type=argparse.FileType('w'),
                        help='Ouput VCF file result.')
    parser.add_argument('-e', '--errors', type=argparse.FileType('w'),
                        help='File containing all LOVD variants that could not be parsed.')
    parser.add_argument('-s', '--source', help='Source from which data is extracted.')

    options = parser.parse_args()
    return options


def main():
    options = parse_args()
    inputFile = options.input
    annotFile_path = options.inAnnot
    vcfFile = options.out
    errorsFile = options.errors
    source = options.source

    # open and store annotation fields in a dictionary
    annotDict = defaultdict()
    with open(annotFile_path) as inAnnotFile:
        for line in inAnnotFile:
            line = line.strip().split('\t')
            annotDict[line[0]] = line[1]

    # print header lines to vcf file
    print('##fileformat=VCFv4.0', file=vcfFile)
    if source == "exLOVD":
        print('##source=exLOVD', file=vcfFile)
    elif source == "LOVD":
        print('##source=LOVD', file=vcfFile)
    else:
        raise ValueError('Source is %s, must be either LOVD or exLOVD' % (source))
    print('##reference=GRCh37', file=vcfFile)
    for annotation, description in annotDict.items():
        print('##INFO=<ID={0},Number=.,Type=String,Description="{1}">'.format(annotation.replace(' ', '_'), description), file=vcfFile)
    print('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO', file=vcfFile)

    # extract INFO field column indicies for annotation terms
    headerline = inputFile.readline().strip().replace(' ', '_').replace('"', '').split('\t')

    fieldIdxDict = defaultdict()
    for index, field in enumerate(headerline):
        fieldIdxDict[field] = index

    seq_fetcher37 = SeqRepoWrapper(assembly_name=SeqRepoWrapper.ASSEMBLY_NAME_hg37)
    hgvs_wrapper = HgvsWrapper.get_instance()

    # extract info from each line of the flat file
    for line in inputFile:
        line = line.replace('"', '')
        INFO_field = list()
        parsedLine = line.strip().split('\t')
        for field in headerline:
            field_index = fieldIdxDict[field]
            field_value = parsedLine[field_index]
            field_value = vcf_files_helper.normalize_field_value(field_value)
            INFO_field.append('{0}={1}'.format(field, field_value))

        # extract hgvs cDNA term for variant and cleanup formatting
        # Sometimes dna_change is in the field cDNA, sometimes it's labeled dna_change.
        if 'cDNA' in fieldIdxDict:
            hgvsName = parsedLine[fieldIdxDict['cDNA']]
        elif 'dna_change' in fieldIdxDict:
            hgvsName = parsedLine[fieldIdxDict['dna_change']].replace(':', '.3:') # TODO fix hack!!
        else:
            sys.exit("ERROR: could not parse hgvs name.")
        if hgvsName == '-':
            print(parsedLine)
            continue
        queryHgvsName = re.sub(r'[^\x00-\x7F]+', '', hgvsName).rstrip().split(';')[0]
        queryHgvsName = re.sub(r'\s', '', queryHgvsName) # remove whitespaces within name
        INFO_field_string = ';'.join(INFO_field)

        v = None

        try:
            v = vcf_files_helper.cdna_str_to_genomic_var(queryHgvsName,
                                                    HgvsWrapper.GRCh37_Assem,
                                                    hgvs_wrapper, seq_fetcher37)
        except HGVSError as e:
            print('Could not parse cdna field ' + str(
                queryHgvsName) + '. Error was ' + str(e), file=errorsFile)

        if not v and source == "LOVD":
            v = from_genomic(parsedLine, fieldIdxDict, hgvs_wrapper, seq_fetcher37, errorsFile)

        if v:
            print('{0}\t{1}\t{2}\t{3}\t{4}\t.\t.\t{5}'.format(v.chr,
                                                                  v.pos,
                                                                  queryHgvsName,
                                                                  v.ref, v.alt,
                                                                  INFO_field_string),
                  file=vcfFile)
        else:
            print('Could not process line. Skipping variant: ' + str(
                queryHgvsName), file=errorsFile)


def from_genomic(parsedLine, fieldIdxDict, hgvs_wrapper, seq_fetcher37, errorsFile):
    acc = 'NC_0000' + str(parsedLine[fieldIdxDict['chromosome']].replace('chr', '')) + '.10'

    var_str = acc + ":" + parsedLine[fieldIdxDict['gDNA']]

    try:
        var_hgvs = hgvs_wrapper.hgvs_parser.parse(var_str)
        var_hgvs_norm = hgvs.normalizer.Normalizer(hgvs_wrapper.hgvs_dp, shuffle_direction=5).normalize(var_hgvs)
        return VCFVariant.from_hgvs_obj(var_hgvs_norm, seq_fetcher37)
    except HGVSError as e:
        print('Could not parse genomic field ' + str(var_str) + '. Error was ' + str(e), file=errorsFile)

    return None


if __name__ == "__main__":
    main()
