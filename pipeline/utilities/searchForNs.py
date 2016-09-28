"""
parse vcf output files for instances of n's

NOTE: requires that all source files to be checked
are in a single directory which must be named
and passed in as input_dir

"""

import vcf
import argparse
import pdb

SOURCE_FILES = {
                '10KG13': 'chr13_brca2_1000g_GRCh37.vcf',
                '10KG17': 'chr17_brca1_1000g_GRCh37.vcf',
                'BIC': 'bic_brca12.sorted.hg38.vcf',
                'CLINVAR': 'ClinVarBrca.vcf',
                'ENIGMA': 'ENIGMA_last_updated.tsv',
                'ESP': 'esp.brca12.sorted.hg38.vcf',
                'EXAC': 'exac.brca12.sorted.hg38.vcf',
                'EXLOVD': 'exLOVD_brca12.sorted.hg38.vcf',
                'SHAREDLOVD': 'sharedLOVD_brca12.sorted.hg38.vcf'
               }

SOURCE_FILE_LINES_BEFORE_VARIANTS = {
                'BIC': 44,
                'CLINVAR': 15,
                'ESP': 34,
                'EXLOVD': 25,
                'SHAREDLOVD': 25
               }


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("--input_dir",
                        help="Directory with vcf_files")

    args = parser.parse_args()

    f_out = open('ns.txt', "w")
    for source in SOURCE_FILES:
        search_source_for_ns(args, source, f_out)

    f_out.close()


def skip_to_variants(f_in, source):
    for i in range(SOURCE_FILE_LINES_BEFORE_VARIANTS[source]):
        f_in.next()


def search_source_for_ns(args, source, f_out):
    if source == '10KG13':
        ns = search_10KG13(args, source)
        f_out.write(ns)
    elif source == '10KG17':
        ns = search_10KG17(args, source)
        f_out.write(ns)
    elif source == 'BIC':
        ns = search_BIC(args, source)
        f_out.write(ns)
    elif source == 'CLINVAR':
        ns = search_CLINVAR(args, source)
        f_out.write(ns)
    elif source == 'ENIGMA':
        ns = search_ENIGMA(args, source)
        f_out.write(ns)
    elif source == 'ESP':
        ns = search_ESP(args, source)
        f_out.write(ns)
    elif source == 'EXAC':
        ns = search_EXAC(args, source)
        f_out.write(ns)
    elif source == 'EXLOVD':
        ns = search_EXLOVD(args, source)
        f_out.write(ns)
    elif source == 'SHAREDLOVD':
        ns = search_SHAREDLOVD(args, source)
        f_out.write(ns)


def has_n(input_to_check):
    if not isinstance(input_to_check, basestring):
        input_to_check = str(input_to_check[0])
    input_to_check = input_to_check.lower()

    allowed_ns = ['ins', 'none', 'variant', 'protein', 'coding', 'nonsense_mediated_decay', 'deletion', 'hgnc', 'ensp',
                  'intron', 'ensg', 'tran', 'down', 'gene', 'n4bp', 'snv', 'uncertain', 'significance', 'benign', 'not',
                  'nm']

    for allowed_word in allowed_ns:
        if allowed_word in input_to_check:
            input_to_check = input_to_check.replace(allowed_word, '')

    if "n" in input_to_check:
        print "**************"
        print input_to_check
        print "**************"
        return input_to_check


def check_for_n_with_vcf_reader(args, source):
    ns = ''
    vcf_reader = vcf.Reader(open(args.input_dir + SOURCE_FILES[source], 'r'))
    for variant in vcf_reader:
        if has_n(variant.REF) or has_n(variant.ALT):
            print 'found'
            print variant
            line = str(source) + ': %s' % (variant)
            ns += (line + '\n')

    return ns


def search_10KG13(args, source):
    print "Searching 10kg13..."
    ns = check_for_n_with_vcf_reader(args, source)
    return ns


def search_10KG17(args, source):
    print "Searching 10kg17..."
    ns = check_for_n_with_vcf_reader(args, source)
    return ns


def search_BIC(args, source):
    print "Searching BIC..."
    ns = ''
    f_in = open(args.input_dir + SOURCE_FILES[source], 'r')
    skip_to_variants(f_in, source)
    current_line_num = SOURCE_FILE_LINES_BEFORE_VARIANTS[source] + 1
    for variant in f_in:
        print current_line_num

        first_hgvs = variant.split()[2]
        HGVS_cDNA = variant.split('HGVS_cDNA=')[1].split(';')[0]
        HGVS_Genomic_hg19 = variant.split('HGVS_Genomic_(hg19)')[1].split(';')[0]

        if has_n(first_hgvs) or has_n(HGVS_cDNA) or has_n(HGVS_Genomic_hg19):
            print 'found'
            print variant
            line = str(source) + ': %s' % (variant)
            ns += (line + '\n')

        current_line_num += 1

    f_in.close()
    return ns


def search_CLINVAR(args, source):
    print "Searching CLINVAR..."
    ns = ''
    f_in = open(args.input_dir + SOURCE_FILES[source], 'r')
    skip_to_variants(f_in, source)
    current_line_num = SOURCE_FILE_LINES_BEFORE_VARIANTS[source] + 1
    for variant in f_in:
        print current_line_num

        HGVS = variant.split('HGVS=')[1].split(';')[0][1:]
        Genomic_Coordinate = variant.split('Genomic_Coordinate=')[1].split(';')[0]
        if has_n(HGVS) or has_n(Genomic_Coordinate):
            print 'found'
            print variant
            line = str(source) + ': %s' % (variant)
            ns += (line + '\n')

        current_line_num += 1

    f_in.close()

    ns += (check_for_n_with_vcf_reader(args, source) + '\n')
    return ns


def search_ESP(args, source):
    print "Searching ESP..."
    ns = ''
    f_in = open(args.input_dir + SOURCE_FILES[source], 'r')
    skip_to_variants(f_in, source)
    current_line_num = SOURCE_FILE_LINES_BEFORE_VARIANTS[source] + 1
    for variant in f_in:
        print current_line_num

        HGVS_cDNA = variant.split('HGVS_CDNA_VAR=')[1].split(';')[0][1:]

        if has_n(HGVS_cDNA):
            print 'found'
            print variant
            line = str(source) + ': %s' % (variant)
            ns += (line + '\n')

        current_line_num += 1

    f_in.close()
    ns += (check_for_n_with_vcf_reader(args, source) + '\n')

    return ns


def search_ENIGMA(args, source):
    print "Searching ENIGMA..."
    ns = ''

    f_in = open(args.input_dir + SOURCE_FILES[source], 'r')
    columns = f_in.next().split('\t')

    for variant in f_in:
        props = variant.split('\t')
        Genomic_Coordinate = props[columns.index('Genomic_Coordinate')]
        Reference_sequence = props[columns.index('Reference_sequence')]
        HGVS_cDNA = props[columns.index('HGVS_cDNA')]
        BIC_Nomenclature = props[columns.index('BIC_Nomenclature')]
        if has_n(Genomic_Coordinate) or has_n(Reference_sequence) or has_n(HGVS_cDNA) or has_n(BIC_Nomenclature):
            print 'found'
            print variant
            line = str(source) + ': %s' % (variant)
            ns += (line + '\n')

    f_in.close()
    return ns


def search_EXAC(args, source):
    print "Searching EXAC..."
    ns = ''
    ns += (check_for_n_with_vcf_reader(args, source) + '\n')
    return ns


def search_EXLOVD(args, source):
    print "Searching EXLOVD..."
    ns = ''
    f_in = open(args.input_dir + SOURCE_FILES[source], 'r')
    skip_to_variants(f_in, source)
    current_line_num = SOURCE_FILE_LINES_BEFORE_VARIANTS[source] + 1
    for variant in f_in:
        print current_line_num

        first_hgvs = variant.split('NM_')[1].split()[0][1:]
        dna_change = variant.split('dna_change=')[1].split(';')[0][1:]
        bic_dna_change = variant.split('bic_dna_change=')[1].split(';')[0]
        if has_n(first_hgvs) or has_n(dna_change) or has_n(bic_dna_change):
            pdb.set_trace()
            print 'found'
            print variant
            line = str(source) + ': %s' % (variant)
            ns += (line + '\n')

        current_line_num += 1

    f_in.close()
    return ns


def search_SHAREDLOVD(args, source):
    print "Searching SHAREDLOVD..."
    ns = ''
    f_in = open(args.input_dir + SOURCE_FILES[source], 'r')
    skip_to_variants(f_in, source)
    current_line_num = SOURCE_FILE_LINES_BEFORE_VARIANTS[source] + 1
    for variant in f_in:
        print current_line_num

        first_hgvs = variant.split('NM_')[1].split()[0][1:]
        dna_change = variant.split('dna_change=')[1].split(';')[0][1:]
        dna_change_genomic = variant.split('dna_change_genomic=')[1].split(';')[0]
        if has_n(first_hgvs) or has_n(dna_change) or has_n(dna_change_genomic):
            pdb.set_trace()
            print 'found'
            print variant
            line = str(source) + ': %s' % (variant)
            ns += (line + '\n')

        current_line_num += 1

    f_in.close()
    return ns


if __name__ == "__main__":
    main()
