#!/usr/bin/env python
"""
From ESP, extract data for the BRCA1 and BRCA2 regions.  Either produce data
in the following format:

chrom_pos_ref_alt Freq

e.g. 13_32890339_G_C0.000199681

where freq is one of EA or AA (European Ancestry or African-American Ancestry)

or if the --full option is given, echo the full VCF record.
"""
import argparse
import vcf

EMPTY = '-'


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("inputVcf")
    parser.add_argument("-s", "--start")
    parser.add_argument("-e", "--end")
    parser.add_argument("-a", "--ancestry")
    parser.add_argument("-f", "--full")
    parser.add_argument("-o", "--output")
    args = parser.parse_args()

    start = int(args.start)
    end = int(args.end)
    reader = vcf.Reader(open(args.inputVcf, 'r'))
    if args.full:
        writer = vcf.Writer(open(args.output, "w"), reader)
    for record in reader:
        if record.INFO.has_key("GRCh38_POSITION"):
            tokens = record.INFO["GRCh38_POSITION"][0].split(":")
            if len(tokens) > 1:
                chrom = tokens[0]
                pos = tokens[1]
                if int(pos) >= start and int(pos) <= end:
                    if args.full:
                        record.CHROM = chrom
                        record.POS = pos
                        (eaAlleleFrequency, aaAlleleFrequency, alleleFrequency) = breakUpESPAlleleFrequencies(record.INFO["MAF"])
                        record.INFO['BX_EAAF'] = eaAlleleFrequency
                        record.INFO['BX_AAAF'] = aaAlleleFrequency
                        record.INFO['BX_AF'] = alleleFrequency
                        writer.write_record(record)
                    else:
                        if args.ancestry == "EA":
                            maf = record.INFO["MAF"][0]
                        elif args.ancestry == "AA":
                            maf = record.INFO["MAF"][1]
                        for alt in record.ALT:
                            print "%s_%s_%s_%s %s" % (record.CHROM,
                                                      record.POS,
                                                      record.REF,
                                                      alt, maf)


def breakUpESPAlleleFrequencies(mafArray):
    eaAlleleFrequency = EMPTY
    aaAlleleFrequency = EMPTY
    alleleFrequency = EMPTY
    if len(mafArray) > 2:
        alleleFrequency = "%s" % (float(mafArray[2]) / 100)
    if len(mafArray) > 1:
        aaAlleleFrequency = "%s" % (float(mafArray[1]) / 100)
    if len(mafArray) > 0:
        eaAlleleFrequency = "%s" % (float(mafArray[0]) / 100)
    return (eaAlleleFrequency, aaAlleleFrequency, alleleFrequency)


if __name__ == "__main__":
    # execute only if run as a script
    main()
