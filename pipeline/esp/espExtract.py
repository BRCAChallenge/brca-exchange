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
import pysam

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
    reader = pysam.VariantFile(args.inputVcf)
    if args.full:
        writer = pysam.VariantFile(args.output, 'w', header=reader.header)
    for record in reader:
        if "GRCh38_POSITION" in record.info:
            tokens = record.info["GRCh38_POSITION"][0].split(":")
            if len(tokens) > 1:
                chrom = tokens[0]
                pos = tokens[1]
                if int(pos) >= start and int(pos) <= end:
                    if args.full:
                        record.chrom = chrom
                        record.pos = pos
                        (eaAlleleFrequency, aaAlleleFrequency, alleleFrequency) = breakUpESPAlleleFrequencies(record.info["MAF"])
                        record.info['BX_EAAF'] = eaAlleleFrequency
                        record.info['BX_AAAF'] = aaAlleleFrequency
                        record.info['BX_AF'] = alleleFrequency
                        writer.write(record)
                    else:
                        if args.ancestry == "EA":
                            maf = record.info["MAF"][0]
                        elif args.ancestry == "AA":
                            maf = record.info["MAF"][1]
                        for alt in record.alts:
                            print("%s_%s_%s_%s %s" % (record.chrom,
                                                      record.pos,
                                                      record.ref,
                                                      alt, maf))
    reader.close()
    if args.full:
        writer.close()


def breakUpESPAlleleFrequencies(mafArray):
    eaAlleleFrequency = EMPTY
    aaAlleleFrequency = EMPTY
    alleleFrequency = EMPTY
    if len(mafArray) > 2:
        alleleFrequency = f"{float(mafArray[2])/100:.{6}}"
    if len(mafArray) > 1:
        aaAlleleFrequency = f"{float(mafArray[1])/100:.{6}}"
    if len(mafArray) > 0:
        eaAlleleFrequency = f"{float(mafArray[0])/100:.{6}}"
    return (eaAlleleFrequency, aaAlleleFrequency, alleleFrequency)


if __name__ == "__main__":
    # execute only if run as a script
    main()
