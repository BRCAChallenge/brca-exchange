#!/usr/bin/env python
from collections import namedtuple, OrderedDict
import argparse


def _add_urls(s, url=None):
    """ transform a list of URLs to hrefs """
    lines = []
    for part in s.split(","):
        part = part.strip()
        if part == "":
            continue
        if part.startswith("http"):
            label = part.split("/")[-1]
            if "=" in label:
                label = label.split("=")[-1]
            part = "<a href='%s'>%s</a>" % (part, label)
            lines.append(part)
        else:
            if url == None:
                lines.append(part)
            else:
                part = "<a href='%s%s'>%s</a>" % (url, part, part)
                lines.append(part)

    return ", ".join(lines)


def _get_parser():
    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--input", help="Path to built_with_change_types.tsv file",
                        default="output/release/built_with_change_types.tsv")

    parser.add_argument("-o19", "--output-hg19", help="Output BED file with hg19",
                        default="brcaExchange.hg19.bed")
    parser.add_argument("-o38", "--output-hg38", help="Output BED file with hg38",
                        default="brcaExchange.hg38.bed")

    parser.add_argument("-a", "--auto-sql-file", help="Field definitions in AutoSQL format",
                        default="brcaExchange.as")
    return parser


def _write_auto_sql_file(as_path):
    with open(as_path, "w") as asFh:

        sql = """table brcaExchanges
        " These data are in BigBed bed9 format, and include selected fields from https://brcaexchange.org"
        (
        string chrom;      "Chromosome (or contig, scaffold, etc.)"
        uint   chromStart; "Start position in chromosome"
        uint   chromEnd;   "End position in chromosome"
        string name;       "Name of item"
        uint   score;      "Score from 0-1000"
        char[1] strand;    "+ or -"
        uint thickStart;   "Start of where display should be thick (start codon)"
        uint thickEnd;     "End of where display should be thick (stop codon)"
        uint reserved;     "Used as itemRgb as of 2004-11-22"        
        string outlink;    "Link to the variant in BRCA Exchange"
        string symbol;     "Gene Symbol"
        string cdna_hgvs;       "Variant ID in cDNA HGVS nomenclature"
        string CA_ID;       "ClinGen Allele Registry ID"
        string Clinical_significance_ENIGMA;      "Clinical Significance as curated by the ENIGMA VCEP"
        string _mouseOver; "mouse over field hidden"
        )
        """

        asFh.write(sql)


        print("wrote as file to {}".format(asFh.name))


def main():
    parser = _get_parser()

    args = parser.parse_args()

    with open(args.input, 'r') as ifh, open(args.output_hg19, 'w') as ofh19, open(args.output_hg38, 'w') as ofh38:
        print("Reading %s..." % ifh.name)

        headers = ifh.readline().rstrip("\n").rstrip("\r").strip().split("\t")
        rowRec = namedtuple("rec", headers)
        include_cols = ["Chr", "Pos", "pyhgvs_Hg37_Start", "pyhgvs_Hg37_End"]

        _write_auto_sql_file(args.auto_sql_file)

        for line in ifh:
            row = line.rstrip("\n").rstrip("\r").split("\t")
            rec = rowRec(*row)
            rd = OrderedDict(zip(headers, row)) # row as dict

            pat = rec.Clinical_significance_ENIGMA.lower()
            if "pathogen" in pat:
                color = "255,0,0"
            elif "benign" in pat:
                color = "0,255,0"
            elif "uncertain" in pat:
                color = "100,100,100"
            else:
                color = "0,0,0"
            out_url = "https://brcaexchange.org/variant/" + rec.CA_ID

            chrom = "chr"+rec.Chr
            score = 0
            strand = "."
            name = rec.pyhgvs_cDNA[0:254]
            if name == "?":
                assert(False)
            #
            # When generating the mouseOver, truncate the HGVS string to 100 characters, to not overhwelm
            # the browser's internal limit of 255 characters.
            mouseOver = ("<b>Variant ID:</b> %s %s<br>" + \
                         "<b>ENIGMA VCEP Clinical Significance:</b> %s<br>" + \
                         "<b>Variant URL:</b> %s<br>") \
                         % (rec.Gene_Symbol, rec.pyhgvs_cDNA[0:100], rec.Clinical_significance_ENIGMA,
                            out_url)

            #Start with the hg19 version
            start = str(int(rec.pyhgvs_Hg37_Start)-1)
            end = rec.pyhgvs_Hg37_End
            thickStart = start
            thickEnd = end
            outRow = [chrom, start, end, name, score, strand, thickStart, thickEnd, color, out_url,
                      rec.Gene_Symbol, rec.pyhgvs_cDNA[0:254], rec.CA_ID,
                      rec.Clinical_significance_ENIGMA, mouseOver]

            outRow = [str(x) for x in outRow]
            ofh19.write("\t".join(outRow)+"\n")

            # Repeat with the hg38 version
            ftLen = int(end)-int(start)
            start = str(int(rec.Hg38_Start)-1)
            end = str(int(start)+ftLen)
            thickStart = start
            thickEnd = end
            outRow = [chrom, start, end, name, score, strand, thickStart, thickEnd, color, out_url,
                      rec.Gene_Symbol, rec.pyhgvs_cDNA[0:254], rec.CA_ID,
                      rec.Clinical_significance_ENIGMA, mouseOver]

            outRow = [str(x) for x in outRow]
            ofh38.write("\t".join(outRow)+"\n")

        print("wrote to %s and %s" % (ofh19.name, ofh38.name))


if __name__ == '__main__':
    main()
