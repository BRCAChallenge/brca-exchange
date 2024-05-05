#!/usr/bin/env python
from collections import namedtuple, OrderedDict
import argparse
import html

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
                        default="popfreq.hg19.bed")
    parser.add_argument("-o38", "--output-hg38", help="Output BED file with hg38",
                        default="popfreq.hg38.bed")

    parser.add_argument("-a", "--auto-sql-file", help="Field definitions in AutoSQL format",
                        default="popfreq.as")
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
        string symbol;     "Gene Symbol"
        string cdna_hgvs;       "Variant ID in cDNA HGVS nomenclature"
        string provisional_evidence_code;      "Provisional ACMG code"
        string provional_code_description; "Accompanying description"
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
            if int(rec.Hg38_End) - int(rec.Hg38_Start) < 60:

                pat = rec.Provisional_evidence_code_popfreq
                if "BA1" in pat:
                    color = "253,231,37"  # Standalone benign, yellow
                elif "BS1 " in pat:
                    color = "160,218,57"  # Benign Strong, light green
                elif "BS1_Supporting" in pat:
                    color = "31,161,135"  # Benign supporting, green
                elif "No code met (below threshold)" in pat:
                    color = "54,92,141"   # Inconclusive, blue"
                elif "No code met" in pat and "(indel)" in pat:
                    color = "54,92,141"   # Likewise inconclusive, blue"
                elif "PM2_Supporting" in pat:
                    color = "68,1,84"     # Supporting pathogenic, purple
                else:
                    color = "128,128,128" # Cannot be called, grey

                chrom = "chr"+rec.Chr
                score = 0
                strand = "."
                name = rec.pyhgvs_cDNA[0:254]
                if name == "?":
                    assert(False)

                #description = html.escape(rec.Provisional_evidence_code_description_popfreq)
                description = "Click on the track item for details on the reference population, observed filter allele frequency, datasets evaluated and thresholds"
                mouseOver = (("<b>Provisional ACMG Evidence Code:</b> %s<br>" + \
                              "<b>Details:</b> %s") \
                             % (rec.Provisional_evidence_code_popfreq, description))[:245] + "<br>"
        

                #Start with the hg19 version
                start = str(int(rec.pyhgvs_Hg37_Start)-1)
                end = rec.pyhgvs_Hg37_End
                thickStart = start
                thickEnd = end
                outRow = [chrom, start, end, name, score, strand, thickStart, thickEnd, color, 
                          rec.Gene_Symbol, rec.pyhgvs_cDNA[0:254],
                          rec.Provisional_evidence_code_popfreq,
                          description[:250],
                          mouseOver]

                outRow = [str(x) for x in outRow]
                ofh19.write("\t".join(outRow)+"\n")

                # Repeat with the hg38 version
                ftLen = int(end)-int(start)
                start = str(int(rec.Hg38_Start)-1)
                end = str(int(start)+ftLen)
                thickStart = start
                thickEnd = end
                outRow = [chrom, start, end, name, score, strand, thickStart, thickEnd, color, 
                          rec.Gene_Symbol, rec.pyhgvs_cDNA[0:254],
                          rec.Provisional_evidence_code_popfreq,
                          description[:250],
                          mouseOver]

                outRow = [str(x) for x in outRow]
                ofh38.write("\t".join(outRow)+"\n")

        print("wrote to %s and %s" % (ofh19.name, ofh38.name))


if __name__ == '__main__':
    main()
