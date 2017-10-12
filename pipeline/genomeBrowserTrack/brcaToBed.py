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


def _write_auto_sql_file(as_path, headers, skip_cols):
    with open(as_path, "w") as asFh:

        sql = """table brcaExchanges
        " bed9 + many additional fields from brcaExchange.org "
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
            string outLink;     "Link to BRCAExchange.org"
        """

        asFh.write(sql)

        for h in headers:
            if h in skip_cols:
                continue
            asFh.write('    lstring %s; "%s"' % (h, h.replace("_", " ")))
            asFh.write('\n')
        asFh.write('    string _mouseOver; "mouse over field hidden"\n')
        asFh.write('    )')

        print("wrote as file to {}".format(asFh.name))


def main():
    parser = _get_parser()

    args = parser.parse_args()

    with open(args.input, 'r') as ifh, open(args.output_hg19, 'w') as ofh, open(args.output_hg38, 'w') as ofh38:
        print("Reading %s..." % ifh.name)

        headers = ifh.readline().rstrip("\n").rstrip("\r").strip().split("\t")
        rowRec = namedtuple("rec", headers)
        skip_cols = ["Chr", "Pos", "pyhgvs_Hg37_Start", "pyhgvs_Hg37_End"]

        _write_auto_sql_file(args.auto_sql_file, headers, skip_cols)

        for line in ifh:
            row = line.rstrip("\n").rstrip("\r").split("\t")
            rec = rowRec(*row)
            rd = OrderedDict(zip(headers, row)) # row as dict

            chrom = "chr"+rec.Chr
            start = str(int(rec.pyhgvs_Hg37_Start)-1)
            end = rec.pyhgvs_Hg37_End
            name = rec.pyhgvs_cDNA.split(":")[-1]
            if name == "?":
                assert(False)
            score = 0
            strand = "."
            thickStart = start
            thickEnd = end

            pat = rec.Pathogenicity_all.lower()
            if "pathogen" in pat:
                color = "255,0,0"
            elif "benign" in pat:
                color = "0,0,255"
            elif "uncertain" in pat:
                color = "100,100,100"
            else:
                color = "0,0,0"

            # special case requested by Melissa
            if "discordant" in rec.Discordant.lower():
                color = "255,160,0"

            outLinkId = rec.pyhgvs_Genomic_Coordinate_38
            if len(outLinkId)>255:
                outLinkId = ""
            outRow = [chrom, start, end, name, score, strand, thickStart, thickEnd, color, outLinkId]

            for h in headers:
                if h in skip_cols:
                    continue
                val = rd[h]
                if val == "-":
                    val = ""
                val = val.strip().strip(",").strip()

                if h == "Source":
                    val = val.replace(",", ", ")
                if h == "Assertion_method_citation_ENIGMA" and val != "":
                    val = "<a href='%s'>%s</a>" % (val, val.split("/")[-1])
                if h == "Source_URL":
                    val = _add_urls(val)
                if h == "SCV_ClinVar":
                    val = _add_urls(val, url="https://www.ncbi.nlm.nih.gov/clinvar/?term=")
                if h == "Submitter_ClinVar":
                    val = val.replace("_", " ")
                if "," in val:
                    val = val.replace(",", ", ")

                outRow.append(val)

            outRow = [str(x) for x in outRow]

            mouseOvers = []
            if rec.Pathogenicity_all != "":
                mouseOvers.append(rec.Pathogenicity_all)
            if rec.pyhgvs_cDNA != "-":
                mouseOvers.append(rec.pyhgvs_cDNA)
            if rec.Discordant != "Concordant":
                mouseOvers.append(rec.Discordant)
            mouseOver = ", ".join(mouseOvers)
            outRow.append(mouseOver)

            ofh.write("\t".join(outRow)+"\n")

            # write out a the hg38 version of this line
            ftLen = int(end)-int(start)
            start = str(int(rec.Hg38_Start)-1)
            end = str(int(start)+ftLen)
            thickStart = start
            thickEnd = end
            outRow[1] = start
            outRow[2] = end
            outRow[6] = thickStart
            outRow[7] = thickEnd

            ofh38.write("\t".join(outRow)+"\n")

        print("wrote to %s and %s" % (ofh.name, ofh38.name))


if __name__ == '__main__':
    main()
