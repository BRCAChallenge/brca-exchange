#!/usr/bin/env python
from collections import namedtuple, OrderedDict
import html
import genomeBrowserUtils




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
        string protein_hgvs;    "Variant ID in protein HGVS nomenclature"
        string CA_ID;       "ClinGen Allele Registry ID"
        string provisional_evidence_code;      "Provisional ACMG code"
        string provional_code_description; "Accompanying description"
        string _mouseOver; "mouse over field hidden"
        )
        """

        asFh.write(sql)
        print("wrote as file to {}".format(asFh.name))

        
def write_track_item(rec, start, end, output_fp):
    chrom = "chr"+rec.Chr
    score = 0
    strand = "."
    name = rec.pyhgvs_cDNA[0:254]
    if name == "?":
        assert(False)
    thickStart = start
    thickEnd = end
    acmgCode = rec.Provisional_Evidence_Code_Popfreq
    color = genomeBrowserUtils.acmgCodeToColor(acmgCode)
    out_url = "https://brcaexchange.org/variant/" + rec.CA_ID
    #                                                                                                 
    # When generating the mouseOver, truncate the strings to 50 characters each,
    # to not overhwelm the browser's internal limit of 255 characters.                                
    description = "Click on the track item for more details"
    mouseOver = (("<b>Provisional ACMG Evidence Code:</b> %s<br>" + \
                  "<b>Details:</b> %s") \
                 % (acmgCode, description))[:245] + "<br>"
    outRow = [chrom, start, end, name, score, strand, thickStart, thickEnd, color, out_url,
              rec.Gene_Symbol,
              genomeBrowserUtils.displayString(rec.pyhgvs_cDNA[0:254]),
              genomeBrowserUtils.displayString(rec.pyhgvs_Protein[0:254]),
              genomeBrowserUtils.displayString(rec.CA_ID),
              acmgCode, description[:254],
              mouseOver]
    outRow = [str(x) for x in outRow]
    output_fp.write("\t".join(outRow)+"\n")

    
def main():
    args = genomeBrowserUtils._get_args()

    with open(args.input, 'r') as ifh:
        ofhg19v = open(args.output_hg19_var, 'w')
        ofhg38v = open(args.output_hg38_var, 'w')
        ofhg19sv = open(args.output_hg19_sv, 'w')
        ofhg38sv =  open(args.output_hg38_sv, 'w')
        print("Reading %s..." % ifh.name)

        headers = ifh.readline().rstrip("\n").rstrip("\r").strip().split("\t")
        rowRec = namedtuple("rec", headers)

        _write_auto_sql_file(args.auto_sql_file)

        for line in ifh:
            row = line.rstrip("\n").rstrip("\r").split("\t")
            rec = rowRec(*row)
            rd = OrderedDict(zip(headers, row)) # row as dict
            if int(rec.Hg38_End) - int(rec.Hg38_Start) < args.length_threshold:
                write_track_item(rec, str(int(rec.pyhgvs_Hg37_Start)-1), rec.pyhgvs_Hg37_End, ofhg19v)
                write_track_item(rec, str(int(rec.Hg38_Start)-1), rec.Hg38_End, ofhg38v)
            else:
                write_track_item(rec, str(int(rec.pyhgvs_Hg37_Start)-1), rec.pyhgvs_Hg37_End, ofhg19sv)
                write_track_item(rec, str(int(rec.Hg38_Start)-1), rec.Hg38_End, ofhg38sv)

        print("wrote to %s and %s" % (ofhg19v.name, ofhg38v.name))


if __name__ == '__main__':
    main()
