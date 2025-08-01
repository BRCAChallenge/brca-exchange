#!/usr/bin/env python
from collections import OrderedDict
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
        string Clinical_significance_ENIGMA;      "Clinical Significance as curated by the ENIGMA VCEP"
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
    color = genomeBrowserUtils.pathogenicityToColor(rec.Clinical_significance_ENIGMA)
    out_url = "https://brcaexchange.org/variant/" + rec.CA_ID
    #
    # When generating the mouseOver, truncate the cDNA and protein HGVS string to 50 characters each, 
    # to not overhwelm the browser's internal limit of 255 characters.
    mouseOver = ("<b>Gene:</b> %s<br>" + \
                 "<b>HGVS cDNA:</b> %s<br>" + \
                 "<b>HGVS Protein:</b> %s<br>" + \
                 "<b>VCEP Curation:</b> %s<br>" + \
                 "<b>URL:</b> %s<br>") \
                 % (rec.Gene_Symbol, rec.pyhgvs_cDNA[0:25], rec.pyhgvs_Protein[0:25],
                    rec.Clinical_significance_ENIGMA,
                    out_url)
    outRow = [chrom, start, end, name, score, strand, thickStart, thickEnd, color, out_url,
              rec.Gene_Symbol,
              genomeBrowserUtils.displayString(rec.pyhgvs_cDNA[0:254]),
              genomeBrowserUtils.displayString(rec.pyhgvs_Protein[0:254]),
              genomeBrowserUtils.displayString(rec.CA_ID),
              rec.Clinical_significance_ENIGMA, mouseOver]
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
        _write_auto_sql_file(args.auto_sql_file)

        for line in ifh:
            row = line.rstrip("\n").rstrip("\r").split("\t")
            rd = OrderedDict(zip(headers, row)) # row as dict

            if int(rd["Hg38_End"]) - int(rd["Hg38_Start"]) < args.length_threshold:
                write_track_item(rec, str(int(rd["pyhgvs_Hg37_Start"])-1), rd["pyhgvs_Hg37_End"], ofhg19v)
                write_track_item(rec, str(int(rd["Hg38_Start"])-1), rd["Hg38_End"], ofhg38v)
            else:
                write_track_item(rec, str(int(rd["pyhgvs_Hg37_Start"])-1), rd["pyhgvs_Hg37_End"], ofhg19sv)
                write_track_item(rec, str(int(rd["Hg38_Start"])-1), rd["Hg38_End"], ofhg38sv)
                

        print("wrote to %s, %s, %s and %s" % (ofhg19v.name, ofhg38v.name, ofhg19sv.name, ofhg38sv.name))


if __name__ == '__main__':
    main()
