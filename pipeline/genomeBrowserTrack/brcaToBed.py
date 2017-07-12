from collections import namedtuple, OrderedDict

def addUrls(s, url=None):
    " transform a list of URLs to hrefs "
    lines = []
    for part in s.split(","):
        part = part.strip()
        if part=="":
            continue
        if part.startswith("http"):
            label = part.split("/")[-1]
            if "=" in label:
                label = label.split("=")[-1]
            part = "<a href='%s'>%s</a>" % (part, label)
            lines.append(part)
        else:
            if url==None:
                lines.append(part)
            else:
                part = "<a href='%s%s'>%s</a>" % (url, part, part)
                lines.append(part)

    return ", ".join(lines)

ifh = open("output/release/built_with_change_types.tsv")

ofh = open("brcaExchange.bed", "w")
asFh = open("brcaExchange.as", "w")

asFh.write("""table brcaExchanges
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
""")

headers = ifh.readline().rstrip("\n").rstrip("\r").strip().split("\t")

rowRec = namedtuple("rec", headers)

skipCols = ["Chr", "Pos", "pyhgvs_Hg37_Start", "pyhgvs_Hg37_End"]

for h in headers:
    if h in skipCols:
        continue
    asFh.write('    lstring %s; "%s"' % (h, h.replace("_", " ")))
    asFh.write('\n')
asFh.write('    string _mouseOver; "mouse over field hidden"\n')
asFh.write('    )')

for line in ifh:
    row = line.rstrip("\n").rstrip("\r").split("\t")
    rec = rowRec(*row)
    rd = OrderedDict(zip(headers, row)) # row as dict

    chrom = "chr"+rec.Chr
    start = str(int(rec.pyhgvs_Hg37_Start)-1)
    end = rec.pyhgvs_Hg37_End

    #name = rec.Protein_Change
    #if name=="-":
        #name = rec.pyhgvs_Protein.split(".")[-1].strip("(").strip(")")
    #if name=="-":
        #name = rec.HGVS_cDNA
    #if name=="?" or name=="=":
    name = rec.pyhgvs_cDNA.split(":")[-1]
    if name=="?":
        assert(False)
    
    #score = len(rec.Source.split(","))
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
        if h in skipCols:
            continue
        val = rd[h]
        if val=="-":
            val = ""
        val = val.strip().strip(",").strip()

        if h=="Source":
            val = val.replace(",", ", ")
        if h=="Assertion_method_citation_ENIGMA" and val!="":
            val = "<a href='%s'>%s</a>" % (val, val.split("/")[-1])
        if h=="Source_URL":
            val = addUrls(val)
        if h=="SCV_ClinVar":
            val = addUrls(val, url="https://www.ncbi.nlm.nih.gov/clinvar/?term=")
        if h=="Submitter_ClinVar":
            val = val.replace("_", " ")
        if "," in val:
            val = val.replace(",", ", ")

        outRow.append( val )

    outRow = [str(x) for x in outRow]

    mouseOvers = []
    if rec.Pathogenicity_all!="":
        mouseOvers.append(rec.Pathogenicity_all)
    if rec.pyhgvs_cDNA!="-":
        mouseOvers.append(rec.pyhgvs_cDNA)
    if rec.Discordant!="Concordant":
        mouseOvers.append(rec.Discordant)
    mouseOver = ", ".join(mouseOvers)
    outRow.append(mouseOver)

    ofh.write("\t".join(outRow)+"\n")

print "wrote to %s and %s" % (ofh, asFh)

