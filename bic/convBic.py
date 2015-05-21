# convert from NHGRI BIC dump format to VCF
import logging

def convFile(fname, chrom, printHead, rows):
    skipCount = 0
    for line in open(fname):
        fs = line.rstrip("\n").split("\t")
        if line.startswith("Accession"):
            # print out VCF header
            fs = [f.replace("#", "") for f in fs]
            fieldIds = [f.replace(" ", "") for f in fs]
            fieldIds = [f.replace("(", "") for f in fieldIds]
            fieldIds = [f.replace(")", "") for f in fieldIds]
            fieldNames = fs
            if not printHead:
                continue
            print "##fileformat=VCFv4.0"
            print "##source=BIC"
            print "##reference=hg19"

            for fieldId, fieldName in zip(fieldIds, fieldNames):
                print '##INFO=<ID=%s,Number=1,Type=String,Description="%s">' % (fieldId, fieldName)
            print "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
            continue

        acc = fs[0]
        hgvs = fs[9]
        #ex g.41277277G>A
        if not ">" in hgvs:
            skipCount += 1
            continue
        hgvs = hgvs.replace("g.","")
        part1, altAll = hgvs.split(">")
        refAll = part1[-1]
        pos = part1[:-1]
        infos = []
        for fieldId, field in zip(fieldIds, fs):
            fieldId = fieldId.replace(";", "")
            field = field.replace(";", "")
            infos.append("%s=%s" % (fieldId, field.replace(" ", "_")))
        row = [int(chrom), int(pos), acc, refAll, altAll, ".", ".", ";".join(infos)]
        rows.append(row)
    logging.error("%d skipped, cannot import indels yet" % skipCount)

        #acc, ex, nt, codon, baseChange, aaChange, Designation, hgvsCdna, hgvsProt, hgvsHg19, genotype, dbSnp, mutType, clinImp, category, evidence, depositor

def main():
    rows = []
    convFile("/hive/groups/cgl/brca/phase1/bic/brca1_data.txt", "17", True, rows)
    convFile("/hive/groups/cgl/brca/phase1/bic/brca2_data.txt", "13", False, rows)
    rows.sort()
    for row in rows:
        row = [str(x) for x in row]
        print "\t".join(row)

main()
