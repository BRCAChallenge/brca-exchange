#!/usr/bin/env python

import argparse
import csv
import os
import re
import sys
# Ideally, I'd import BRCA1 and BRCA2 from variant_merging.py instead of copying it,
# but when I tried this, somehow it messed up my argument parsing.
#from variant_merging import BRCA1, BRCA2  


def variantToAllele(variant, gene, length):
    """Given a variant in genomic HGVS representation, return the alternative allele string,
    consisting of the ALT portion of the HGVS identifier plus flanking bases"""
    (chrom, rawPosition, ref, alt) = re.split("[:|>]", variant)
    if re.search("\.", rawPosition):
        position = rawPosition.split(".")[1]
    else:
        position = rawPosition
    relativeStart = int(position) - gene["hg38"]["start"] - 1
    prefixStart = relativeStart - length
    prefix = gene["hg38"]["sequence"][prefixStart:(prefixStart + length)]
    suffixStart = relativeStart
    suffixLength = len(ref) + length
    suffix = gene["hg38"]["sequence"][suffixStart:suffixStart+suffixLength]
    if not re.search("^"+ref, suffix):
        sys.exit("Error: reference allele %s from variant %s not found in variant %s|%s" \
                     % (ref, variant, prefix, suffix))
    if ref == "-" or ref == ".":
        ref = ""
    if alt == "-" or alt == ".":
        alt = ""
    alleleString = prefix + re.sub("^" + ref, alt, suffix)
    return alleleString

def allelesForOutputVariants(dataset, reference, length, brca1, brca2):
    """Given a tab-delimited output dataset, return the set of alleles represented by the 
    output variants.  Return a dictionary in which the keys are the allele strings, and the
    values are the corresponding rows of the output dataset

   sample command line:
   testForVariantSet.py --input ENIGMA_combined_hg38.tsv --format tsv --dataset built.tsv \
        --reference brca/pipeline-resources/ --length 70
    """
    alleleData = dict()
    for row in csv.DictReader(open(dataset, "r"), delimiter='\t'):
        if row["Gene_Symbol"] == "BRCA1":
            geneThisRow = brca1
        elif row["Gene_Symbol"] == "BRCA2":
            geneThisRow = brca2
        else:
            sys.exit("Unexpected gene in the output file: %s" % row["Gene_Symbol"])
        variants = row["pyhgvs_Genomic_Coordinate_38"]
        for thisVariant in variants.split(","):
            allele = variantToAllele(thisVariant, geneThisRow, length)
            alleleData[allele] = row
    return alleleData
    

def main():
    argparser = argparse.ArgumentParser()
    argparser.add_argument("-i", "--input", help="Input data (variants we're testing for)")
    argparser.add_argument("-f", "--format", help="input file format: tsv (default) or vcf",
                        default="tsv")
    argparser.add_argument("-d", "--dataset", help="Merged dataset")
    argparser.add_argument('-r', "--reference", help="reference data directory",
                        default="/home/brca/pipeline-data/pipeline-resources/")
    argparser.add_argument("-l", "--length", help="Length flanking seqeunce", 
                           type=int, default=70)
    argparser.add_argument("-v", "--verbose", action="count", default=False, 
                        help="determines logging")
    args = argparser.parse_args()

    #
    # This is copied verbatim from variant_merging.py.  I don't like doing this, but
    # my more sophisticated attempts to import the data left me with messed-up
    # argument parsing.
    BRCA1 = {"hg38": {"start": 43000000,
                      "sequence": open(args.reference + "brca1_hg38.txt", "r").read()},
             "hg19": {"start": 41100000,
                      "sequence": open(args.reference + "brca1_hg19.txt", "r").read()}}
    BRCA2 = {"hg38": {"start": 32300000,
                      "sequence": open(args.reference + "brca2_hg38.txt", "r").read()},
             "hg19": {"start": 32800000,
                      "sequence": open(args.reference + "brca2_hg19.txt", "r").read()}}
    
    outputVariantData = allelesForOutputVariants(args.dataset, args.reference, args.length+5, 
                                                 BRCA1, BRCA2)
    if args.format == "tsv":
        inputVariantSet = csv.DictReader(open(args.input, "r"), delimiter='\t')
        for row in inputVariantSet:
            variant = row["Genomic_Coordinate"]
            if not re.search("None", variant):
                if row["Gene_symbol"] == "BRCA1":
                    allele = variantToAllele(variant, BRCA1, args.length)
                elif row["Gene_symbol"] == "BRCA2":
                    allele = variantToAllele(variant, BRCA2, args.length)
                else:
                    sys.exit("unregocnized gene %s in the input data" % row["Gene_symbol"])
                alleleFound = False
                for outputAllele in list(outputVariantData.keys()):
                    if re.search(allele, outputAllele):
                        alleleFound = True
                        matchingVariantData = outputVariantData[outputAllele]
                        print("Variant %s found in output variant(s) %s from source(s) %s" % \
                            (variant, matchingVariantData["pyhgvs_Genomic_Coordinate_38"],
                             matchingVariantData["Source"]))
                if not alleleFound:       
                    print("Error: variant %s with allele %s not found" % (variant, allele))


if __name__ == "__main__":
    main()
