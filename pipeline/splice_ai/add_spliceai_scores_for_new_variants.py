#!/usr/bin/env python
"""add_spliceai_scores_for_new_variants: identify the variants which are new
since the last release, run spliceAI on those variants only, and merge these
new scores with the scores from the previous release to generate the new
spliceAI scores.  

Note these new scores may contain scores on variants which have been
deleted.  In practice, this is not an issue, since the scores for deleted
variants will be ignored in the merge step.
"""

import argparse
import os
import shutil
import subprocess
import vcf


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-a", "--all_variants_vcf",
                        help="VCF with all variants to be scored")
    parser.add_argument("-b", "--batch_size", default="1000",
                        help="Max batch size for spliceAI runs")
    parser.add_argument("-d", "--depth",
                        help="Search depth for spliceAI")
    parser.add_argument("-f", "--genome_fa_file",
                        help="Pathname for the genome fa file (e.g. hg38.fa)")
    parser.add_argument("-g", "--genome_name",
                        help="Name of the genome")
    parser.add_argument("-o", "--output_vcf",
                        help="Output VCF file with all variants scored")
    parser.add_argument("-s", "--scored_variants_vcf",
                        help="VCF with spliceAI scores for scored variants")
    parser.add_argument('-t', "--temp_dir",
                        help="Temp directory")
    args = parser.parse_args()
    return(args)


def variant_id_string(record):
    variant = "%s:%s:%s:%s" % (record.CHROM, str(record.POS),
                               record.REF, record.ALT)
    return(variant)
    

def list_all_variants_in_vcf(this_vcf):
    variants_in_vcf = {}
    reader = vcf.Reader(open(this_vcf, "r"))
    for record in reader:
        variants_in_vcf[variant_id_string(record)] = 1
    return(variants_in_vcf)

def vcf_unscored_variants(all_variants_vcf, scored_variants_vcf,
                          unscored_variants_vcf,
                          max_batch_size):
    variant_count = 0
    scored_variants = list_all_variants_in_vcf(scored_variants_vcf)
    with open(all_variants_vcf, "r") as rfp:
        reader = vcf.Reader(rfp) 
        with open(unscored_variants_vcf, "w") as wfp:
            writer = vcf.Writer(wfp, reader)
            for record in reader:
                variant = variant_id_string(record)
                if not variant in scored_variants:
                    writer.write_record(record)
                    variant_count += 1
                    if variant_count >= max_batch_size:
                        break
    print(variant_count, "records written")
    return(variant_count)

def run_spliceai(unscored_vcf, newly_scored_vcf,
                 genome_fa_file, genome_name, depth, debug=True):
    spliceai_cmd = ["spliceai", "-I", unscored_vcf, "-O", newly_scored_vcf,
                    "-R", genome_fa_file, "-A", genome_name, "-D", depth]
    if debug:
        print("About to execute", spliceai_cmd)
    subprocess.run(spliceai_cmd)


def merge_scored_vcf(scored_vcf, newly_scored_vcf, output_vcf, debug=True):
    merge_cmd = ["vcf-concat", scored_vcf, newly_scored_vcf]
    if debug:
        print("About to run", merge_cmd)
    with open(output_vcf, "w") as fp:
        subprocess.run(merge_cmd, stdout=fp)
    
def main():
    args = parse_args()
    all_variants_scored = False
    print(args.temp_dir)
    scored_vcf = args.temp_dir + "/scored_variants.vcf"
    unscored_vcf = args.temp_dir + "/unscored_variants.vcf"
    newly_scored_vcf = args.temp_dir + "/newly_score_varaints.vcf"
    shutil.copy2(args.scored_variants_vcf, scored_vcf)
    while not all_variants_scored:
        #
        # Build a VCF file of unscored variants up to the max batch size,
        # run spliceAI on them, and merge their records with those of
        # the previously scored variants
        unscored_variant_count = vcf_unscored_variants(args.all_variants_vcf,
                                                       scored_vcf,
                                                       unscored_vcf,
                                                       int(args.batch_size))
        if unscored_variant_count == 0:
            all_variants_scored = True
        else:
            run_spliceai(unscored_vcf, newly_scored_vcf, args.genome_fa_file,
                         args.genome_name, args.depth)
            merge_scored_vcf(scored_vcf, newly_scored_vcf, args.output_vcf)
            #
            # Copy the merged output file, with the old and new scores, to the
            # temp input scored file, and run the loop again to see if any
            # variants have not yet been scored
            shutil.copy2(args.output_vcf, scored_vcf)
    #
    # Cleanup
    os.remove(scored_vcf)
    os.remove(unscored_vcf)
    os.remove(newly_scored_vcf)
            
            
if __name__ == "__main__":
    main()
