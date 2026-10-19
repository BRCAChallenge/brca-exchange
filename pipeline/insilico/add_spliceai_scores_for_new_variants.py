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
import sys
import pysam

from common.hgvs_utils import HgvsWrapper

HGVS_38_COL = 'pyhgvs_Genomic_Coordinate_38'

# Score against the MANE Select transcripts rather than SpliceAI's packaged
# GENCODE V24 table: SpliceAI N-masks everything outside the annotated
# transcript, so the packaged BRCA1 record (1,334 bases short at the 3' end)
# distorts every variant within ~5 kb of it.  See make_spliceai_annotation.py.
DEFAULT_ANNOTATION = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                  'spliceai_annotations',
                                  'brca_mane_grch38.txt')
SPLICEAI_RUNNER = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                               'spliceai_runner.py')


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
    parser.add_argument("-g", "--annotation", default=DEFAULT_ANNOTATION,
                        help="Gene annotation for spliceAI (-A): 'grch37', "
                             "'grch38', or a path to an annotation table. "
                             "Defaults to the MANE Select table for BRCA1/2.")
    parser.add_argument("-o", "--output_vcf",
                        help="Output VCF file with all variants scored")
    parser.add_argument("-p", "--precision", default="3",
                        help="Decimal places for the spliceAI delta scores")
    parser.add_argument("-s", "--scored_variants_vcf",
                        help="VCF with spliceAI scores for scored variants")
    parser.add_argument('-t', "--temp_dir",
                        help="Temp directory")
    args = parser.parse_args()
    return(args)


def variant_id_string(record):
    return HgvsWrapper.get_instance().genomic_hgvs(
        record.chrom, record.pos, record.ref, record.alts[0])
    

def list_all_variants_in_vcf(this_vcf):
    variants_in_vcf = {}
    reader = pysam.VariantFile(this_vcf)
    for record in reader:
        variants_in_vcf[variant_id_string(record)] = 1
    reader.close()
    return(variants_in_vcf)

def vcf_unscored_variants(all_variants_vcf, scored_variants_vcf,
                          unscored_variants_vcf,
                          max_batch_size):
    variant_count = 0
    scored_variants = list_all_variants_in_vcf(scored_variants_vcf)
    reader = pysam.VariantFile(all_variants_vcf)
    writer = pysam.VariantFile(unscored_variants_vcf, 'w', header=reader.header)
    for record in reader:
        variant = variant_id_string(record)
        if not variant in scored_variants:
            writer.write(record)
            variant_count += 1
            if variant_count >= max_batch_size:
                break
    reader.close()
    writer.close()
    print(variant_count, "records written")
    return(variant_count)

def run_spliceai(unscored_vcf, newly_scored_vcf,
                 genome_fa_file, annotation, depth, precision, debug=True):
    # spliceai_runner.py is the spliceAI CLI with the delta-score precision
    # raised from spliceAI's hard-coded two decimal places.
    spliceai_cmd = [sys.executable, SPLICEAI_RUNNER,
                    "-I", unscored_vcf, "-O", newly_scored_vcf,
                    "-R", genome_fa_file, "-A", annotation, "-D", depth,
                    "--precision", precision]
    if debug:
        print("About to execute", spliceai_cmd)
    subprocess.run(spliceai_cmd, check=True)


def merge_scored_vcf(scored_vcf, newly_scored_vcf, output_vcf, debug=True):
    merge_cmd = ["vcf-concat", scored_vcf, newly_scored_vcf]
    if debug:
        print("About to run", merge_cmd)
    with open(output_vcf, "w") as fp:
        subprocess.run(merge_cmd, stdout=fp, check=True)
    
def main():
    args = parse_args()
    all_variants_scored = False
    print(args.temp_dir)
    scored_vcf = args.temp_dir + "/scored_variants.vcf"
    unscored_vcf = args.temp_dir + "/unscored_variants.vcf"
    newly_scored_vcf = args.temp_dir + "/newly_score_variants.vcf"
    if args.scored_variants_vcf:
        shutil.copy2(args.scored_variants_vcf, scored_vcf)
    else:
        # No previous release to seed from (e.g. a from-scratch run): start
        # from an empty "already scored" VCF, so every variant in
        # all_variants_vcf is treated as unscored and run through spliceAI.
        reader = pysam.VariantFile(args.all_variants_vcf)
        seed = pysam.VariantFile(scored_vcf, 'w', header=reader.header)
        seed.close()
        reader.close()
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
                         args.annotation, args.depth, args.precision)
            merge_scored_vcf(scored_vcf, newly_scored_vcf, args.output_vcf)
            #
            # Copy the merged output file, with the old and new scores, to the
            # temp input scored file, and run the loop again to see if any
            # variants have not yet been scored
            shutil.copy2(args.output_vcf, scored_vcf)
    #
    # scored_vcf holds every variant's score at this point -- unchanged from
    # the previous release if nothing needed scoring, or the last merged
    # result otherwise. Either way, it's the final output.
    shutil.copy2(scored_vcf, args.output_vcf)
    #
    # Cleanup
    os.remove(scored_vcf)
    if os.path.exists(unscored_vcf):
        os.remove(unscored_vcf)
    if os.path.exists(newly_scored_vcf):
        os.remove(newly_scored_vcf)
            
            
if __name__ == "__main__":
    main()
