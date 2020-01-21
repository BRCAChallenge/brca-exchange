#!/usr/bin/env python
"""
this scripts takes the enigma variant list and merge vcf files in a folder into
the exisitng enigma variants:
"""
import argparse
import csv
import logging
import os
import pickle
import re
import subprocess
from copy import deepcopy
from numbers import Number
from shutil import copy

import vcf

import aggregate_reports
from common import seq_utils, config
import utilities
import variant_equivalence
from variant_merging_constants import *

DISCARDED_REPORTS_WRITER = None

def options(parser):
    parser.add_argument("-i", "--input", help="Input VCF directory",
                        default="/home/brca/pipeline-data/pipeline-input/")
    parser.add_argument("-o", "--output",
                        default="/home/brca/pipeline-data/pipeline-output/")
    parser.add_argument("-c", "--config")
    parser.add_argument('-a', "--artifacts_dir", help='Artifacts directory with pipeline artifact files.')
    parser.add_argument("-v", "--verbose", action="count", default=False, help="determines logging")


def main():
    global DISCARDED_REPORTS_WRITER

    parser = argparse.ArgumentParser()
    options(parser)

    args = parser.parse_args()

    gene_config_df = config.load_config(args.config)

    gene_regions_dict = config.extract_gene_regions_dict(gene_config_df, 'start_hg38_legacy_variants', 'end_hg38_legacy_variants')

    gene_regions_trees = seq_utils.build_interval_trees_by_chr(gene_regions_dict.keys(), lambda c,s,e: None)

    genome_regions_symbol_dict = config.get_genome_regions_symbol_dict(gene_config_df, 'start_hg38_legacy_variants', 'end_hg38_legacy_variants')

    seq_provider = seq_utils.SeqRepoWrapper(regions_preload=gene_regions_dict.keys())

    if args.verbose:
        logging_level = logging.DEBUG
    else:
        logging_level = logging.CRITICAL

    log_file_path = args.artifacts_dir + "variant_merging.log"
    logging.basicConfig(filename=log_file_path, filemode="w", level=logging_level,
                        format=' %(asctime)s %(filename)-15s %(message)s')

    discarded_reports_file = open(args.artifacts_dir + "discarded_reports.tsv", "w")

    fieldnames = ['Report_id', 'Source', 'Reason', 'Variant']

    DISCARDED_REPORTS_WRITER = csv.DictWriter(discarded_reports_file, delimiter="\t", fieldnames=fieldnames)
    DISCARDED_REPORTS_WRITER.writeheader()

    # merge repeats within data sources before merging between data sources
    source_dict, columns, variants = preprocessing(args.input, args.output, seq_provider, gene_regions_trees)

    # merges repeats from different data sources, adds necessary columns and data
    print "\n------------merging different datasets------------------------------"
    for source_name, file in source_dict.iteritems():
        (columns, variants) = add_new_source(columns, variants, source_name,
                                             file, FIELD_DICT[source_name], genome_regions_symbol_dict)

    # standardizes genomic coordinates for variants
    print "\n------------standardizing genomic coordinates-------------"
    variants = variant_standardize(columns, seq_provider, gene_regions_trees, variants=variants)

    # compare dna sequence results of variants and merge if equivalent
    print "------------dna sequence comparison merge-------------------------------"
    variants = string_comparison_merge(variants, seq_provider)

    # write final output to file
    write_new_tsv(args.output + "merged.tsv", columns, variants)

    # copy enigma file to artifacts directory along with other ready files
    copy(args.input + ENIGMA_FILE, args.output)

    # write reports to reports file
    aggregate_reports.write_reports_tsv(args.output + "reports.tsv", columns, args.output, genome_regions_symbol_dict)

    discarded_reports_file.close()

    print "final number of variants: %d" % len(variants)
    print "Done"


def variant_standardize(columns, seq_provider, gene_regions_tree, variants="pickle"):
    """standardize variants such
    1. "-" in ref or alt is removed, and a leading base is added, e.g. ->T is changed to N > NT
    2. remove trailing same bases: e.g. AGGGG > TGGGG is changed to A>T
    3. remove leading same baes: e.g. position 100, AAT > AAG is changed to position 102 T>G
    """

    global DISCARDED_REPORTS_WRITER

    # Get indexes of all BX_ID columns by source.
    bx_id_column_indexes = get_bx_id_column_indexes(columns)

    if variants == "pickle":
        with open("temp_variants.pkl", "r") as fv:
            variants = pickle.loads(fv.read())
        fv.close()
    variants_to_remove = list()
    variants_to_add = {}
    for ev, items in variants.iteritems():
        bx_ids_for_variant = get_bx_ids_for_variant(bx_id_column_indexes, items)
        chr = items[COLUMN_VCF_CHR]
        pos = items[COLUMN_VCF_POS]
        ref = items[COLUMN_VCF_REF]
        alt = items[COLUMN_VCF_ALT]
        if ref == "None":
            ref = ""
        if alt == "None":
            alt = ""
        if re.search("^-", ref) or re.search("^-", alt):
            (chr, pos, ref, alt) = add_leading_base(chr, pos, ref, alt, seq_provider)
        if len(ref) < 1 or len(alt) < 1:
            (chr, pos, ref, alt) = add_leading_base(chr, pos, ref, alt, seq_provider)
        (chr, pos, ref, alt) = trim_bases(chr, pos, ref, alt)

        hgvs = "chr%s:g.%s:%s>%s" % (str(chr), str(pos), ref, alt)

        # If the reference is wrong, remove the variant
        reason_for_discard = ""

        if is_outside_boundaries(chr, pos, gene_regions_tree):
            reason_for_discard = "Reference outside boundaries"
        elif not ref_correct(chr, pos, ref, alt, seq_provider):
            reason_for_discard = "Incorrect Reference"
        elif variant_is_false(ref, alt):
            reason_for_discard = "Variant ref and alt are the same"

        if reason_for_discard:
            variants_to_remove = prepare_variant_for_removal_and_log(ev, hgvs, items, bx_ids_for_variant, reason_for_discard, variants_to_remove)
            continue

        items[COLUMN_VCF_POS] = pos
        items[COLUMN_VCF_REF] = ref
        items[COLUMN_VCF_ALT] = alt
        newHgvs = "chr%s:g.%s:%s>%s" % (str(chr), str(pos), ref, str(alt))

        if newHgvs != ev:
            logging.debug("Changed genomic coordinate representation, replacing %s with %s", ev, newHgvs)
            variants_to_remove.append(ev)
            variants_to_add = add_variant_to_dict(variants_to_add, newHgvs, items)

    variants = remove_bad_variants(variants_to_remove, variants)
    variants = add_and_merge_new_variant_representations(variants_to_add, variants)

    return variants


def remove_bad_variants(variants_to_remove, variants):
    for old_variant in variants_to_remove:
        del variants[old_variant]
    return variants


def add_and_merge_new_variant_representations(variants_to_add, variants):
    for genomic_coordinate, values in variants_to_add.iteritems():
        variants = add_variant_to_dict(variants, genomic_coordinate, values)
    return variants


def add_variant_to_dict(variant_dict, genomic_coordinate, values):
    # If the variant is already in the dictionary, merge them together.
    if genomic_coordinate in variant_dict:
        existing_variant = variant_dict[genomic_coordinate]
        equivalent_variant = values
        logging.info("Merging equivalent variants \n %s and \n %s", existing_variant, equivalent_variant)
        assert(len(existing_variant) == len(equivalent_variant))

        # merge properties from equivalent variant into existing variant
        for i, existing_variant_property in enumerate(existing_variant):

            # skip if dealing with chr, pos, ref, or alt since one representation is enough
            if i == COLUMN_VCF_CHR or i == COLUMN_VCF_POS or i == COLUMN_VCF_REF or i == COLUMN_VCF_ALT:
                continue

            # get same property from equivalent variant
            equivalent_variant_property = normalize_values(equivalent_variant[i])
            existing_variant_property = normalize_values(existing_variant_property)
            existing_variant[i] = existing_variant_property

            # move on if they're equal or if equivalent variant property is blank
            if equivalent_variant_property == existing_variant_property or equivalent_variant_property == "-":
                continue

            # if the old value is blank, replace it with the new value
            if existing_variant_property == "-":
                existing_variant[i] = equivalent_variant_property
            else:
                # combine properties into a list
                if type(existing_variant_property) != list:
                    merged_properties = [existing_variant_property]
                else:
                    merged_properties = existing_variant_property

                assert type(merged_properties) == list

                if type(equivalent_variant_property) == list:
                    for prop in equivalent_variant_property:
                        if prop not in merged_properties:
                            merged_properties.append(prop)
                elif equivalent_variant_property not in merged_properties:
                    merged_properties.append(equivalent_variant_property)

                # replace existing data with updates
                existing_variant[i] = merged_properties
                logging.debug("Merged properties: %s", merged_properties)

        logging.debug('Merged output: \n %s', existing_variant)
    else:
        variant_dict[genomic_coordinate] = values

    return variant_dict


def get_bx_ids_for_variant(bx_id_column_indexes, items):
    bx_ids_for_variant = {}
    for key in bx_id_column_indexes.keys():
        bx_ids_for_variant[key] = items[bx_id_column_indexes[key]]
    return bx_ids_for_variant


def get_bx_id_column_indexes(columns):
    bx_id_column_indexes = {}
    for i, column in enumerate(columns):
        if "BX_ID" in column:
            bx_id_column_indexes[column] = i
    return bx_id_column_indexes


def normalize_values(value):
    # standardize data representation by denoting empty as '-' and stripping whitespace off strings
    if value is None or value == "":
        value = DEFAULT_CONTENTS
        return value

    if value == ['-'] or value == []:
        return [DEFAULT_CONTENTS]

    if isinstance(value, int) or isinstance(value, float):
        value = str(value)

    if isinstance(value, basestring):
        value = value.strip()
    else:
        # handle lists
        normalized_values = []
        for v in value:
            if v is None or v == "-" or v == "":
                continue
            else:
                if isinstance(v, basestring):
                    v = v.strip()
                if isinstance(v, int) or isinstance(v, float):
                    v = str(v)
                if v not in normalized_values:
                    normalized_values.append(v)
        value = normalized_values
    return value


def trim_bases(chr, pos, ref, alt):
    #ref, alt = v.split(":")[2].split(">")
    if len(ref) <= 1 or len(alt) <= 1:
        return (chr, pos, ref, alt)
    else:
        (ref, alt) = trim_trailing(ref, alt)
        #v = ":".join(v.split(":")[0:2] + ["{0}>{1}".format(ref, alt)])
        (chr, pos, ref, alt) = trim_leading(chr, pos, ref, alt)
        return (chr, pos, ref, alt)


def trim_trailing(ref, alt):
    if len(ref) <= 1 or len(alt) <= 1:
        return ref, alt
    elif ref[-1] != alt[-1]:
        return ref, alt
    else:
        ref = ref[:-1]
        alt = alt[:-1]
        return trim_trailing(ref, alt)


def trim_leading(chr, pos, ref, alt):
    pos = int(pos)
    if len(ref) == 1 or len(alt) == 1 or ref[0] != alt[0]:
        return (chr, pos, ref, alt)
    else:
        ref = ref[1:]
        alt = alt[1:]
        return trim_leading(chr, str(pos+1), ref, alt)


def add_leading_base(chr, pos, ref, alt, seq_provider):
    pos = int(pos)
    empty_ref = False
    empty_alt = False
    if utilities.isEmpty(ref):
        ref = ""
        empty_ref = True
    if utilities.isEmpty(alt):
        alt = ""
        empty_alt = True

    seq = seq_provider.get_seq(int(chr), pos - 1, 2)
    seq_pos = 1

    if empty_ref is True and empty_alt is True:
        raise Exception("both ref and alt are empty")
    elif empty_ref is True:
        # If the ref is empty, get the base at the position and append it to ref and alt
        leading_base = seq[seq_pos]
        return (chr, str(pos), leading_base + ref, leading_base + alt)
    elif empty_alt is True:
        # If the alt is empty, get the base at the position just before where the deletion happens
        # and append it to the ref and alt
        leading_base = seq[seq_pos - 1]
        return (chr, str(pos - 1), leading_base + ref, leading_base + alt)
    else:
        raise Exception("add leading base called but both ref and alt were provided!")


def variant_is_false(ref, alt):
    # If ref and alt are the same, the variant is considered bad data
    return ref == alt


def string_comparison_merge(variants, seq_wrapper):
    # makes sure the input genomic coordinate strings are unique (no dupes)
    assert (len(variants.keys()) == len(set(variants.keys())))

    logging.info('Calculating all equivalent variants without pickle dump.')
    vcf_variant_dict = { k : VCFVariant(int(v[COLUMN_VCF_CHR]),
                                        int(v[COLUMN_VCF_POS]),
                                        v[COLUMN_VCF_REF],
                                        v[COLUMN_VCF_ALT]) for k, v in variants.iteritems() }

    whole_seq_provider = seq_utils.WholeSeqSeqProvider(seq_wrapper)

    equivalence = variant_equivalence.find_equivalent_variants_whole_seq(vcf_variant_dict, whole_seq_provider)

    n_before_merge = 0
    for each in equivalence:
        n_before_merge += len(each)
    n_after_merge = len(equivalence)
    logging.info('Before merge: %s', str(n_before_merge))
    logging.info('After merge: %s', str(n_after_merge))
    print "%d equivalent variants are merged into %d unique variants" %(
          n_before_merge, n_after_merge)
    for equivalent_v in equivalence:
        #
        # equivalent_v contains a set of variants found to be equivalent.
        # The next step is to merge data for these variants, which will
        # end up in the array merged_row.
        merged_row = []
        for each_v in equivalent_v:
            if len(merged_row) == 0:
                #
                # If this is the first variant in the equivalence set, initialize
                # the merged data to the data for this variant.
                merged_row = variants[each_v]
            else:
                for index, values_merged_so_far in enumerate(merged_row):
                    values_to_be_merged = variants[each_v][index]
                    # If either the new value or the old value is non-blank, use it.
                    if values_merged_so_far == DEFAULT_CONTENTS and values_to_be_merged != DEFAULT_CONTENTS:
                        merged_row[index] = values_to_be_merged
                    elif values_merged_so_far != DEFAULT_CONTENTS and values_to_be_merged == DEFAULT_CONTENTS:
                        merged_row[index] = values_merged_so_far
                        # Skip over the VCF columns.  We're going to assume that one
                        # equivalence of them is enough, which will simplify life for
                        # downstream methods.
                    elif index == COLUMN_VCF_CHR or index == COLUMN_VCF_POS or index == COLUMN_VCF_REF or index == COLUMN_VCF_ALT:
                        continue
                    else:
                        # If both the new value and old value are non-blank and different,
                        # generate a list that contains both new and old values.
                        # If the old value is already a list, append the new value.
                        # If the old value is not a list, create a list containing the old
                        # value and append the new value.
                        if values_merged_so_far != values_to_be_merged:
                            if type(values_merged_so_far) != list:
                                values_merged_so_far = [values_merged_so_far]
                            if values_to_be_merged not in values_merged_so_far:
                                if type(values_to_be_merged) == list:
                                    values_merged_so_far.extend(values_to_be_merged)
                                else:
                                    values_merged_so_far.append(values_to_be_merged)
                            merged_row[index] = values_merged_so_far

            # Remove each variant in the equivalence set from the hash of
            # variants.  Later on, we'll add an entry for the entire equivalence
            # set.
            variants.pop(each_v)
        variants[",".join(list(equivalent_v))] = merged_row
    return variants


def preprocessing(input_dir, output_dir, seq_provider, gene_regions_trees):
    # Preprocessing variants:
    source_dict = {
                   "1000_Genomes": GENOME1K_FILE + "for_pipeline",
                   "ClinVar": CLINVAR_FILE,
                   "LOVD": LOVD_FILE,
                   "exLOVD": EX_LOVD_FILE,
                   "ExAC": EXAC_FILE,
                   "ESP": ESP_FILE,
                   "BIC": BIC_FILE,
                   "GnomAD": GNOMAD_FILE,
                   "Findlay_BRCA1_Ring_Function_Scores": FINDLAY_BRCA1_RING_FUNCTION_SCORES_FIELDS_FILE
                   }
    print "\n" + input_dir + ":"
    print "---------------------------------------------------------"
    print "ENIGMA: {0}".format(ENIGMA_FILE)
    for source_name, file_name in source_dict.iteritems():
        print source_name, ":", file_name
    print "\n------------preprocessing--------------------------------"
    print "remove sample columns and two erroneous rows from 1000 Genome file"
    f_1000G = open(input_dir+ GENOME1K_FILE + "for_pipeline", "w")
    subprocess.call(
       ["bash", "1000g_preprocess.sh", input_dir + GENOME1K_FILE], stdout=f_1000G)

    # merge multiple variant per vcf into multiple lines
    for source_name, file_name in source_dict.iteritems():
        print "convert to one variant per line in ", source_name
        f_in = open(input_dir + file_name, "r")
        f_out = open(output_dir+ source_name + ".vcf", "w")
        # Individual reports (lines in VCF/TSV) are given ids as part of the one_variant_transform method.
        one_variant_transform(f_in, f_out, source_name)
        f_in.close()
        f_out.close()

        print "merge repetitive variants within ", source_name
        f_in = open(output_dir + source_name + ".vcf", "r")
        f_out = open(output_dir + source_name + "ready.vcf", "w")
        repeat_merging(f_in, f_out)
        source_dict[source_name] = f_out.name

    print "-------check if genomic coordinates are correct----------"
    (columns, variants) = save_enigma_to_dict(input_dir + ENIGMA_FILE, output_dir, seq_provider, gene_regions_trees)

    new_source_dict = {}
    for source_name, file_name in source_dict.iteritems():
        f = open(file_name, "r")
        d_wrong = output_dir + "wrong_genome_coors/"
        if not os.path.exists(d_wrong):
            os.makedirs(d_wrong)
        f_wrong = open(output_dir + "wrong_genome_coors/" +
                       source_name + "_wrong_genome_coor.vcf", "w")
        f_right = open(output_dir+ "right" + source_name, "w")

        new_source_dict[source_name] = f_right.name

        vcf_reader = vcf.Reader(f, strict_whitespace=True)
        vcf_wrong_writer = vcf.Writer(f_wrong, vcf_reader)
        vcf_right_writer = vcf.Writer(f_right, vcf_reader)
        n_wrong, n_total = 0, 0
        for record in vcf_reader:
            ref = record.REF.replace("-", "")
            v = [record.CHROM, record.POS, ref, "dummy"]
            if not ref_correct(record.CHROM, record.POS, record.REF, record.ALT, seq_provider) or is_outside_boundaries(record.CHROM, record.POS, gene_regions_trees):
                logging.warning("Reference incorrect for Chrom: %s, Pos: %s, Ref: %s, and Alt: %s",
                                record.CHROM, record.POS, record.REF, record.ALT)
                vcf_wrong_writer.write_record(record)
                n_wrong += 1
            else:
                vcf_right_writer.write_record(record)
            n_total += 1
        f_right.close()
        f_wrong.close()
        print "in {0}, wrong: {1}, total: {2}".format(source_name, n_wrong, n_total)

    return new_source_dict, columns, variants


def repeat_merging(f_in, f_out):
    """takes a vcf file, collapses repetitive variant rows and write out
        to a new vcf file (without header)"""
    vcf_reader = vcf.Reader(f_in, strict_whitespace=True)
    variant_dict = {}  # str -> Record
    num_repeats = 0
    for record in vcf_reader:
        genome_coor = "chr{0}:{1}:{2}>{3}".format(
            record.CHROM, str(record.POS), record.REF, record.ALT[0])
        if genome_coor not in variant_dict.keys():
            variant_dict[genome_coor] = deepcopy(record)
        else:
            num_repeats += 1
            for key in record.INFO:
                if key not in variant_dict[genome_coor].INFO.keys():
                    variant_dict[genome_coor].INFO[key] = deepcopy(record.INFO[key])
                else:
                    new_value = deepcopy(record.INFO[key])
                    old_value = deepcopy(variant_dict[genome_coor].INFO[key])

                    if type(new_value) != list:
                        new_value = [new_value]
                    if type(old_value) != list:
                        old_value = [old_value]

                    # This if statement is crucial to not mess up text fields
                    # containing ',' and hence being treated as separate fields.
                    # The list(set(new_value + old_value)) statement below would
                    # garble it otherwise.
                    if new_value == old_value and key != "individuals":
                        continue
                    else:
                        # FIXME: is there a better name for this? it seems it now only
                        # applies to scv to ensure the order is the same,
                        # but we don't hold this concern for other list fields...
                        if key in LIST_TYPE_FIELDS:
                            merged_value = list(new_value + old_value)
                        # The "individuals" values from LOVD submissions are
                        # added together when merging variants.
                        elif key == "individuals":
                            merged_value = [str(int(new_value[0]) + int(old_value[0]))]
                        else:
                            merged_value = list(set(new_value + old_value))

                        # Remove empty strings from list
                        merged_value = filter(None, merged_value)
                        variant_dict[genome_coor].INFO[key] = deepcopy(merged_value)
    print "number of repeat records: ", num_repeats, "\n"
    vcf_writer = vcf.Writer(f_out, vcf_reader)
    for record in variant_dict.values():
        vcf_writer.write_record(record)
    f_in.close()
    f_out.close()


def get_header(f):
    header = ""
    for line in f:
        if "#" in line:
            header += line
    return header


def one_variant_transform(f_in, f_out, source_name):
    """takes a vcf file, read each row, if the ALT field contains more than
       one item, create multiple variant row based on that row. also adds
       ids to all individual reports (each line in the vcf). writes new vcf"""
    vcf_reader = vcf.Reader(f_in, strict_whitespace=True)
    vcf_writer = vcf.Writer(f_out, vcf_reader)
    count = 1
    for record in vcf_reader:
        n = len(record.ALT)
        if n == 1:
            if source_name == "ExAC":
                record = append_exac_allele_frequencies(record)
            record.INFO['BX_ID'] = count
            count += 1
            vcf_writer.write_record(record)
        else:
            for i in range(n):
                new_record = deepcopy(record)
                new_record.ALT = [deepcopy(record.ALT[i])]
                new_record.INFO['BX_ID'] = count
                count += 1
                for key in record.INFO.keys():
                    value = deepcopy(record.INFO[key])
                    if type(value) == list and len(value) == n:
                        new_record.INFO[key] = [value[i]]
                if source_name == "ExAC":
                    new_record = append_exac_allele_frequencies(record, new_record, i)
                vcf_writer.write_record(new_record)


def append_exac_allele_frequencies(record, new_record=None, i=None):
    if new_record is None:
        for subpopulation in EXAC_SUBPOPULATIONS:
            # calculate allele frequencies for each subpopulation
            allele_count = record.INFO[("AC_" + subpopulation)]
            allele_number = record.INFO[("AN_" + subpopulation)]
            allele_frequency = "-"
            if len(allele_count) > 0 and allele_number != 0:
                allele_frequency = float(allele_count[0]) / float(allele_number)
                allele_frequency = str(utilities.round_sigfigs(allele_frequency, 3))
            record.INFO[("AF_" + subpopulation)] = allele_frequency
        return record
    else:
        new_record.INFO['AF'] = record.INFO['AF'][i]
        for subpopulation in EXAC_SUBPOPULATIONS:
            allele_count = record.INFO[("AC_" + subpopulation)][i]
            allele_number = record.INFO[("AN_" + subpopulation)]
            allele_frequency = "-"
            if allele_number != 0:
                allele_frequency = float(allele_count) / float(allele_number)
                allele_frequency = str(utilities.round_sigfigs(allele_frequency, 3))
            new_record.INFO[("AF_" + subpopulation)] = allele_frequency
        return new_record


def write_new_tsv(filename, columns, variants):
    merged_file = open(filename, "w")
    merged_file.write("\t".join(columns)+"\n")
    for key, variant in sorted(variants.iteritems()):
        if len(variant) != len(columns):
            raise Exception("mismatching number of columns in head and row")
        for ii in range(len(variant)):
            if type(variant[ii]) == list:
                comma_delimited_string = ",".join(str(xx) for xx in variant[ii])
                variant[ii] = comma_delimited_string
            elif type(variant[ii]) == int:
                variant[ii] = str(variant[ii])
        merged_file.write("\t".join(variant)+"\n")
    merged_file.close()


def add_new_source(columns, variants, source, source_file, source_dict, genome_regions_symbol_dict):
    print "adding {0} into merged file.....".format(source)
    old_column_num = len(columns)
    for column_title in source_dict.keys():
        columns.append(column_title+"_{0}".format(source))
    vcf_reader = vcf.Reader(open(source_file, 'r'), strict_whitespace=True)
    overlap = 0
    variants_num = 0
    for record in vcf_reader:
        variants_num += 1
        genome_coor = ("chr" + str(record.CHROM) + ":g." + str(record.POS) + ":" +
                       record.REF + ">" + str(record.ALT[0]))
        if genome_coor in variants.keys():
            overlap += 1
            if type(variants[genome_coor][COLUMN_SOURCE]) != list:
                variants[genome_coor][COLUMN_SOURCE] = [variants[genome_coor][COLUMN_SOURCE]]
            variants[genome_coor][COLUMN_SOURCE].append(source)
        else:
            variants[genome_coor] = associate_chr_pos_ref_alt_with_item(record, old_column_num, source, genome_coor, genome_regions_symbol_dict)
        for value in source_dict.values():
            try:
                variants[genome_coor].append(record.INFO[value])
            except KeyError:
                logging.warning("KeyError appending VCF record.INFO[value] to variant. Variant: %s \n Record.INFO: %s \n value: %s", variants[genome_coor], record.INFO, value)
                if source == "BIC":
                    variants[genome_coor].append(DEFAULT_CONTENTS)
                    logging.debug("Could not find value %s for source %s in variant %s, inserting default content %s instead.", value, source, DEFAULT_CONTENTS)
                else:
                    raise Exception("There was a problem appending a value for %s to variant %s" % (value, variants[genome_coor]))
    # for those enigma record that doesn't have a hit with new genome coordinate
    # add extra cells of "-" to the end of old record
    for value in variants.values():
        if len(value) != len(columns):
            value += [DEFAULT_CONTENTS] * len(source_dict)
    print "number of variants in " + source + " is ", variants_num
    print "overlap with previous dataset: ", overlap
    print "number of total variants with the addition of " + source + " is: ", len(variants), "\n"
    for index, value in variants.iteritems():
        if len(value) != len(columns):
            raise Exception("mismatching number of columns in head and row")
    return (columns, variants)


def associate_chr_pos_ref_alt_with_item(line, column_num, source, genome_coor, genome_regions_symbol_dict):
    # places genomic coordinate data in correct positions to align with relevant columns in output tsv file.
    item = ['-'] * column_num
    item[COLUMN_SOURCE] = source

    item[COLUMN_GENOMIC_HGVS] = genome_coor
    item[COLUMN_VCF_CHR] = line.CHROM
    item[COLUMN_VCF_POS] = line.POS
    item[COLUMN_VCF_REF] = line.REF
    item[COLUMN_VCF_ALT] = str(line.ALT[0])

    chr_tree = genome_regions_symbol_dict.get(int(item[COLUMN_VCF_CHR]))

    if not chr_tree:
        raise Exception(
            "Did find data for chromosome {}".format(int(item[COLUMN_VCF_CHR])))

    symbols = list(chr_tree.at(int(item[COLUMN_VCF_POS])))
    assert len(symbols) == 1, "Expect exactly one symbol at a given position, but got {} for chr {} pos {}".format(len(symbols), item[COLUMN_VCF_CHR], item[COLUMN_VCF_POS])

    _, _, item[COLUMN_GENE] = symbols[0] # don't care about start and end of genomic position

    return item


def associate_chr_pos_ref_alt_with_enigma_item(line):
    # places source and genomic coordinate data in correct positions to align with enigma columns
    items = line.strip().split("\t")
    items.insert(COLUMN_SOURCE, "ENIGMA")
    v = items[COLUMN_GENOMIC_HGVS].replace("-", "").replace("chr", "").replace(">", ":")
    (chrom, pos, ref, alt) = v.split(":")
    items.insert(COLUMN_VCF_CHR, chrom)
    items.insert(COLUMN_VCF_POS, pos)
    items.insert(COLUMN_VCF_REF, ref)
    items.insert(COLUMN_VCF_ALT, alt)
    for ii in range(len(items)):
        if items[ii] is None or items[ii] == '':
            items[ii] = DEFAULT_CONTENTS
    return (items, chrom, pos, ref, alt)


def add_columns_to_enigma_data(line):
    # adds necessary columns to enigma data
    columns = line.strip().split("\t")
    columns = [c + "_ENIGMA" for c in columns if c != "Genomic_Coordinate"]
    columns.insert(COLUMN_SOURCE, "Source")
    columns.insert(COLUMN_GENOMIC_HGVS, "Genomic_Coordinate")
    columns.insert(COLUMN_VCF_CHR, "Chr")
    columns.insert(COLUMN_VCF_POS, "Pos")
    columns.insert(COLUMN_VCF_REF, "Ref")
    columns.insert(COLUMN_VCF_ALT, "Alt")
    return columns


def save_enigma_to_dict(path, output_dir, seq_provider, gene_regions_trees):
    global DISCARDED_REPORTS_WRITER

    enigma_file = open(path, "r")
    variants = dict()
    line_num = 0
    f_wrong = open(output_dir + "ENIGMA_wrong_genome.txt", "w")
    n_wrong, n_total = 0, 0
    bx_id_column_index = None
    for line in enigma_file:
        line_num += 1
        if line_num == 1:
            columns = add_columns_to_enigma_data(line)
            for i, column in enumerate(columns):
                if "BX_ID" in column:
                    bx_id_column_index = i
            f_wrong.write(line)
        else:
            (items, chrom, pos, ref, alt) = associate_chr_pos_ref_alt_with_enigma_item(line)
            bx_id = items[bx_id_column_index]
            hgvs = "chr%s:g.%s:%s>%s" % (str(chrom), str(pos), ref, alt)

            if ref_correct(chrom, pos, ref, alt, seq_provider) and not is_outside_boundaries(chrom, pos, gene_regions_trees):
                variants = add_variant_to_dict(variants, hgvs, items)
            elif pos == 'None':
                logging.warning("Position is none for Enigma report, throwing away: %s", line)
                log_discarded_reports("ENIGMA", bx_id, hgvs, "None position")
                n_wrong += 1
                f_wrong.write(line)
            else:
                logging.warning("Ref incorrect for Enigma report, throwing away: %s", line)
                log_discarded_reports("ENIGMA", bx_id, hgvs, "Incorrect Reference. Is outside Boundaries {}".format(is_outside_boundaries(chrom, pos, gene_regions_trees)))
                n_wrong += 1
                f_wrong.write(line)

            n_total += 1

    f_wrong.close()
    print "in ENIGMA, wrong: {0}, total: {1}".format(n_wrong, n_total)
    return (columns, variants)


def is_outside_boundaries(c, pos, gene_regions_trees):
    c = int(c)
    pos = int(pos)

    if c not in gene_regions_trees.keys():
        raise ValueError("No relevant genes on chromosome {}".format(c))

    chr_regions = gene_regions_trees[c]
    return len(chr_regions.at(pos)) == 0


def ref_correct(chr, pos, ref, alt, seq_provider):
    if pos == "None":
        return False

    pos = int(pos)

    genomeRef = seq_provider.get_seq_at(int(chr), pos, len(ref))

    if (not ref.upper().startswith(genomeRef.upper())):
        logging.warning("genomeref not equal ref for: chr, pos, ref, len(ref), genomeref, len(genomeref), alt: %s, %s, %s, %s, %s, %s, %s",
                        chr, pos, ref.upper(), len(ref), genomeRef.upper(), len(genomeRef), alt)
        return False
    else:
        return True


def prepare_variant_for_removal_and_log(original_hgvs, normalized_hgvs, items, bx_ids_for_variant, reason_for_discard, variants_to_remove):
    if reason_for_discard == "Incorrect Reference":
        logging.warning("Ref incorrect using %s", normalized_hgvs)
        logging.warning("Original variant representation of incorrect ref variant before add_leading_base: %s", str(items))
    elif reason_for_discard == "Variant ref and alt are the same":
        logging.warning("Variant ref and alt are the same for variant %s", normalized_hgvs)
        logging.warning("Original variant representation: %s", str(items))
    else:
        logging.warning("Bad data for variant: %s", normalized_hgvs)
        logging.warning("Original variant representation: %s", str(items))

    for key in bx_ids_for_variant.keys():
        reports = bx_ids_for_variant[key]
        if utilities.isEmpty(reports):
            continue
        else:
            prefix = "BX_ID_"
            source = key[len(prefix):]
            log_discarded_reports(source, reports, normalized_hgvs, reason_for_discard)
    variants_to_remove.append(original_hgvs)
    return variants_to_remove


def log_discarded_reports(source, reports, hgvs, reason):
    # if reports is a list, log each report individually
    if not isinstance(reports, basestring) and not isinstance(reports, Number):
        for report in reports:
            log_discarded_report(source, report, hgvs, reason)
    else:
        report = reports
        log_discarded_report(source, report, hgvs, reason)


def log_discarded_report(source, report, hgvs, reason):
    report = int(report)
    logging.warning("Report discarded: %s \n Source: %s \n Reason for discard: %s \n Variant: %s", report, source, reason, hgvs)
    DISCARDED_REPORTS_WRITER.writerow({'Report_id': report, 'Source': source, 'Reason': reason, 'Variant': hgvs})


if __name__ == "__main__":
    main()
