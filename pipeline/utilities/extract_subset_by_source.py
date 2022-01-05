#!/usr/bin/env python                                                           
'''     
Given a variants_output.tsv file and the corresponding metadata file 
variants_output_field_metadata.tsv, and given also a list of one or more data 
sources to exclude, generate a subset of variants_output.tsv that excludes
the columns specific to the indicated sources, exclusing also any variants that
are specific to those sources, and generates a subset variants_output.tsv with
corresponding variants_output_field_metadata.tsv.
''' 

import argparse
import csv
import collections
import re
import sys


csv.field_size_limit(sys.maxsize)

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--variants_input", help="input variant tsv file",
                        default="bb.tsv")
    parser.add_argument("--metadata_input", help="input metadata file",
                        default="field_metadata.tsv")
    parser.add_argument("--variants_output", help="output variant tsv file",
                        default="variants_output_selected.tsv")
    parser.add_argument("--metadata_output", help="output metadata file",
                        default="variants_output_selected_field_metadata.tsv")
    parser.add_argument("-e", "--exclude", action='append',
                        help="Data source to exclude")
    parser.add_argument('-d', '--debug', default=False)
    args = parser.parse_args()
    return args


def select_output_fields(all_fields, exclude):
    """
    Given the set of all available fields, select the ones to retain for 
    output by excluding any fields for which the suffix matches the name of
    some data source to exclude
    """
    fields_to_remove = list()
    output_fields = all_fields.copy()
    for field in all_fields:
        field = re.sub("(\s)+$", "", field)  
        for source in exclude:
            search_string = source.lower() + "$"
            if re.search(search_string, field.lower()):
                print("excluding field", field, "due to", source)
                fields_to_remove.append(field)
    for field in fields_to_remove:
        if field in output_fields:
            output_fields.remove(field)
    return output_fields


def remove_excluded_sources(source_str, exclude):
    """
    Given the list of input sources for a variant, remove any that are
    on the list to be excluded
    """
    sources_to_remove = list()
    sources = re.split(',', source_str)
    for item in exclude:
        if item in sources:
            sources_to_remove.append(item)
    for item in sources_to_remove:
        sources.remove(item)
    remaining_sources = ','.join(sources)
    return(remaining_sources)

def subset_this_variant(input_variant, exclude, fieldnames,
                        source_colname="Source"):
    """
    For a single variant, generate a subset with the columns to output
    """
    included_sources = remove_excluded_sources(input_variant[source_colname],
                                               exclude)
    if len(included_sources) == 0:
        return(None)
    else:
        output_variant = dict((item,item) for item in fieldnames)
        for item in fieldnames:
            output_variant[item] = input_variant[item]
        output_variant[source_colname] = included_sources
        return(output_variant)


def write_variant_subset(variants_in, variants_out, exclude,
                         source_colname="Source"):
    """
    For each variant in the TSV file, generate the output subset.
    In the source field, remove any excluded sources.  If no sources are
    left, go onto the next line.  Otherwise, remove the columns to be
    excluded in the final output.  
    """
    variants_out.writeheader()
    for variant in variants_in:
        output_variant = subset_this_variant(variant, exclude,
                                             variants_out.fieldnames)
        if output_variant is not None:
            variants_out.writerow(output_variant)

            
def generate_metadata_subset(meta_in, meta_out, exclude):
    """
    Generate a subset of the fields metadata file that omits any fields that
    pertain to an excluded source
    """
    first_line = True
    for line in meta_in:
        if first_line:
            print("writing first line")
            meta_out.write(line)
            first_line = False
        else:
            skip_line = False
            fields = re.split('\t', line)
            for item in exclude:
                search_string = "_" + item.lower() + '$'
                if re.search(search_string, fields[0].lower()):
                    skip_line = True
                    print("skipping", fields[0], "due to", item)
            if not skip_line:
                meta_out.write(line)

def main():
    args = parse_args()
    with open(args.variants_input, 'r') as vi:
        with open(args.variants_output, 'w') as vo:
            variants_in = csv.DictReader(vi, delimiter="\t")
            output_columns = select_output_fields(variants_in.fieldnames,
                                                  args.exclude)
            variants_out = csv.DictWriter(vo, delimiter="\t",
                                          fieldnames=output_columns)
            write_variant_subset(variants_in, variants_out, args.exclude)
    with open(args.metadata_input, 'r') as m_in:
        with open(args.metadata_output, 'w') as m_out:
            generate_metadata_subset(m_in, m_out, args.exclude)


if __name__ == "__main__":
    main()

