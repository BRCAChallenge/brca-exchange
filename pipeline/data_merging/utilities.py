
from math import floor, log10
import re
from data_merging.variant_merging_constants import *

def round_sigfigs(num, sig_figs):
    if num != 0:
        return round(num, -int(floor(log10(abs(num))) - (sig_figs - 1)))
    else:
        return 0  # Can't take the log of 0


def isEmpty(value):
    return value == '-' or value is None or value == [] or value == ['-'] or value == ''

def is_outside_boundaries(c, pos, gene_regions_trees):
    # Return a flag indicating if the chromosome contains none  of the
    # gene regions of interest, or if the chromosome contains a region
    # of interest but the variant falls outside of it.
    c = int(re.sub("^chr", "", c))
    pos = int(pos)
    if c not in gene_regions_trees.keys():
        return(True)
    chr_regions = gene_regions_trees[c]
    return len(chr_regions.at(pos)) == 0

def associate_chr_pos_ref_alt_with_item(line, column_num, source, genome_coor, genome_regions_symbol_dict):
    # Given a VCF record, initialize the output dictionary with default values, and columns
    # from VCF fields which are needed for all variants
    item = ['-'] * column_num
    item[COLUMN_SOURCE] = source
    item[COLUMN_GENOMIC_HGVS] = genome_coor
    item[COLUMN_VCF_CHR] = line.CHROM
    item[COLUMN_VCF_POS] = line.POS
    item[COLUMN_VCF_REF] = line.REF
    item[COLUMN_VCF_ALT] = str(line.ALT[0])
    symbol = chrom_pos_to_symbol(int(re.sub("^chr", "", line.CHROM)), int(line.POS),
                                      genome_regions_symbol_dict)
    item[COLUMN_GENE] = symbol
    return item


def associate_chr_pos_ref_alt_with_enigma_item(line):
    # Initialize a variant dict from an ENIGMA record
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

def chrom_pos_to_symbol(chrom, pos, genome_regions_symbol_dict):
    # Given a coordinate, return the gene symbol
    chr_tree = genome_regions_symbol_dict.get(int(chrom))
    if not chr_tree:
        raise Exception(
            "In gene symbol dict, did find data for chromosome {}".format(chrom))
    symbols = list(chr_tree.at(int(pos)))
    assert len(symbols) == 1, "Expect exactly one symbol at a given position, but got {} for chr {} pos {}".format(len(symbols), chrom, pos)
    (start_pos, end_pos, gene_symbol) = symbols[0]
    return(gene_symbol)

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

