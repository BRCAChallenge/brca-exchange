
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
    item[COLUMN_VCF_CHR] = str(line.chrom)
    item[COLUMN_VCF_POS] = str(line.pos)
    item[COLUMN_VCF_REF] = line.ref
    item[COLUMN_VCF_ALT] = str(line.alts[0])
    symbol = chrom_pos_to_symbol(int(re.sub("^chr", "", str(line.chrom))), int(line.pos),
                                      genome_regions_symbol_dict)
    item[COLUMN_GENE] = symbol
    return item


def associate_chr_pos_ref_alt_with_enigma_vcf_record(record, info_field_names):
    # Initialize a variant list from an ENIGMA VCF record.
    # Returns the same structure as the old TSV-based function so that
    # all downstream code (variant_merging, aggregate_reports) is unchanged.
    chrom = str(record.chrom)
    pos   = str(record.pos)   # 1-based VCF position
    ref   = record.ref
    alt   = str(record.alts[0])
    genomic_coor = "chr{}:{}:{}>{}".format(chrom, pos, ref, alt)

    # Build the INFO-field portion of the items list in field-declaration order.
    info_values = []
    for field in info_field_names:
        val = record.info.get(field)
        if val is None:
            info_values.append(DEFAULT_CONTENTS)
        elif isinstance(val, tuple):
            info_values.append(
                ','.join(str(v) if v is not None else DEFAULT_CONTENTS for v in val))
        else:
            s = str(val)
            info_values.append(s if s not in ('', 'None') else DEFAULT_CONTENTS)

    items = info_values
    items.insert(COLUMN_SOURCE,      "ENIGMA")
    items.insert(COLUMN_GENOMIC_HGVS, genomic_coor)
    items.insert(COLUMN_VCF_CHR,     chrom)
    items.insert(COLUMN_VCF_POS,     pos)
    items.insert(COLUMN_VCF_REF,     ref)
    items.insert(COLUMN_VCF_ALT,     alt)
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

def add_columns_to_enigma_data_vcf(info_field_names):
    # Build the ENIGMA column list from the VCF INFO field names.
    # Produces the same column structure as the old TSV-based function so that
    # all downstream consumers (variant_merging, aggregate_reports) are unchanged.
    columns = [f + "_ENIGMA" for f in info_field_names]
    columns.insert(COLUMN_SOURCE,      "Source")
    columns.insert(COLUMN_GENOMIC_HGVS, "Genomic_Coordinate")
    columns.insert(COLUMN_VCF_CHR,     "Chr")
    columns.insert(COLUMN_VCF_POS,     "Pos")
    columns.insert(COLUMN_VCF_REF,     "Ref")
    columns.insert(COLUMN_VCF_ALT,     "Alt")
    return columns

