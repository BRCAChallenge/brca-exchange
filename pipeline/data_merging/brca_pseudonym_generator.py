import csv
import logging
import subprocess
import tempfile
from typing import Dict, List, Iterable, Optional
from os import path
import sys

import click
import hgvs.assemblymapper
import hgvs.dataproviders.uta
import hgvs.normalizer
import hgvs.projector
from hgvs.exceptions import HGVSError
from hgvs.normalizer import Normalizer
from hgvs.sequencevariant import SequenceVariant

from common import config
from common import utils
from common import ucsc
from common.hgvs_utils import HgvsWrapper
from common.variant_utils import VCFVariant

ALT_COL = 'Alt'
CHR_COL = 'Chr'
GENE_SYMBOL_COL = 'Gene_Symbol'
HG38_END_COL = 'Hg38_End'
HG38_START_COL = 'Hg38_Start'
HGVS_CDNA_COL = 'HGVS_cDNA'
POS_COL = 'Pos'
PYHGVS_CDNA_COL = 'pyhgvs_cDNA'
PYHGVS_GENOMIC_COORDINATE_37_COL = 'pyhgvs_Genomic_Coordinate_37'
PYHGVS_GENOMIC_COORDINATE_38_COL = 'pyhgvs_Genomic_Coordinate_38'
PYHGVS_HG37_END_COL = 'pyhgvs_Hg37_End'
PYHGVS_HG37_START_COL = 'pyhgvs_Hg37_Start'
PYHGVS_PROTEIN_COL = 'pyhgvs_Protein'
GENOMIC_HGVS_HG37_COL = 'Genomic_HGVS_37'
GENOMIC_HGVS_HG38_COL = 'Genomic_HGVS_38'
REFERENCE_SEQUENCE_COL = 'Reference_Sequence'
REF_COL = 'Ref'
SYNONYMS_COL = 'Synonyms'

# temporary fields
VAR_OBJ_FIELD = 'var_objs'
NEW_SYNONYMS_FIELD = 'new_syns'
TMP_HGVS_HG38 = 'tmp_Genomic_HGVS_38'
TMP_HGVS_HG37 = 'tmp_Genomic_HGVS_37'
TMP_HGVS_HG38_LEFT_ALIGNED = 'tmp_Genomic_HGVS_38_left'
TMP_HGVS_HG37_LEFT_ALIGNED = 'tmp_Genomic_HGVS_37_left'
TMP_CDNA_FROM_SOURCE = "tmp_hgvs_cdna_source"
TMP_CDNA_UNORM_FIELD = 'tmp_HGVS_CDNA_FIELD_unorm'
TMP_CDNA_NORM_FIELD = 'tmp_HGVS_CDNA_FIELD_norm'
TMP_CDNA_NORM_LEFT_ALINGED_FIELD = 'tmp_HGVS_CDNA_FIELD_left'
TMP_PROTEIN_LEFT_ALINGED_FIELD = 'tmp_Protein_Field_left'
TMP_PREFERRED_TRANSCRIPT_FIELD = 'tmp_Preferred_Transcript'


def _normalize_genomic_coordinates(hgvs_obj: Optional[SequenceVariant], strand: str, hgvs_norm_3: Normalizer, hgvs_norm_5: Normalizer):
    normalizer = hgvs_norm_3 if strand == config.POSITIVE_STRAND else hgvs_norm_5
    if hgvs_obj is None:
        return None
    try:
        return normalizer.normalize(hgvs_obj)
    except HGVSError as e:
        logging.warning("Issue normalizing genomic coordinates {}: {}".format(hgvs_obj, e))
    return None


def cdna_from_cdna_field(row: dict, cdna_ac_dict: Dict[str, str], hgvs_proc: HgvsWrapper):
    """Attempt to extract cDNA representation from already existing field (i.e. from the original source) if it looks reasonable"""
    if row[HGVS_CDNA_COL] and row[HGVS_CDNA_COL].startswith('c.') and row[REFERENCE_SEQUENCE_COL] == cdna_ac_dict[
        row[GENE_SYMBOL_COL]]:
        c = row[HGVS_CDNA_COL]
        if ',' in c:
            c = c.split(',')[0]

        return hgvs_proc.hgvs_parser.parse(row[REFERENCE_SEQUENCE_COL] + ":" + c)

    return None


def convert_to_hg37(vars: Iterable[VCFVariant], brca_resources_dir: str):
    """Lifting over hg38 variants to hg37 using crossmap

    Not using hgvs library since it doesn't handle intronic variants.
    """

    def pseudo_vcf_entry(v):
        entries = [v.chr, v.pos, '.', v.ref, v.alt, '', '', '']
        return '\t'.join([str(s) for s in entries])

    lst = [pseudo_vcf_entry(v) for v in vars]
    vcf_tmp = tempfile.mktemp('.vcf')
    with open(vcf_tmp, 'w') as f:
        f.write('\n'.join(lst))
    vcf_tmp_out = tempfile.mktemp('.vcf')
    args = ["CrossMap.py", "vcf",
            brca_resources_dir + "/hg38ToHg19.over.chain.gz",
            vcf_tmp,
            brca_resources_dir + "/hg19.fa",
            vcf_tmp_out]
    logging.info("Running CrossMap.py to convert to hg19")
    sp = subprocess.Popen(args)
    out, err = sp.communicate()
    if out:
        logging.info("standard output of subprocess: {}".format(out))
    if err:
        logging.info("standard output of subprocess: {}".format(err))
    vcf_out_lines = open(vcf_tmp_out, 'r').readlines()
    if path.exists(vcf_tmp_out + '.unmap'):
        vcf_out_failed_lines = open(vcf_tmp_out + '.unmap', 'r').readlines()
        return ([VCFVariant(v[0], int(v[1]), v[3], v[4]) for v in [l.strip().split('\t') for l in vcf_out_lines]],
                [VCFVariant(v[0], int(v[1]), v[3], v[4]) for v in [l.strip().split('\t') for l in vcf_out_failed_lines]])
    else:
        return ([VCFVariant(v[0], int(v[1]), v[3], v[4]) for v in [l.strip().split('\t') for l in vcf_out_lines]], [])


def handle_failed_hg37_translations(rows, var_objs_hg37, var_objs_hg37_failed, tmp_hgvs_hg37_values, hgvs_proc):
    # set None values for any hg37 coordinates that cannot be derived and add to lists in proper order
    hg38_to_idx = {row[GENOMIC_HGVS_HG38_COL]: i for i, row in enumerate(rows)}
    for v in var_objs_hg37_failed:
        hgvs_obj = v.to_hgvs_obj(hgvs_proc.contig_maps[HgvsWrapper.GRCh38_Assem])
        row_number = hg38_to_idx[str(hgvs_obj)]
        var_objs_hg37.insert(row_number, None)
        tmp_hgvs_hg37_values.insert(row_number, None)
        logging.info("Could not compute hg37 representation of internal for {}".format(str(hgvs_obj)))
    return (var_objs_hg37, tmp_hgvs_hg37_values)


def get_synonyms(row: dict, hgvs_proc: HgvsWrapper, syn_ac_dict: Dict[str, List[str]]):
    """ Determine other representations a variant may be known as """

    # take the representations from the source and all 'left' aligned (non-standard) representations
    synonyms = [str(row[TMP_CDNA_FROM_SOURCE]),
                str(row[TMP_HGVS_HG37_LEFT_ALIGNED]),
                str(row[TMP_HGVS_HG38_LEFT_ALIGNED]),
                str(row[TMP_CDNA_NORM_LEFT_ALINGED_FIELD]),
                str(row[TMP_PROTEIN_LEFT_ALINGED_FIELD])]

    if not row[GENE_SYMBOL_COL]:
        return []
    # calculate other representations wrt other accensions supported by the hgvs library
    for _, _, _, dst, alt_ac, method in hgvs_proc.hgvs_dp.get_tx_for_gene(row[GENE_SYMBOL_COL]):
        if row[GENE_SYMBOL_COL] not in syn_ac_dict:
            continue
        accessions = syn_ac_dict[row[GENE_SYMBOL_COL]]
        if dst in accessions:
            cdna_rep_list = [row[TMP_CDNA_NORM_FIELD]]
            if row[TMP_CDNA_NORM_FIELD] != row[TMP_CDNA_NORM_LEFT_ALINGED_FIELD]:
                cdna_rep_list = [row[TMP_CDNA_NORM_FIELD], row[TMP_CDNA_NORM_LEFT_ALINGED_FIELD]]
            for vc in cdna_rep_list:
                if not vc:
                    continue
                try:
                    pj = hgvs.projector.Projector(hdp=hgvs_proc.hgvs_dp,
                                                  alt_ac=alt_ac,
                                                  src_ac=vc.ac,
                                                  dst_ac=dst, dst_alt_aln_method=method)
                    vp = pj.project_variant_forward(vc)
                    synonyms.append(vp)
                except HGVSError as e:
                    logging.info("Exception in synonym handling " + str(vc) + " from " + str(vc.ac) + " to " +
                                 str(dst) + " using " + str(method) + " via " + str(alt_ac) + " : " + str(
                        e) + " " + str(e.__class__.__name__))
    return list({str(s) for s in synonyms}) # making sure, every representation appears only once


def _merge_and_clean_synonyms(row: dict):
    """ Merging with synonmys already determined from sources and cleaning up"""
    orig_list = [s for s in row[SYNONYMS_COL].split(',') if s]  # filter away ''
    combined = set(orig_list + row[NEW_SYNONYMS_FIELD])
    list_sorted = sorted(list(combined))
    # remove redundant representations, which we are already included in other columns
    reps_other_cols = {row[GENOMIC_HGVS_HG37_COL], row[GENOMIC_HGVS_HG38_COL], row[PYHGVS_CDNA_COL], row[PYHGVS_PROTEIN_COL],
                       row[HGVS_CDNA_COL]}
    list_sorted_cleaned = [v for v in list_sorted if v not in reps_other_cols and v != str(None) and v != '-']
    return ','.join(list_sorted_cleaned)


@click.command()
@click.argument('input', type=click.Path(readable=True))
@click.argument('output', type=click.Path(writable=True))
@click.option('--log-path', default='pseudonym_generator.log', help="Log file path")
@click.option("--config-file", required=True, help="path to gene configuration file")
@click.option('--resources', help="path to directory containing reference sequences")
def main(input, output, log_path, config_file, resources):
    csv.field_size_limit(sys.maxsize)
    utils.setup_logfile(log_path)

    cfg_df = config.load_config(config_file)
    syn_ac_dict = {r[config.SYMBOL_COL]: r[config.SYNONYM_AC_COL].split(';') for _, r in cfg_df.iterrows()}
    cdna_default_ac_dict = {r[config.SYMBOL_COL]: r[config.HGVS_CDNA_DEFAULT_AC] for _, r in cfg_df.iterrows()}
    strand_dict = {r[config.SYMBOL_COL]: r[config.STRAND_COL] for _, r in cfg_df.iterrows()}

    hgvs_proc = HgvsWrapper()

    # normalizer pairs for right-shift: positive-strand uses shuffle=3, negative-strand uses shuffle=5
    hgvs_norm_3 = hgvs.normalizer.Normalizer(hgvs_proc.hgvs_dp, shuffle_direction=3)
    hgvs_norm_5 = hgvs.normalizer.Normalizer(hgvs_proc.hgvs_dp, shuffle_direction=5)
    # normalizer pairs for left-align: directions reversed
    hgvs_norm_3_left = hgvs.normalizer.Normalizer(hgvs_proc.hgvs_dp, shuffle_direction=5)
    hgvs_norm_5_left = hgvs.normalizer.Normalizer(hgvs_proc.hgvs_dp, shuffle_direction=3)

    logging.info("Loading data from {}".format(input))
    with open(input, 'r', newline='') as f:
        reader = csv.DictReader(f, delimiter='\t')
        rows = list(reader)

    for row in rows:
        row[VAR_OBJ_FIELD] = VCFVariant(row[CHR_COL], int(row[POS_COL]), row[REF_COL], row[ALT_COL])

    logging.info("Converting variants to hgvs objects")
    for row in rows:
        row[TMP_HGVS_HG38] = row[VAR_OBJ_FIELD].to_hgvs_obj(hgvs_proc.contig_maps[HgvsWrapper.GRCh38_Assem])

    logging.info("Normalize genomic representation")
    for row in rows:
        row[TMP_HGVS_HG38] = _normalize_genomic_coordinates(
            row[TMP_HGVS_HG38], strand_dict.get(row[GENE_SYMBOL_COL]), hgvs_norm_3, hgvs_norm_5)
    for row in rows:
        row[TMP_HGVS_HG38_LEFT_ALIGNED] = _normalize_genomic_coordinates(
            row[TMP_HGVS_HG38], strand_dict.get(row[GENE_SYMBOL_COL]), hgvs_norm_3_left, hgvs_norm_5_left)

    for row in rows:
        row[GENOMIC_HGVS_HG38_COL] = str(row[TMP_HGVS_HG38])

    logging.info("Compute hg37 representation of internal representation")
    var_objs_hg37, var_objs_hg37_failed = convert_to_hg37([row[VAR_OBJ_FIELD] for row in rows], resources)
    tmp_hgvs_hg37_values = [v.to_hgvs_obj(hgvs_proc.contig_maps[HgvsWrapper.GRCh37_Assem]) for v in var_objs_hg37]
    var_objs_hg37, tmp_hgvs_hg37_values = handle_failed_hg37_translations(
        rows, var_objs_hg37, var_objs_hg37_failed, tmp_hgvs_hg37_values, hgvs_proc)

    for i, row in enumerate(rows):
        row[TMP_HGVS_HG37] = tmp_hgvs_hg37_values[i]

    logging.info("Compute hg37 normalized representation of internal")
    # normalizing again for the hg37 representation. An alternative would be to convert the normalized hg38 representation to hg37.
    # If we use crossmap, we would need a way to convert the VCF like representation back to an hgvs object, which we currently
    # are unable to do properly. That is, we can use VCFVariant.to_hgvs_obj, however, structural variants will be converted
    # to delins, losing information if a variant was e.g. a del, ins, or dup.
    for row in rows:
        row[GENOMIC_HGVS_HG37_COL] = _normalize_genomic_coordinates(
            row[TMP_HGVS_HG37], strand_dict.get(row[GENE_SYMBOL_COL]), hgvs_norm_3, hgvs_norm_5)
    for row in rows:
        row[GENOMIC_HGVS_HG37_COL] = str(row[GENOMIC_HGVS_HG37_COL])
    for row in rows:
        row[TMP_HGVS_HG37_LEFT_ALIGNED] = _normalize_genomic_coordinates(
            row[TMP_HGVS_HG37], strand_dict.get(row[GENE_SYMBOL_COL]), hgvs_norm_3_left, hgvs_norm_5_left)

    logging.info("Compute cDNA representation")
    for row in rows:
        row[TMP_PREFERRED_TRANSCRIPT_FIELD] = cdna_default_ac_dict.get(row[GENE_SYMBOL_COL], 'Not Found')

    for row in rows:
        row[TMP_CDNA_NORM_FIELD] = hgvs_proc.genomic_to_cdna(
            row[TMP_HGVS_HG38], target_transcript=row[TMP_PREFERRED_TRANSCRIPT_FIELD])
        row[TMP_CDNA_NORM_LEFT_ALINGED_FIELD] = hgvs_proc.genomic_to_cdna(
            row[TMP_HGVS_HG38_LEFT_ALIGNED], target_transcript=row[TMP_PREFERRED_TRANSCRIPT_FIELD])

    # extract cdna from source if it could not be computed
    for row in rows:
        row[TMP_CDNA_FROM_SOURCE] = row[HGVS_CDNA_COL]  # "backup" to be used later during synonym computation
        if not row[TMP_CDNA_NORM_FIELD]:
            row[TMP_CDNA_NORM_FIELD] = cdna_from_cdna_field(row, cdna_default_ac_dict, hgvs_proc)

    #### CDNA and Genomic HGVS conversions
    for row in rows:
        row[PYHGVS_CDNA_COL] = str(row[TMP_CDNA_NORM_FIELD])

    for row in rows:
        if row[PYHGVS_CDNA_COL].startswith("NM_"):
            parts = row[PYHGVS_CDNA_COL].split(':')
            row[REFERENCE_SEQUENCE_COL] = parts[0]
            row[HGVS_CDNA_COL] = parts[1]
        else:
            # still setting a reference sequence for downstream steps, even though no cDNA could be determined
            row[REFERENCE_SEQUENCE_COL] = cdna_default_ac_dict[row[GENE_SYMBOL_COL]]
            row[HGVS_CDNA_COL] = '-'

    #### Internal Genomic Coordinates
    for i, row in enumerate(rows):
        v = var_objs_hg37[i]
        row[PYHGVS_GENOMIC_COORDINATE_38_COL] = str(row[VAR_OBJ_FIELD])
        row[PYHGVS_GENOMIC_COORDINATE_37_COL] = str(v)
        # handles missing hg37 coordinates
        row[PYHGVS_HG37_START_COL] = v.pos if v else None
        row[PYHGVS_HG37_END_COL] = (v.pos + (int(row[HG38_END_COL]) - int(row[HG38_START_COL]))) if v else None

    #### Protein
    logging.info("Protein Conversion")
    for row in rows:
        row[PYHGVS_PROTEIN_COL] = str(hgvs_proc.cdna_to_protein(row[TMP_CDNA_NORM_FIELD]))
        row[TMP_PROTEIN_LEFT_ALINGED_FIELD] = str(hgvs_proc.cdna_to_protein(row[TMP_CDNA_NORM_LEFT_ALINGED_FIELD]))

    #### Synonyms
    logging.info("Compute Synonyms")
    for row in rows:
        row[NEW_SYNONYMS_FIELD] = get_synonyms(row, hgvs_proc, syn_ac_dict)

    for row in rows:
        row[SYNONYMS_COL] = (row[SYNONYMS_COL] or '').strip()

    # merge existing synonyms with generated ones and sort them
    for row in rows:
        row[SYNONYMS_COL] = _merge_and_clean_synonyms(row)

    #### Writing out
    # removing temporary fields
    tmp_fields = [k for k in (rows[0] if rows else {}) if k.startswith('tmp_')]
    fields_to_remove = set([VAR_OBJ_FIELD, NEW_SYNONYMS_FIELD] + tmp_fields)
    for row in rows:
        for field in fields_to_remove:
            row.pop(field, None)

    logging.info(f"Writing out to {output}")
    fieldnames = list(rows[0].keys()) if rows else []
    with open(output, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
