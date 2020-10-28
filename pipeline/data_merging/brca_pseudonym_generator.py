import logging
import subprocess
import tempfile
from typing import Dict, List, Iterable

import click
import hgvs.assemblymapper
import hgvs.dataproviders.uta
import hgvs.projector
import pandas as pd
from hgvs.exceptions import HGVSError
from hgvs.normalizer import Normalizer
from hgvs.sequencevariant import SequenceVariant

from common import config
from common import utils
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


def _normalize_genomic_coordinates(hgvs_obj: SequenceVariant, strand: str, hgvs_norm_3: Normalizer, hgvs_norm_5: Normalizer):
    normalizer = hgvs_norm_3 if strand == config.POSITIVE_STRAND else hgvs_norm_5

    try:
        return normalizer.normalize(hgvs_obj)
    except HGVSError as e:
        logging.warning("Issue normalizing genomic coordinates {}: {}".format(hgvs_obj, e))
    return None


def _normalize_genomic_fnc(src_col: str, target_col: str, right_shift: bool, strand_dict: Dict[str, str]):
    """Returns a function to normalize genomic coordinates

    Used for parallelizing normalization
    """
    def _ret(df_part):
        hgvs_proc = HgvsWrapper()  # create it again in subprocess to avoid pickle issues when otherwise copying it to subprocess
        hgvs_norm_3 = hgvs.normalizer.Normalizer(hgvs_proc.hgvs_dp, shuffle_direction=3 if right_shift else 5)
        hgvs_norm_5 = hgvs.normalizer.Normalizer(hgvs_proc.hgvs_dp, shuffle_direction=5 if right_shift else 3)

        df_part[target_col] = df_part.apply(lambda r: _normalize_genomic_coordinates(r[src_col],
                                                                                     strand_dict.get(
                                                                                         r[GENE_SYMBOL_COL]),
                                                                                     hgvs_norm_3,
                                                                                     hgvs_norm_5), axis=1)
        return df_part

    return _ret


def cdna_from_cdna_field(row: pd.Series, cdna_ac_dict: Dict[str, str], hgvs_proc: HgvsWrapper):
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

    return [VCFVariant(v[0], int(v[1]), v[3], v[4]) for v in
            [l.strip().split('\t') for l in vcf_out_lines]]


def get_synonyms(row: pd.Series, hgvs_proc: HgvsWrapper, syn_ac_dict: Dict[str, List[str]]):
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


def _merge_and_clean_synonyms(row: pd.Series):
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
@click.option('--log-path', default='pseudonym_generator.log', help="Log file pth")
@click.option("--config-file", required=True, help="path to gene configuration file")
@click.option('--resources', help="path to directory containing reference sequences")
@click.option('--processes', type=int, help='Number of processes to use for parallelization', default=8)
def main(input, output, log_path, config_file, resources, processes):
    utils.setup_logfile(log_path)

    cfg_df = config.load_config(config_file)

    syn_ac_dict = {r[config.SYMBOL_COL]: r[config.SYNONYM_AC_COL].split(';') for _, r in cfg_df.iterrows()}
    cdna_default_ac_dict = {r[config.SYMBOL_COL]: r[config.HGVS_CDNA_DEFAULT_AC] for _, r in cfg_df.iterrows()}
    strand_dict = { r[config.SYMBOL_COL] : r[config.STRAND_COL] for _, r in cfg_df.iterrows() }

    hgvs_proc = HgvsWrapper()

    logging.info("Loading data from {}".format(input))
    df = pd.read_csv(input, sep='\t')

    df[VAR_OBJ_FIELD] = df.apply(lambda r: VCFVariant(r[CHR_COL], r[POS_COL], r[REF_COL], r[ALT_COL]), axis=1)

    logging.info("Converting variants to hgvs objects")

    df[TMP_HGVS_HG38] = df[VAR_OBJ_FIELD].apply(
        lambda v: v.to_hgvs_obj(hgvs_proc.contig_maps[HgvsWrapper.GRCh38_Assem]))

    logging.info("Normalize genomic representation")

    df = utils.parallelize_dataframe(df, _normalize_genomic_fnc(TMP_HGVS_HG38,
                                                                TMP_HGVS_HG38, True,
                                                                strand_dict), processes)
    df = utils.parallelize_dataframe(df, _normalize_genomic_fnc(TMP_HGVS_HG38,
                                                                TMP_HGVS_HG38_LEFT_ALIGNED, False,
                                                                strand_dict), processes)

    df[GENOMIC_HGVS_HG38_COL] = df[TMP_HGVS_HG38].apply(str)

    logging.info("Compute hg37 representation of internal representation")
    var_objs_hg37 = convert_to_hg37(df[VAR_OBJ_FIELD], resources)

    logging.info("Compute hg37 normalized representation of internal")
    # normalizing again for the hg37 representation. An alternative would be to convert the normalized hg38 representation to hg37.
    # If we use crossmap, we would need a way to convert the VCF like representation back to an hgvs object, which we currently
    # are unable to do properly. That is, we can use VCFVariant.to_hgvs_obj, however, structural variants will be converted
    # to delins, losing information if a variant was e.g. a del, ins, or dup.

    df[TMP_HGVS_HG37] = pd.Series(
        [v.to_hgvs_obj(hgvs_proc.contig_maps[HgvsWrapper.GRCh37_Assem]) for v in var_objs_hg37])

    df = utils.parallelize_dataframe(df, _normalize_genomic_fnc(TMP_HGVS_HG37,
                                                                GENOMIC_HGVS_HG37_COL, True, strand_dict), processes)
    df[GENOMIC_HGVS_HG37_COL] = df[GENOMIC_HGVS_HG37_COL].apply(str)

    df = utils.parallelize_dataframe(df, _normalize_genomic_fnc(TMP_HGVS_HG37,
                                                                TMP_HGVS_HG37_LEFT_ALIGNED, False, strand_dict), processes)

    logging.info("Compute cDNA representation")

    def _compute_cdna(df_part):
        hgvs_proc = HgvsWrapper()  # create it to avoid pickle issues when copying to subprocess
        df_part[TMP_CDNA_NORM_FIELD] = df_part[TMP_HGVS_HG38].apply(
            lambda hgvs_obj: hgvs_proc.genomic_to_cdna(hgvs_obj))
        df_part[TMP_CDNA_NORM_LEFT_ALINGED_FIELD] = df_part[TMP_HGVS_HG38_LEFT_ALIGNED].apply(
            lambda hgvs_obj: hgvs_proc.genomic_to_cdna(hgvs_obj))
        return df_part

    df = utils.parallelize_dataframe(df, _compute_cdna, processes)

    # extract cdna from source if it could not be computed
    df[TMP_CDNA_FROM_SOURCE] = df[HGVS_CDNA_COL]  # "backup" to be used later during synonym computation
    df[TMP_CDNA_NORM_FIELD] = df.apply(
        lambda r: cdna_from_cdna_field(r, cdna_default_ac_dict, hgvs_proc) if not r[TMP_CDNA_NORM_FIELD] else r[
            TMP_CDNA_NORM_FIELD], axis=1)

    #### CDNA and Genomic HGVS conversions
    df[PYHGVS_CDNA_COL] = df[TMP_CDNA_NORM_FIELD].apply(str)

    available_cdna = df[PYHGVS_CDNA_COL].str.startswith("NM_")
    df.loc[available_cdna, REFERENCE_SEQUENCE_COL] = df.loc[available_cdna, PYHGVS_CDNA_COL].str.split(':').apply(
        lambda l: l[0])
    df.loc[available_cdna, HGVS_CDNA_COL] = df.loc[available_cdna, PYHGVS_CDNA_COL].str.split(':').apply(lambda l: l[1])

    # still setting a reference sequence for downstream steps, even though no cDNA could be determined
    df.loc[~available_cdna, REFERENCE_SEQUENCE_COL] = df.loc[~available_cdna, GENE_SYMBOL_COL].apply(
        lambda g: cdna_default_ac_dict[g])
    df.loc[~available_cdna, HGVS_CDNA_COL] = '-'

    #### Internal Genomic Coordinates
    df[PYHGVS_GENOMIC_COORDINATE_38_COL] = df[VAR_OBJ_FIELD].apply(lambda v: str(v))

    df[PYHGVS_GENOMIC_COORDINATE_37_COL] = pd.Series([str(v) for v in var_objs_hg37])
    df[PYHGVS_HG37_START_COL] = pd.Series([v.pos for v in var_objs_hg37])
    df[PYHGVS_HG37_END_COL] = df[PYHGVS_HG37_START_COL] + (df[HG38_END_COL] - df[HG38_START_COL])

    #### Protein
    logging.info("Protein Conversion")

    def _compute_proteins(df_part):
        hgvs_proc = HgvsWrapper()
        df_part[PYHGVS_PROTEIN_COL] = df_part[TMP_CDNA_NORM_FIELD].apply(lambda x: str(hgvs_proc.cdna_to_protein(x)))
        df_part[TMP_PROTEIN_LEFT_ALINGED_FIELD] = df_part[TMP_CDNA_NORM_LEFT_ALINGED_FIELD].apply(
            lambda x: str(hgvs_proc.cdna_to_protein(x)))
        return df_part

    df = utils.parallelize_dataframe(df, _compute_proteins, processes)

    #### Synonyms
    logging.info("Compute Synonyms")

    def _compute_synonyms(df_part):
        hgvs_proc = HgvsWrapper()
        df_part[NEW_SYNONYMS_FIELD] = df_part.apply(lambda s: get_synonyms(s, hgvs_proc, syn_ac_dict), axis=1)
        return df_part

    df = utils.parallelize_dataframe(df, _compute_synonyms, processes)

    df[SYNONYMS_COL] = df[SYNONYMS_COL].fillna('').str.strip()

    # merge existing synonyms with generated ones and sort them
    df[SYNONYMS_COL] = df.apply(_merge_and_clean_synonyms, axis=1)

    #### Writing out
    # removing temporary fields
    tmp_fields = [c for c in df.columns if c.startswith('tmp_')]
    df = df.drop(columns=[VAR_OBJ_FIELD, NEW_SYNONYMS_FIELD] + tmp_fields)

    logging.info(f"Writing out to {output}")
    df.to_csv(output, sep='\t', index=False)


if __name__ == "__main__":
    main()
