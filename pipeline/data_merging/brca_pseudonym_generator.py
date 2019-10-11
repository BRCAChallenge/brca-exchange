import logging
import os
import pickle
import subprocess
import tempfile

import click
import hgvs.assemblymapper
import hgvs.dataproviders.uta
import hgvs.normalizer
import hgvs.parser
import hgvs.projector
import hgvs.validator

import pandas as pd
from hgvs.exceptions import HGVSError

from common import config
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
REFERENCE_SEQUENCE_COL = 'Reference_Sequence'
REF_COL = 'Ref'
SYNONYMS_COL = 'Synonyms'

# temporary fields
VAR_OBJ_FIELD = 'var_objs'
NEW_SYNONYMS_FIELD = 'new_syns'
TMP_CDNA_UNORM_FIELD = 'tmp_HGVS_CDNA_FIELD_unorm'
TMP_CDNA_NORM_FIELD = 'tmp_HGVS_CDNA_FIELD_norm'


def _get_cdna(df, pkl, hgvs_proc, cdna_ac_dict, normalize):
    def cdna_from_cdna_field(x):
        if x[HGVS_CDNA_COL] and x[HGVS_CDNA_COL].startswith('c.') and x[REFERENCE_SEQUENCE_COL] == cdna_ac_dict[x[GENE_SYMBOL_COL]]:
            c = x[HGVS_CDNA_COL]
            if ',' in c:
                c = c.split(',')[0]

            return hgvs_proc.hgvs_parser.parse(x[REFERENCE_SEQUENCE_COL] + ":" + c)

        return None

    def compute_hgvs(x):
        v = VCFVariant(x[CHR_COL], x[POS_COL], x[REF_COL], x[ALT_COL])
        v = v.to_hgvs_obj(hgvs_proc.contig_maps[HgvsWrapper.GRCh38_Assem])

        if normalize:
            vn = hgvs_proc.normalizing(v)
            v = vn if vn else v

        return hgvs_proc.to_cdna(v)

    def from_field_or_compute(row):
        computed = compute_hgvs(row)

        if not computed:
            return cdna_from_cdna_field(row)
        return computed

    var_objs = df[VAR_OBJ_FIELD]

    if pkl and os.path.exists(pkl):
        pickle_dict = pickle.load(open(pkl, 'rb'))
        s_cdna = pd.Series([ pickle_dict[str(v)]  for v in var_objs  ])
    else:
        s_cdna = df.apply(from_field_or_compute, axis=1)

        if pkl:
            pickle_dict = {str(v): c for (c, v) in zip(s_cdna, var_objs)}
            pickle.dump(pickle_dict, open(pkl, 'wb'))

    return s_cdna


def convert_to_hg37(vars, brca_resources_dir):
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


def get_synonyms(x, hgvs_proc, syn_ac_dict):
    synonyms = []

    for _, _, _, dst, alt_ac, method in hgvs_proc.hgvs_dp.get_tx_for_gene(x[GENE_SYMBOL_COL]):
        if x[GENE_SYMBOL_COL] not in syn_ac_dict:
            continue

        accessions = syn_ac_dict[x[GENE_SYMBOL_COL]]

        if dst in accessions:
            for vc in [x[TMP_CDNA_NORM_FIELD]]:
                if not vc:
                    continue

                try:
                    pj = hgvs.projector.Projector(hdp=hgvs_proc.hgvs_dp,
                                                  alt_ac=alt_ac,
                                                  src_ac=vc.ac,
                                                  dst_ac=dst, dst_alt_aln_method=method)


                    vp = pj.project_variant_forward(vc)
                    synonyms.append(vp)
                    vp_norm = hgvs_proc.normalizing(vp)
                    if vp_norm:
                        if vp_norm not in synonyms:
                            logging.info("Found new synonym! " + str(vp_norm) + " " + str(vp) + " " + str(x[PYHGVS_GENOMIC_COORDINATE_38_COL]))
                        synonyms.append(vp_norm)
                except HGVSError as e:
                    logging.info("Exception in synonym handling " + str(vc) + " from " + str(vc.ac) + " to " + str(dst) + " using " + str(method) + " via " + str(alt_ac) + " : " + str(e) + " " + str(e.__class__))

    return list({str(s) for s in synonyms})


def _merge_synonyms(x):
    orig_list = [s for s in x[SYNONYMS_COL].split(',') if s] # filter away ''

    combined = set(orig_list + x[NEW_SYNONYMS_FIELD])
    list_sorted = sorted(list(combined))
    return ','.join(list_sorted)


@click.command()
@click.argument('input', click.Path(readable=True))
@click.argument('output', click.Path(writable=True))
@click.option('--log-path', default='pseudonym_generator.log', help="Log file pth")
@click.option("--pkl", help="Saving HGVS cDNA objects to save time during development")
@click.option("--config-file", required=True, help="path to gene configuration file")
@click.option('--resources', help="path to directory containing reference sequences")
def main(input, output, pkl, log_path, config_file, resources):
    logging.basicConfig(filename=log_path, filemode="w", level=logging.INFO,
                        format=' %(asctime)s %(filename)-15s %(message)s')

    cfg_df = config.load_config(config_file)

    syn_ac_dict = { x[config.SYMBOL_COL] : x[config.SYNONYM_AC_COL].split(';') for _, x in cfg_df.iterrows()}
    cdna_default_ac_dict = { x[config.SYMBOL_COL] : x[config.HGVS_CDNA_DEFAULT_AC] for _, x in cfg_df.iterrows()}

    hgvs_proc = HgvsWrapper()

    df = pd.read_csv(input, sep='\t')

    df[VAR_OBJ_FIELD] = df.apply(lambda x: VCFVariant(x[CHR_COL], x[POS_COL], x[REF_COL], x[ALT_COL]), axis=1)

    #### CDNA conversions
    df[TMP_CDNA_NORM_FIELD] = _get_cdna(df, pkl, hgvs_proc, cdna_default_ac_dict, normalize=True)
    df[PYHGVS_CDNA_COL] = df[TMP_CDNA_NORM_FIELD].apply(str)

    available_cdna = df[PYHGVS_CDNA_COL].str.startswith("NM_")
    df.loc[available_cdna, REFERENCE_SEQUENCE_COL] = df.loc[available_cdna, PYHGVS_CDNA_COL].str.split(':').apply(lambda l: l[0])
    df.loc[available_cdna, HGVS_CDNA_COL] = df.loc[available_cdna, PYHGVS_CDNA_COL].str.split(':').apply(lambda l: l[1])

    # still setting a reference sequence for downstream steps, even though no cDNA could be determined
    df.loc[~available_cdna, REFERENCE_SEQUENCE_COL] = df.loc[~available_cdna, GENE_SYMBOL_COL].apply(lambda g: cdna_default_ac_dict[g])
    df.loc[~available_cdna, HGVS_CDNA_COL] = '-'

    #### Genomic Coordinates
    df[PYHGVS_GENOMIC_COORDINATE_38_COL] = df[VAR_OBJ_FIELD].apply(lambda v: str(v))

    var_objs_hg37 = convert_to_hg37(df[VAR_OBJ_FIELD], resources)
    df[PYHGVS_GENOMIC_COORDINATE_37_COL] = pd.Series([str(v) for v in var_objs_hg37])

    df[PYHGVS_HG37_START_COL] = pd.Series([v.pos for v in var_objs_hg37])
    df[PYHGVS_HG37_END_COL] = df[PYHGVS_HG37_START_COL] + (df[HG38_END_COL] - df[HG38_START_COL])

    #### Protein
    df[PYHGVS_PROTEIN_COL] = (df[TMP_CDNA_NORM_FIELD].
                              apply(lambda HGVS_CDNA_FIELD: str(hgvs_proc.to_protein(HGVS_CDNA_FIELD))))

    #### Synonyms
    df[NEW_SYNONYMS_FIELD] = df.apply(lambda s: get_synonyms(s, hgvs_proc, syn_ac_dict), axis=1)

    df[SYNONYMS_COL] = df[SYNONYMS_COL].fillna('').str.strip()

    # merge existing synonyms with generated ones and sort them
    df[SYNONYMS_COL] = df.apply(_merge_synonyms, axis=1)

    #### Writing out
    # cleaning up temporary fields
    df = df.drop(columns=[VAR_OBJ_FIELD, NEW_SYNONYMS_FIELD, TMP_CDNA_NORM_FIELD])

    df.to_csv(output, sep='\t', index=False)


if __name__ == "__main__":
    main()
