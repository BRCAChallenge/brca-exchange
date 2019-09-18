import logging
import os
import pickle
import subprocess
import tempfile
from multiprocessing.pool import ThreadPool

import click
import hgvs.assemblymapper
import hgvs.dataproviders.uta
import hgvs.normalizer
import hgvs.parser
import hgvs.projector
import hgvs.validator
import numpy as np
import pandas as pd
from hgvs.exceptions import HGVSError

from common import config
from common.hgvs_utils import HgvsWrapper
from common.variant_utils import VCFVariant

SYNONYMS_FIELD = 'Synonyms'

# temporary fields
VAR_OBJ_FIELD = 'var_objs'
NEW_SYNONYMS_FIELD = 'new_syns'
TMP_CDNA_UNORM_FIELD = 'tmp_hgvs_cdna_unorm'
TMP_CDNA_NORM_FIELD = 'tmp_hgvs_cdna_norm'

# TODO: use! (check pathos)
def parallelize_dataframe(df, func, n_cores=4):
    df_split = np.array_split(df, n_cores)

    pool = ThreadPool(n_cores)
    df = pd.concat(pool.map(func, df_split))
    pool.close()
    pool.join()
    return df


def _get_cdna(df, pkl, hgvs_proc):
    def cdna_from_cdna_field(x):
        if x['HGVS_cDNA'] and x['HGVS_cDNA'].startswith('c.'):
            c = x['HGVS_cDNA']
            if ',' in c:
                c = c.split(',')[0]

            return hgvs_proc.hgvs_parser.parse(x["Reference_Sequence"] + ":" + c)

        return None

    def compute_hgvs(x):
        v = VCFVariant(x['Chr'], x['Pos'], x['Ref'], x['Alt'])
        return hgvs_proc.to_cdna(v.to_hgvs_obj(hgvs_proc.contig_maps[HgvsWrapper.GRCh38_Assem]))

    def from_field_or_compute(row):
        r = cdna_from_cdna_field(row)
        if not r:
            r = compute_hgvs(row)
        return r

    var_objs = df[VAR_OBJ_FIELD]

    if pkl and os.path.exists(pkl):
        pickle_dict = pickle.load(open(pkl, 'rb'))
        s_cdna = pd.Series([ pickle_dict[str(v)]  for v in var_objs  ])
    else:
        s_cdna = parallelize_dataframe(df, lambda dfx: dfx.apply(from_field_or_compute, axis=1), 4)

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

    sp = subprocess.Popen(args) # TODO: make cleaner
    out, err = sp.communicate()
    if out:
        print("standard output of subprocess:")
        print out
    if err:
        print("standard error of subprocess:")
        print err

    vcf_out_lines = open(vcf_tmp_out, 'r').readlines()

    return [VCFVariant(v[0], int(v[1]), v[3], v[4]) for v in
     [l.strip().split('\t') for l in vcf_out_lines]]


def get_synonyms(x, hgvs_proc, syn_ac_dict):
    synonyms = []

    for _, _, _, dst, alt_ac, method in hgvs_proc.hgvs_dp.get_tx_for_gene(x['Gene_Symbol']):
        if x['Gene_Symbol'] not in syn_ac_dict:
            continue

        accessions = syn_ac_dict[x['Gene_Symbol']]

        if(dst in accessions):
            for vc in [x[TMP_CDNA_UNORM_FIELD]]:
                if not vc:
                    continue
                # TODO: need to optimize?

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
                            logging.info("Found new synonym! " + str(vp_norm) + " " + str(vp) + " " + str(x['pyhgvs_Genomic_Coordinate_38']))
                        synonyms.append(vp_norm)
                except HGVSError as e:
                    logging.info("Exception in synonym handling " + str(vc) + " from " + str(vc.ac) + " to " + str(dst) + " using " + str(method) + " via " + str(alt_ac) + " : " + str(e) + " " + str(e.__class__))

    return list({str(s) for s in synonyms})


def _merge_synonyms(x):
    orig_list = [s for s in x[SYNONYMS_FIELD].split(',') if s] # filter away ''

    combined = set(orig_list + x[NEW_SYNONYMS_FIELD])
    list_sorted = sorted(list(combined))
    return ','.join(list_sorted)


@click.command()
@click.argument('input', click.Path(readable=True))
@click.argument('output', click.Path(writable=True))
@click.option('--log-path', default='pseudonym_generator.log')
@click.option("--pkl") # TODO: keep?
@click.option("--config-file", required=True)
@click.option('--resources')
def main(input, output, pkl, log_path, config_file, resources):
    logging.basicConfig(filename=log_path, filemode="w", level=logging.INFO,
                        format=' %(asctime)s %(filename)-15s %(message)s')

    cfg_df = config.load_config(config_file)

    syn_ac_dict = { x[config.SYMBOL_COL] : x[config.SYNONYM_AC_COL].split(';') for _, x in cfg_df.iterrows()}
    cdna_default_ac_dict = { x[config.SYMBOL_COL] : x[config.HGVS_CDNA_DEFAULT_AC] for _, x in cfg_df.iterrows()}

    hgvs_proc = HgvsWrapper()

    df = pd.read_csv(input, sep='\t')

    df = df.iloc[0:1000] # TODO: remove

    df[VAR_OBJ_FIELD] = df.apply(lambda x: VCFVariant(x['Chr'], x['Pos'], x['Ref'], x['Alt']), axis=1)

    # CDNA conversions
    df[TMP_CDNA_UNORM_FIELD] = _get_cdna(df, pkl, hgvs_proc)
    df[TMP_CDNA_NORM_FIELD] = df[TMP_CDNA_UNORM_FIELD].apply(hgvs_proc.normalizing)
    df['pyhgvs_cDNA'] = df[TMP_CDNA_UNORM_FIELD].apply(str)

    available_cdna = df['pyhgvs_cDNA'].str.startswith("NM_")
    df.loc[available_cdna, 'Reference_Sequence'] = df.loc[available_cdna, 'pyhgvs_cDNA'].str.split(':').apply(lambda l: l[0])
    df.loc[available_cdna, 'HGVS_cDNA'] = df.loc[available_cdna, 'pyhgvs_cDNA'].str.split(':').apply(lambda l: l[1])

    # still setting a reference sequence for downstream steps, even though no cDNA could be determined
    df.loc[~available_cdna, 'Reference_Sequence'] = df.loc[~available_cdna, 'Gene_Symbol'].apply(lambda g: cdna_default_ac_dict[g])
    df.loc[~available_cdna, 'HGVS_cDNA'] = '-'

    # Genomic Coordinates
    df['pyhgvs_Genomic_Coordinate_38'] = df[VAR_OBJ_FIELD].apply(lambda v: str(v))

    var_objs_hg37 = convert_to_hg37(df[VAR_OBJ_FIELD], resources)
    df['pyhgvs_Genomic_Coordinate_37'] = pd.Series([str(v) for v in var_objs_hg37])

    df['pyhgvs_Hg37_Start'] = pd.Series([v.pos for v in var_objs_hg37])
    df['pyhgvs_Hg37_End'] = df['pyhgvs_Hg37_Start'] + (df['Hg38_End'] - df['Hg38_Start'])

    # ## Protein
    df['pyhgvs_Protein'] = (df[TMP_CDNA_NORM_FIELD].
                            apply(lambda hgvs_cdna: hgvs_proc._to_protein(hgvs_cdna)))

    df[NEW_SYNONYMS_FIELD] = df.apply(lambda s: get_synonyms(s, hgvs_proc, syn_ac_dict), axis=1)

    df[SYNONYMS_FIELD] = df[SYNONYMS_FIELD].fillna('').str.strip()

    # TODO: rename to generated_syonynms?
    # merge existing synonyms with generated ones and sort them
    df[SYNONYMS_FIELD] = df.apply(_merge_synonyms, axis=1)

    # cleaning up temporary fields
    df = df.drop(columns=[VAR_OBJ_FIELD, NEW_SYNONYMS_FIELD, TMP_CDNA_UNORM_FIELD, TMP_CDNA_NORM_FIELD])

    df.to_csv(output, sep='\t', index=False)


if __name__ == "__main__":
    main()
