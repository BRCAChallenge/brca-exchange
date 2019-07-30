
import os
import sys
import pandas as pd
import hgvs.parser
import hgvs.dataproviders.uta
import hgvs.assemblymapper
import hgvs.normalizer
import hgvs.projector
import hgvs.validator
from hgvs.exceptions import HGVSError
import bioutils.assemblies
import logging
import pickle
import subprocess
import numpy as np
import click
# TODO: fix
sys.path.append(os.path.realpath(os.path.join(__file__, '../..')))

print(sys.path)

import types
from common.types import VCFVariant

# TODO: find better name
class HgvsProcessor:
    GRCh38_Assem = 'GRCh38'
    GRCh37_Assem = 'GRCh37'

    # TODO: refactor, put into proper utils
    def __init__(self):
        # TODO: should be parameter
        self.hgvs_dp = hgvs.dataproviders.uta.connect(
        'postgresql://anonymous@localhost:50827/uta/uta_20170629')
        self.hgvs_parser = hgvs.parser.Parser()
        self.hgvs_norm = hgvs.normalizer.Normalizer(self.hgvs_dp)

        assemblies = [self.GRCh37_Assem, self.GRCh38_Assem]

        self.hgvs_ams = { a : (hgvs.assemblymapper.AssemblyMapper(self.hgvs_dp, assembly_name=a, normalize=False, prevalidation_level=None)) for a in assemblies }

        self.hgvs_ams_normalizing = {a: (
            hgvs.assemblymapper.AssemblyMapper(self.hgvs_dp, assembly_name=a,
                                               normalize=True)) for a in
                         assemblies}

        all_assemblies = bioutils.assemblies.get_assemblies()
        # TODO: separate object?
        self.contig_maps = {a: ({
            s['name']: s['refseq_ac']
            for s in all_assemblies[a]['sequences'] if
            s['refseq_ac'] is not None
        }) for a in assemblies
        }

    def to_cdna(self, hgvs_obj, assembly = GRCh38_Assem):
        am = self.hgvs_ams[assembly]

        try:
            tr = am.relevant_transcripts(hgvs_obj)

            if tr:
                return am.g_to_c(hgvs_obj, tr[0])

        except HGVSError as e:
            logging.info("CDNA conversion issues " + str(e))

        return None


def _to_protein(hgvs_cdna, hgvs_proc):
    if not hgvs_cdna:
        return None

    try:
        #var_c1_norm = hgvs_proc.hgvs_norm.normalize(
        #    hgvs_cdna)  # TODO: for better error msgs?
        return str(hgvs_proc.hgvs_ams[HgvsProcessor.GRCh38_Assem].c_to_p(hgvs_cdna))
    except hgvs.exceptions.HGVSError as e:
        logging.info("Protein conversion issues with " + str(hgvs_cdna) + ": " + str(e))

    return str(None)


def determine_cdna():
    # attempt to convert, if not present try to substitue with data from other field

    #
    pass



def convert_to_hg37(vars):
    def pseudo_vcf_entry(v):
        entries = [v.chr, v.pos, '.', v.ref, v.alt, '', '', '']
        return '\t'.join([str(s) for s in entries])

    lst = [pseudo_vcf_entry(v) for v in vars]

    vcf_tmp = '/tmp/myvcf_all.vcf' # TODO better bat
    with open(vcf_tmp, 'w') as f:
        f.write('\n'.join(lst))

    vcf_tmp_out = '/tmp/myvcf_all_out.vcf'
    brca_resources_dir = '/Users/marc/brca/nobackup/enigma_wdir/resources'
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

    return [VCFVariant(v[0], v[1], v[3], v[4]) for v in
     [l.strip().split('\t') for l in vcf_out_lines]]



@click.command()
@click.argument('input', click.Path(readable=True))
@click.argument('output', click.Path(writable=True))
@click.option("--pkl") # TODO: keep?
def main(input, output, pkl):
    log_file_path = "/Users/marc/brca/nobackup/pseudonym_generator.log"
    logging.basicConfig(filename=log_file_path, filemode="w", level=logging.INFO,
                        format=' %(asctime)s %(filename)-15s %(message)s')

    # TODO: parametrize?
    os.environ[
        'HGVS_SEQREPO_DIR'] = '/Users/marc/brca/nobackup/enigma_wdir/resources/seq_repo/latest'

    df = pd.read_csv(input,
                     sep='\t')

    #df.loc[df['Synonyms'].isna(), 'Synonyms'] = '-'

    hgvs_proc = HgvsProcessor()

    def proc(x):
        v = VCFVariant(x['Chr'], x['Pos'], x['Ref'], x['Alt'])
        return str(v), hgvs_proc.to_cdna(v.to_hgvs_obj(hgvs_proc.contig_maps[HgvsProcessor.GRCh38_Assem]))


    d = pickle.load(open(pkl))
    #print(len(d))

    # calcualte CDNA only for missing variants?

    # TODO: how to log stuff?

    # do column wise, one function per attribute


    # conversion to protein

    df_sub = df #df.iloc[0:200]

    #d = {k: v for (k, v) in df_sub.apply(proc, axis=1)}
    #pickle.dump(d, open(pkl, 'w'))

    def get_cdna(x):
        v = VCFVariant(x['Chr'], x['Pos'], x['Ref'], x['Alt'])
        return d[str(v)]

    df_sub['tmp_hgvs_cdna_unorm'] = df_sub.apply(get_cdna, axis=1)

    def cdna_from_cdna_field(x):
        if x['HGVS_cDNA'] and x['HGVS_cDNA'] != '-':
            c = x['HGVS_cDNA']
            if ',' in c:
                c = c.split(',')[0]

            return hgvs_proc.hgvs_parser.parse(x["Reference_Sequence"] + ":" + c)
        return None

    df_sub['tmp_hgvs_cdna_unorm'] = df_sub.apply(lambda x: x['tmp_hgvs_cdna_unorm'] if x[
        'tmp_hgvs_cdna_unorm'] else None, axis=1)

    # TODO: paralellize?
    #df_sub['tmp_hgvs_cdna'] = df_sub.apply(lambda x: proc(x)[1], axis=1)

    def normalizing(v):
        if v:
            try:
                return hgvs_proc.hgvs_norm.normalize(v)
            except Exception as e:
                logging.info("Issues with normalizing " + str(v) + ": " + str(e))
        else:
            None

    df_sub['tmp_hgvs_cdna_norm'] = df_sub['tmp_hgvs_cdna_unorm'].apply(normalizing)

    df_sub['tmp_hgvs_g37'] = df_sub['tmp_hgvs_cdna_unorm'].apply(lambda hgvs_cdna: hgvs_proc.hgvs_ams[HgvsProcessor.GRCh37_Assem].c_to_g(hgvs_cdna) if hgvs_cdna else None)


    df_sub['tmp_hgvs_cdna_completed'] = df_sub.apply(lambda x:  x['tmp_hgvs_cdna_unorm'] if x['tmp_hgvs_cdna_unorm'] else cdna_from_cdna_field(x), axis=1)

    df_sub['pyhgvs_cDNA'] = df_sub['tmp_hgvs_cdna_completed'].apply(str)

    df_sub.loc[:, 'pyhgvs_Genomic_Coordinate_38'] = df_sub.apply(
        lambda x: str(VCFVariant(x['Chr'], x['Pos'], x['Ref'], x['Alt'])),
        axis=1)

    # conversion to 37
    #df_sub.loc[:, 'pyhgvs_Genomic_Coordinate_37'] = df_sub['tmp_hgvs_g37'].apply(lambda h: str(VCFVariant.from_hgvs_obj(h)) if h else None)

    cc = convert_to_hg37(df_sub.apply(lambda x: VCFVariant(x['Chr'], x['Pos'], x['Ref'], x['Alt']), axis=1))

    print(cc[1:10])
    print(len(cc))
    df_sub.loc[:, 'pyhgvs_Genomic_Coordinate_37'] = [str(v) for v in cc]

    #df_sub.loc[:, 'pyhgvs_Hg37_Start'] = df_sub['tmp_hgvs_g37'].apply(
    #    lambda h: h.posedit.pos.start.base if h else None)

    #df_sub.loc[:, 'pyhgvs_Hg37_End'] = df_sub['tmp_hgvs_g37'].apply(
    #    lambda h: h.posedit.pos.end.base if h else None)

    df_sub.loc[:, 'pyhgvs_Protein'] = (df_sub['tmp_hgvs_cdna_unorm'].
                                       #fillna(df_sub['tmp_hgvs_cdna']). # TODO: does this help?
                                       apply(lambda hgvs_cdna: _to_protein(hgvs_cdna, hgvs_proc)))


    # filter on U and NM

    # TODO: generalize
    brca1_trans = {'NM_007294.2', 'NM_007300.3', 'NM_007299.3',
                              'NM_007298.3', 'NM_007297.3', 'U14680.1'}
    brca2_trans = {'U43746.1'}

    def get_synonyms(x):
        synonyms = []

        for _, _, _, dst, alt_ac, method in hgvs_proc.hgvs_dp.get_tx_for_gene(x['Gene_Symbol']):
            #print(alt_ac +  " " + x['hgvs_cdna'].ac)

            if(dst in brca1_trans or dst in brca2_trans):

                for vc in [x['tmp_hgvs_cdna_completed']]:
                    if not vc:
                        continue
                    # TODO: need to optimize?
                    pj = hgvs.projector.Projector(hdp=hgvs_proc.hgvs_dp,
                                                  alt_ac=alt_ac,
                                                  src_ac=vc.ac,
                                                  dst_ac=dst, dst_alt_aln_method=method)

                    try:
                        vp = pj.project_variant_forward(vc)
                        synonyms.append(vp)
                        vp_norm = normalizing(vp)
                        if vp_norm:
                            if vp_norm not in synonyms:
                                logging.info("Found new synonym! " + str(vp_norm) + " " + str(vp) + " " + str(x['pyhgvs_Genomic_Coordinate_38']))
                            synonyms.append(vp_norm)
                    except Exception as e:
                        logging.info("Exception in synonym handling " + str(vc) + " with " + str(dst) + "using " + str(method) + " via " + str(dst) + " : " + str(e))

        return list({str(s) for s in synonyms})

    df_sub.loc[:, "new_syns"] = df_sub.apply(get_synonyms, axis=1)

    # TODO: sort within?

    df_sub['Synonyms'] = df_sub.apply(lambda x: ','.join(sorted(list(set( (x["Synonyms"].split(',') if (str(x["Synonyms"]) != 'nan' and x["Synonyms"] != '-') else list() )+ x["new_syns"])))), axis=1)


    df_sub.to_csv(output, sep='\t', index=False)


if __name__ == "__main__":
    main()