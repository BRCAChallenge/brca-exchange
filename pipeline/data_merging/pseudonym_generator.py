import argparse
import csv
import logging
import re
import subprocess
import sys
import tempfile
import time
from os import path
from typing import Dict, Iterable, List, Optional

import requests
import hgvs.dataproviders.uta
from hgvs.assemblymapper import AssemblyMapper
from hgvs.exceptions import HGVSError
from hgvs.normalizer import Normalizer
from hgvs.parser import Parser
from hgvs.sequencevariant import SequenceVariant
from hgvs.variantmapper import VariantMapper

from common import config, ucsc, utils
from common.hgvs_utils import HgvsWrapper
from common.variant_utils import VCFVariant

CA_ID_COL = 'CA_ID'
CLINGEN_ALLELE_REGISTRY_ENDPOINT = "http://reg.clinicalgenome.org/alleles?file=hgvs"
MAX_TRIES = 5
CA_ID_RE = re.compile(r'CA[0-9]+')

RESP_ID_KEY = '@id'
RESP_ERROR_TYPE_KEY = 'errorType'
RESP_DESCRIPTION_KEY = 'description'
RESP_MESSAGE_KEY = 'message'

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


def _extract_ca_id(r, orig_hgvs):
    if RESP_ID_KEY not in r:
        if RESP_ERROR_TYPE_KEY in r:
            logging.warning(
                f"Server could not process {orig_hgvs}: "
                f"{r[RESP_ERROR_TYPE_KEY]}. {r[RESP_DESCRIPTION_KEY]}. {r[RESP_MESSAGE_KEY]}")
        return None

    ca_id = r[RESP_ID_KEY].split('/')[-1]
    if CA_ID_RE.match(ca_id):
        return ca_id

    logging.warning(f"Could not extract ca_id for {orig_hgvs} out of {ca_id}")
    return None


def _perform_request(url, data):
    tries = 0
    req = None
    while not req and tries < MAX_TRIES:
        try:
            req = requests.post(url, data)
        except requests.exceptions.RequestException as e:
            logging.error(f"Request to {url} failed. Going to retry")
            logging.exception(e)
            time.sleep(10)
            tries += 1

    if not req:
        logging.error(f"Requests failed {MAX_TRIES} times, exiting.")
        sys.exit(1)

    if req.status_code != 200:
        logging.error(
            f"Request failed, server responded with status code "
            f"{req.status_code} and message {req.json()}")
        sys.exit(1)

    return req.json()


def _extract_grch37_allele(r):
    for allele in r.get('genomicAlleles', []):
        if allele.get('referenceGenome') == 'GRCh37':
            hgvs_list = allele.get('hgvs', [])
            return hgvs_list[0] if hgvs_list else None
    return None


def _query_clingen_allele_registry(hgvs_ids: Iterable[str], url: str, max_chunk_size: int = 50000) -> list:
    """Submit HGVS IDs to the ClinGen Allele Registry and return the raw JSON response list."""
    resp_all = []
    for split in utils.split_list_in_chunks(hgvs_ids, max_chunk_size):
        hgvs_data = '\n'.join(split)
        logging.info(f"Requesting chunk of size {len(split)} from {url}")
        resp_all.extend(_perform_request(url, hgvs_data))
    return resp_all


def _extract_mane_select_cdna(r: dict) -> Optional[str]:
    """Return the MANE Select RefSeq nucleotide HGVS string, or None."""
    for transcript_allele in (r.get('transcriptAlleles') or []):
        mane = transcript_allele.get('MANE') if transcript_allele else None
        if mane and mane.get('maneStatus') == 'MANE Select':
            nucleotide = mane.get('nucleotide')
            refseq_nt = nucleotide.get('RefSeq') if nucleotide else None
            return refseq_nt.get('hgvs') if refseq_nt else None
    return None


def _extract_mane_select_protein(r: dict) -> Optional[str]:
    """Return the MANE Select RefSeq protein HGVS string, or None."""
    for transcript_allele in (r.get('transcriptAlleles') or []):
        mane = transcript_allele.get('MANE') if transcript_allele else None
        if mane and mane.get('maneStatus') == 'MANE Select':
            protein = mane.get('protein')
            refseq_prot = protein.get('RefSeq') if protein else None
            return refseq_prot.get('hgvs') if refseq_prot else None
    return None


def _extract_transcript_allele_hgvs(r: dict) -> List[str]:
    """Return all HGVS strings from transcriptAlleles for a single registry response object.

    Collects every entry in transcriptAlleles[x]['hgvs'] (a list) and the single
    value in transcriptAlleles[x]['proteinEffect']['hgvs'], skipping any that are absent or None.
    """
    results = []
    for transcript_allele in (r.get('transcriptAlleles') or []):
        if not transcript_allele:
            continue
        for hgvs_str in (transcript_allele.get('hgvs') or []):
            if hgvs_str:
                results.append(hgvs_str)
        protein_effect = transcript_allele.get('proteinEffect')
        if protein_effect:
            protein_hgvs = protein_effect.get('hgvs')
            if protein_hgvs:
                results.append(protein_hgvs)
    community_title = r.get('communityStandardTitle')
    if community_title:
        results.extend(community_title)
    return results


def _extract_fields_from_registry_response(resp_all: list, hgvs_ids: Iterable[str]) -> tuple:
    """Extract CA IDs and GRCh37 HGVS strings from a ClinGen Allele Registry response list."""
    ca_ids = [_extract_ca_id(r, h) for r, h in zip(resp_all, hgvs_ids)]
    grch37_hgvs = [_extract_grch37_allele(r) for r in resp_all]
    mane_select_cDNA = [_extract_mane_select_cdna(r) for r in resp_all]
    mane_select_protein = [_extract_mane_select_protein(r) for r in resp_all]
    transcript_allele_hgvs = [_extract_transcript_allele_hgvs(r) for r in resp_all]
    return ca_ids, grch37_hgvs, mane_select_cDNA, mane_select_protein, transcript_allele_hgvs


def _query_server(hgvs_ids: Iterable[str], url: str, max_chunk_size: int = 50000):
    resp_all = _query_clingen_allele_registry(hgvs_ids, url, max_chunk_size)
    return _extract_fields_from_registry_response(resp_all, hgvs_ids)


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
TMP_CDNA_NORM_LEFT_ALIGNED_FIELD = 'tmp_HGVS_CDNA_FIELD_left'
TMP_PROTEIN_LEFT_ALIGNED_FIELD = 'tmp_Protein_Field_left'
TMP_PREFERRED_TRANSCRIPT_FIELD = 'tmp_Preferred_Transcript'


def genomic_hgvs_to_coords(hgvs_string: str, parser: Parser = None) -> tuple:
    """Parse a genomic HGVS string and return (start, end) genomic coordinates.

    Args:
        hgvs_string: A genomic HGVS string, e.g. 'NC_000017.10:g.41246747del'
        parser: An hgvs Parser instance. If None, a new one is created.

    Returns:
        A (start, end) tuple of integer coordinates, or (None, None) on failure.
    """
    if parser is None:
        parser = Parser()
    try:
        var = parser.parse(hgvs_string)
        start = int(str(var.posedit.pos.start))
        end = int(str(var.posedit.pos.end))
        return start, end
    except HGVSError as e:
        logging.warning(f"Could not parse genomic HGVS string '{hgvs_string}': {e}")
        return None, None


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


def get_synonyms(gene_symbol, genomic_hgvs, target_transcripts, hdp, normalizer, assembly_mapper):
    """ Determine other representations a variant may be known as """

    # calculate other representations wrt other accensions supported by the hgvs library
    synonyms = []
    for _, _, _, transcript, alt_ac, method in hdp.get_tx_for_gene(gene_symbol):
        if transcript in target_transcripts:
            try:
                var_c = assembly_mapper.g_to_c(genomic_hgvs, transcript)
                synonyms.append(str(normalizer.normalize(var_c)))
            except HGVSError as e:
                logging.info("Exception in synonym handling for mapping %s to transcript %s" % (str(genomic_hgvs), transcript))
    return list({str(s) for s in synonyms}) # making sure, every representation appears only once


def _filter_and_extend_synonyms(synonyms: List[str], extra: List[str]) -> List[str]:
    """Remove NM_/NP_ accession entries from synonyms, then append extra synonyms."""
    filtered = [s for s in synonyms if not (s.startswith('NM_') or s.startswith('NP_'))]
    return filtered + extra


def _merge_and_clean_synonyms(row: dict, new_synonyms):
    """ Merging with synonmys already determined from sources and cleaning up"""
    orig_list = [s for s in row[SYNONYMS_COL].split(',') if s]  # filter away ''
    combined = set(orig_list + new_synonyms)
    list_sorted = sorted(list(combined))
    # remove redundant representations, which we are already included in other columns
    reps_other_cols = {row[PYHGVS_GENOMIC_COORDINATE_37_COL], row[PYHGVS_GENOMIC_COORDINATE_38_COL],
                       row[PYHGVS_CDNA_COL], row[PYHGVS_PROTEIN_COL],
                       row[HGVS_CDNA_COL]}
    list_sorted_cleaned = [v for v in list_sorted if v not in reps_other_cols and v != str(None) and v != '-']
    return ','.join(list_sorted_cleaned)


def main_old():
    args = parse_args()
    input = args.input
    output = args.output
    log_path = args.log_path
    config_file = args.config_file
    resources = args.resources

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
        row[TMP_CDNA_NORM_LEFT_ALIGNED_FIELD] = hgvs_proc.genomic_to_cdna(
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
        row[TMP_PROTEIN_LEFT_ALIGNED_FIELD] = str(hgvs_proc.cdna_to_protein(row[TMP_CDNA_NORM_LEFT_ALIGNED_FIELD]))

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

def parse_args():
    parser = argparse.ArgumentParser(description="Generate pseudonyms for BRCA variants")
    parser.add_argument('input', help="Input TSV file path")
    parser.add_argument('output', help="Output TSV file path")
    parser.add_argument('--logfile', default='pseudonym_generator.log', help="Log file path")
    parser.add_argument('--configfile', required=True, help="path to gene configuration file")
    parser.add_argument('--resources', help="path to directory containing reference sequences")
    parser.add_argument('--debug', help="Turn on extra debugging info?",
                        action='store_true', default=False)
    return parser.parse_args()


def genomic_hgvs(chr, pos, ref, alt, normalizer, hgvs_proc):
    this_variant = VCFVariant(chr, pos, ref, alt)
    tmp_genomic_hgvs_38 = this_variant.to_hgvs_obj(hgvs_proc.contig_maps[HgvsWrapper.GRCh38_Assem])
    normalizer.normalize(tmp_genomic_hgvs_38)
    return str(tmp_genomic_hgvs_38)


def map_via_seqrepo(this_gene, genomic_hgvs_str, default_cdna, normalizer,
                    hgvs_proc, assembly_mapper_38, assembly_mapper_37, debug=True):
    try:
        tmp_cdna_hgvs = normalizer.normalize(am38.g_to_c(tmp_genomic_hgvs_38,
                                                         default_cdna))
    except (hgvs.exceptions.HGVSInvalidIntervalError,
            hgvs.exceptions.HGVSUnsupportedOperationError):
        if debug:
            print("Interval error.  Setting all the pyhgvs columns to None")
        pyhgvs_cdna = None
        genomic_coordinate_37 = None
        protein = None
    else:
        if debug:
            print("No interval error!  cDNA", str(tmp_cdna_hgvs))
        pyhgvs_cdna = str(tmp_cdna_hgvs)
        target_chrom_37 = hgvs_proc.contig_maps['GRCh37'][str(chrom_ac_dict[this_gene])]
        tmp_genomic_hgvs_37 = normalizer.normalize(assembly_mapper_37.c_to_g(tmp_cdna_hgvs))
        genomic_coordinate_37 = str(tmp_genomic_hgvs_37)
        if debug:
            print("genomic HGVS (37)", str(tmp_genomic_hgvs_37))
        protein = str(assembly_mapper_38.c_to_p(tmp_cdna_hgvs))
        if debug:
            print("Protein HGVS:", protein)
    return pyhgvs_cdna, genomic_coordinate_37, protein


def main():
    args = parse_args()
    csv.field_size_limit(sys.maxsize)
    utils.setup_logfile(args.logfile)

    config_df = config.load_config(args.configfile)
    cdna_default_ac_dict = {r[config.SYMBOL_COL]: r[config.HGVS_CDNA_DEFAULT_AC] for _, r in config_df.iterrows()}
    syn_ac_dict = {r[config.SYMBOL_COL]: r[config.SYNONYM_AC_COL].split(';') for _, r in config_df.iterrows()}
    chrom_ac_dict = {r[config.SYMBOL_COL]: r[config.CHROM_COL] for _, r in config_df.iterrows()}

    hdp = hgvs.dataproviders.uta.connect()
    hp = Parser()
    genomic_normalizer = Normalizer(hdp)
    am38 = AssemblyMapper(hdp, assembly_name="GRCh38", alt_aln_method="splign", normalize=True)
    am37 = AssemblyMapper(hdp, assembly_name="GRCh37", alt_aln_method="splign", normalize=True)
    hgvs_proc = HgvsWrapper()

    with open(args.input, 'r') as input_fp:
        reader = csv.DictReader(input_fp, delimiter='\t')
        new_fieldnames = [
            PYHGVS_GENOMIC_COORDINATE_38_COL, PYHGVS_CDNA_COL, PYHGVS_GENOMIC_COORDINATE_37_COL,
            PYHGVS_HG37_START_COL, PYHGVS_HG37_END_COL, PYHGVS_PROTEIN_COL, CA_ID_COL,
        ]
        rows = list(reader)
        output_fieldnames = reader.fieldnames + new_fieldnames

        hgvs_values = [
            genomic_hgvs(row[CHR_COL], int(row[POS_COL]), row[REF_COL], row[ALT_COL], genomic_normalizer, hgvs_proc)
            for row in rows
        ]

        logging.info("Querying ClinGen Allele Registry for CA IDs")
        ca_ids, grch37_hgvs, mane_select_cdna, mane_select_protein, transcript_allele_hgvs = \
            _query_server([str(v) for v in hgvs_values], CLINGEN_ALLELE_REGISTRY_ENDPOINT)

        for row, ca_id, grch37_hgvs_val, cdna, protein, synonyms in zip(
                rows, ca_ids, grch37_hgvs, mane_select_cdna, mane_select_protein, transcript_allele_hgvs):
            row[CA_ID_COL] = ca_id
            row[PYHGVS_GENOMIC_COORDINATE_37_COL] = grch37_hgvs_val
            row[PYHGVS_CDNA_COL] = cdna
            row[PYHGVS_PROTEIN_COL] = protein
            row[SYNONYMS_COL] = _filter_and_extend_synonyms(row[SYNONYMS_COL], synonyms)

        processed_rows = []
        for row, genomic_hgvs_38 in zip(rows, hgvs_values):
            row[PYHGVS_GENOMIC_COORDINATE_38_COL] = genomic_hgvs_38
            if args.debug:
                print("working on variant", row[PYHGVS_GENOMIC_COORDINATE_38_COL])
            if row[PYHGVS_GENOMIC_COORDINATE_37_COL] is None:
                (row[PYHGVS_CDNA_COL], row[PYHGVS_GENOMIC_COORDINATE_37_COL], row[PYHGVS_PROTEIN_COL]) = \
                    map_via_seqrepo(row['Gene_Symbol'], genomic_hgvs_38,
                                    cdna_default_ac_dict[row['Gene_Symbol']],
                                    genomic_normalizer, hgvs_proc,
                                    am38, am37, debug=args.debug)
                if args.debug:
                    print("Updating synonyms")
                new_synonyms = get_synonyms(
                    row['Gene_Symbol'], genomic_hgvs_38, syn_ac_dict[row['Gene_Symbol']], hdp, gn, am38)
                row[SYNONYMS_COL] = _merge_and_clean_synonyms(row, new_synonyms)

            (row[PYHGVS_HG37_START_COL], row[PYHGVS_HG37_END_COL]) = genomic_hgvs_to_coords(row[PYHGVS_GENOMIC_COORDINATE_37_COL], hp)
            processed_rows.append(row)

        with open(args.output, mode='w') as output_fp:
            writer = csv.DictWriter(output_fp, fieldnames=output_fieldnames, delimiter='\t')
            writer.writeheader()
            writer.writerows(processed_rows)

if __name__ == "__main__":
    main()
