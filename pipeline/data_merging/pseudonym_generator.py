import argparse
import csv
import logging
import re
import sys
import time
from typing import Iterable, List, Optional

import requests
import hgvs.dataproviders.uta
from hgvs.assemblymapper import AssemblyMapper
from hgvs.exceptions import HGVSError
from hgvs.normalizer import Normalizer
from hgvs.parser import Parser

from common import config, utils
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


def get_synonyms(gene_symbol, genomic_hgvs_string, target_transcripts, hdp, parser, normalizer, assembly_mapper):
    """ Determine other representations a variant may be known as """

    # calculate other representations wrt other accensions supported by the hgvs library
    genomic_hgvs = parser.parse(genomic_hgvs_string)
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
    orig_list = [s for s in row[SYNONYMS_COL] if s]  # filter away ''
    combined = set(orig_list + new_synonyms)
    list_sorted = sorted(list(combined))
    # remove redundant representations, which we are already included in other columns
    reps_other_cols = {row[PYHGVS_GENOMIC_COORDINATE_37_COL], row[PYHGVS_GENOMIC_COORDINATE_38_COL],
                       row[PYHGVS_CDNA_COL], row[PYHGVS_PROTEIN_COL],
                       row[HGVS_CDNA_COL]}
    list_sorted_cleaned = [v for v in list_sorted if v not in reps_other_cols and v != str(None) and v != '-']
    return ','.join(list_sorted_cleaned)


def parse_args():
    parser = argparse.ArgumentParser(description="Generate pseudonyms for BRCA variants")
    parser.add_argument('input', help="Input TSV file path")
    parser.add_argument('output', help="Output TSV file path")
    parser.add_argument('--logpath', default='pseudonym_generator.log', help="Log file path")
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


def map_via_seqrepo(this_gene, genomic_hgvs_38, default_cdna, normalizer,
                    hgvs_proc, assembly_mapper_38, assembly_mapper_37,
                    chrom, debug=True):
    genomic_hgvs_38_obj = hgvs_proc.hgvs_parser.parse(genomic_hgvs_38)
    try:
        tmp_cdna_hgvs = normalizer.normalize(assembly_mapper_38.g_to_c(genomic_hgvs_38_obj,
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
        target_chrom_37 = hgvs_proc.contig_maps['GRCh37'][chrom]
        tmp_genomic_hgvs_37 = normalizer.normalize(assembly_mapper_37.c_to_g(tmp_cdna_hgvs))
        genomic_coordinate_37 = str(tmp_genomic_hgvs_37)
        if debug:
            print("genomic HGVS (37)", str(tmp_genomic_hgvs_37))
        protein = str(assembly_mapper_38.c_to_p(tmp_cdna_hgvs))
        if debug:
            print("Protein HGVS:", protein)
    return pyhgvs_cdna, genomic_coordinate_37, protein


def ensure_mane_transcript_cdna(row, mane_transcript, hgvs_proc, normalizer, am38, debug=False):
    """If the cDNA in row does not use the MANE transcript, remap using SeqRepo."""
    if not row[PYHGVS_CDNA_COL] or not row[PYHGVS_CDNA_COL].startswith(mane_transcript) or row[REFERENCE_SEQUENCE_COL] is not mane_transcript:
        genomic_hgvs_38_obj = hgvs_proc.hgvs_parser.parse(str(row[PYHGVS_GENOMIC_COORDINATE_38_COL]))
        try:
            cdna_hgvs = normalizer.normalize(am38.g_to_c(genomic_hgvs_38_obj, mane_transcript))
            row[PYHGVS_CDNA_COL] = str(cdna_hgvs)
            parts = row[PYHGVS_CDNA_COL].split(':')
            row[REFERENCE_SEQUENCE_COL] = parts[0]
            row[HGVS_CDNA_COL] = parts[1]
        except (hgvs.exceptions.HGVSInvalidIntervalError,
                hgvs.exceptions.HGVSUnsupportedOperationError):
            logging.warning(
                f"Could not remap {row[PYHGVS_GENOMIC_COORDINATE_38_COL]} to MANE transcript {mane_transcript}")
            row[REFERENCE_SEQUENCE_COL] = mane_transcript
            row[HGVS_CDNA_COL] = "-"
            row[PYHGVS_CDNA_COL] = "-"


def main():
    args = parse_args()
    csv.field_size_limit(sys.maxsize)
    utils.setup_logfile(args.logpath)

    config_df = config.load_config(args.configfile)
    mane_transcript_dict = {r[config.SYMBOL_COL]: r[config.HGVS_CDNA_DEFAULT_AC] for _, r in config_df.iterrows()}
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
            GENOMIC_HGVS_HG38_COL, GENOMIC_HGVS_HG37_COL,
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
                print("working on variant", row[PYHGVS_GENOMIC_COORDINATE_38_COL], "reference",
                      row[REFERENCE_SEQUENCE_COL], "cDNA", row[HGVS_CDNA_COL],
                      "allele registry", row[CA_ID_COL])
            this_gene = row[GENE_SYMBOL_COL]
            if row[PYHGVS_GENOMIC_COORDINATE_37_COL] is None:
                (row[PYHGVS_CDNA_COL], row[PYHGVS_GENOMIC_COORDINATE_37_COL], row[PYHGVS_PROTEIN_COL]) = \
                    map_via_seqrepo(row['Gene_Symbol'], genomic_hgvs_38,
                                    mane_transcript_dict[this_gene],
                                    genomic_normalizer, hgvs_proc,
                                    am38, am37, str(chrom_ac_dict[this_gene]),
                                    debug=args.debug)
                if args.debug:
                    print("Updating synonyms via SeqRepo")
                new_synonyms = get_synonyms(
                    this_gene, genomic_hgvs_38, syn_ac_dict[this_gene], hdp,
                    hp, genomic_normalizer, am38)
                row[SYNONYMS_COL] = _merge_and_clean_synonyms(row, new_synonyms)
            if row[PYHGVS_GENOMIC_COORDINATE_37_COL]:
                (row[PYHGVS_HG37_START_COL], row[PYHGVS_HG37_END_COL]) = genomic_hgvs_to_coords(row[PYHGVS_GENOMIC_COORDINATE_37_COL], hp)
            else:
                row[PYHGVS_HG37_START_COL] = None
                row[PYHGVS_HG37_END_COL] = None
            ensure_mane_transcript_cdna(row, mane_transcript_dict[this_gene], hgvs_proc, genomic_normalizer,
                                        am38, debug=args.debug)
            row[GENOMIC_HGVS_HG38_COL] = row[PYHGVS_GENOMIC_COORDINATE_38_COL]
            row[GENOMIC_HGVS_HG37_COL] = row[PYHGVS_GENOMIC_COORDINATE_37_COL]
            processed_rows.append(row)

        with open(args.output, mode='w') as output_fp:
            writer = csv.DictWriter(output_fp, fieldnames=output_fieldnames, delimiter='\t')
            writer.writeheader()
            writer.writerows(processed_rows)

if __name__ == "__main__":
    main()
