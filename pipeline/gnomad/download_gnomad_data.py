import json
import requests
import sys
import time
import argparse
import logging


def parse_args():
        parser = argparse.ArgumentParser(description='Download gnomad data and convert to .tsv format.')
        parser.add_argument('-o', '--output', type=argparse.FileType('w'),
                            help='Ouput json file result.')
        parser.add_argument('-l', '--logfile', default='/tmp/download_gnomad_data.log')
        parser.add_argument('-v', '--verbose', action='count', default=False, help='determines logging')
        options = parser.parse_args()
        return options

def gene_to_coordinates(gene_name):
    """ Given a gene name, return the chromosome and start/stop coordinates"""
    gene_query = """
    query GnomadGene($geneName: String!) {
      gene(gene_name: $geneName) {
        chrom
        start
        stop
      }
    }
    """
    gene_variables = {
        "geneName": gene_name,
    }
    headers = { "content-type": "application/json" }
    response = requests.post(
        'http://gnomad.broadinstitute.org/api',
        json={ "query": gene_query, "variables": gene_variables },
        headers=headers)
    parse = json.loads(response.text)
    gene_data = parse['data']['gene']
    chrom = gene_data['chrom']
    start = gene_data['start']
    stop = gene_data['stop']
    return(chrom, start, stop)

def region_to_variants(chrom, start, stop, subset):
    """Given a genomic region and a subset (e.g. non-cancer), return the variants in that
    region detected in that subset"""
    region_query = """
    query GnomadRegion($chrom: String!, $start: Int!, $stop: Int!, $datasetId: DatasetId!) {
      region(chrom: $chrom, start: $start, stop: $stop) {
        start
        stop
        chrom
        variants(dataset: $datasetId) {
          variantId
        }
      }
    }
    """
    region_variables = {
        "chrom": chrom,
        "start": start,
        "stop": stop,
        "datasetId": subset,
    }
    headers = { "content-type": "application/json" }
    response = requests.post(
        'http://gnomad.broadinstitute.org/api',
        json={ "query": region_query, "variables": region_variables },
        headers=headers)
    parse = json.loads(response.text)
    region_data = parse['data']['region']
    variant_ids = map(lambda v: v['variantId'], region_data['variants'])
    return(variant_ids)

def query_one_variant(one_variant_id, transcript, subset):
    """Get the data for one single variant"""
    variant_detail_query = """
    query GnomadVariant($variantId: String!, $datasetId: DatasetId!) {
      variant(variantId: $variantId, dataset: $datasetId) {
        alt
        chrom
        pos
        ref
        ... on GnomadVariantDetails {
          multiNucleotideVariants {
            combined_variant_id
            changes_amino_acids
            n_individuals
            other_constituent_snvs
          }
          exome {
            ac
            an
            ac_hom
            faf95 {
              popmax
              popmax_population
            }
            faf99 {
              popmax
              popmax_population
            }
            filters
            populations {
              id
              ac
              an
              ac_hom
            }
          }
          genome {
            ac
            an
            ac_hom
            faf95 {
              popmax
              popmax_population
            }
            faf99 {
              popmax
              popmax_population
            }
            filters
            populations {
              id
              ac
              an
              ac_hom
            }
          }
          flags
          sortedTranscriptConsequences {
            hgvs
            hgvsc
            hgvsp
            transcript_id
            lof
            lof_flags
            lof_filter
            lof_info
          }
        }
      }
    }
    """
    headers = { "content-type": "application/json" }
    variant_detail_variables = {
          "variantId": one_variant_id,
          "datasetId": subset,
        }
    max_retries = 3
    retry_count = 0
    success = False
    while retry_count < max_retries and not success:
        try:
            response = requests.post(
                'http://gnomad.broadinstitute.org/api',
                json={
                    "query": variant_detail_query,
                    "variables": variant_detail_variables
                },
                headers=headers)
            one_variant_data = json.loads(response.text)
        except json.decoder.JSONDecodeError:
            retry_count += 1
            one_variant_data = None
        else:
            success = True
            one_variant_data = select_transcript(one_variant_data, transcript)
    if one_variant_data == None:
        sys.stderr.write("Could not read data on variant %s\n" % one_variant_id)
    return(one_variant_data)


def select_transcript(variant_data, transcript):
    """Set a new field indicating the data from the selected transcript"""
    consequences_selected_transcript = None
    for item in variant_data["data"]["variant"]["sortedTranscriptConsequences"]:
        if item["transcript_id"] == transcript:
            consequences_selected_transcript = item
    variant_data["data"]["variant"]["consequencesForSelectedTranscript"] = consequences_selected_transcript
    # if variant data can't be found for selected transcript, skip variant
    if consequences_selected_transcript == None:
        return None
    return variant_data


def query_variant_details(variant_ids, transcript, subset, verbose=False):
    """Get the details on the indicated variants"""
    variant_details = []
    for variant_id in variant_ids:
        one_variant_detail = query_one_variant(variant_id, transcript, subset)
        if one_variant_detail is not None:
            one_variant_detail['data']['variant']['variant_id'] = variant_id
            variant_details.append(one_variant_detail['data']['variant'])
        if verbose:
            print('fetched', variant_id)
        time.sleep(0.01)
    return(variant_details)


def get_gnomad_data(gene_name, transcript, subset):
    """Get the variants for the indicated gene, transcript and subset"""
    (chrom, start, stop) = gene_to_coordinates(gene_name)
    variant_ids = region_to_variants(chrom, start, stop, subset)
    variant_details = query_variant_details(variant_ids, transcript, subset)
    return variant_details

if __name__ == "__main__":
    options = parse_args()
    output = options.output
    logfile = options.logfile

    if options.verbose:
        logging_level = logging.DEBUG
    else:
        logging_level = logging.CRITICAL

    logging.basicConfig(filename=logfile, filemode="w", level=logging_level)

    brca1v = get_gnomad_data("BRCA1", "ENST00000357654", "gnomad_r2_1_non_cancer")
    brca2v = get_gnomad_data("BRCA2", "ENST00000544455", "gnomad_r2_1_non_cancer")
    variants = brca1v + brca2v
    json.dump(variants, output)
