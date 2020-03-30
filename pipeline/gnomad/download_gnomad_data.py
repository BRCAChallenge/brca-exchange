#!/usr/bin/env python
import requests
import json
import time
import numpy as np
import pandas as pd
from pandas.io.json import json_normalize
import argparse
from math import floor, log10, isnan


def fetch(jsondata, url="https://gnomad.broadinstitute.org/api"):
    # The server gives a generic error message if the content type isn't
    # explicitly set
    headers = {"Content-Type": "application/json"}
    response = requests.post(url, json=jsondata, headers=headers)
    json = response.json()
    if "errors" in json:
        raise Exception(str(json["errors"]))
    return json

def unique_variant_set(variant_list_of_dicts):
    """ 
    Given a list of dicts, for which the items are variant IDs (as returned by the gnomAD API), 
    generate a set of unique set of variant IDs
    """
    variant_set = set()
    for item in variant_list_of_dicts:
        variant_set.add(list(item.values())[0])
    return variant_set

def transcript_to_variants(transcript_id, dataset):
    """
    Given a transcript, return the list of variants that map to the exons
    of the transcript, and were observed in samples from the indicated
    dataset.
    """
    fmt_graphql = """
    {
        transcript(transcript_id: "%s") {
          variants(dataset: %s) {
            variant_id: variantId
          }                                                                     
        }                                                                       
      }                                                                         
    """
    req_variantlist = {
        "query": fmt_graphql % (transcript_id, dataset),
        "variables": {}
    }
    response = fetch(req_variantlist)
    return unique_variant_set(response["data"]["transcript"]["variants"])

def gene_to_coords(gene_name):
    """                                                                         
    Given a gene symbol, return the coordinates.                                
    """
    graphql_query = """                                                         
    {                                                                           
        gene(gene_name: "%s") {                         
            chrom                                                               
            start                                                               
            stop                                                                
        }                                                                       
    }"""
    graphql_request = {
        "query": graphql_query % (gene_name),
        "variables": {}
    }
    response = fetch(graphql_request)
    return response["data"]["gene"]

def coords_to_variants(chrom, start, stop, dataset_id):
    region_query = """
        query GnomadRegion($chrom: String!, $start: Int!, $stop: Int!, 
                           $dataset_id: DatasetId!) {
            region(chrom: $chrom, start: $start, stop: $stop) {
                variants(dataset: $dataset_id) {
                   variantId
                }
            }
         }
         """
    region_variables = {
        "chrom": chrom,
        "start": start,
        "stop": stop,
        "dataset_id": dataset_id
        }
    headers = { "content-type": "application/json" }
    response = requests.post(
        'http://gnomad.broadinstitute.org/api',
        json={ "query": region_query, "variables": region_variables },
        headers=headers)
    parse = json.loads(response.text)
    variants = parse['data']['region']['variants']
    return(unique_variant_set(variants))

def gene_to_region_variants(gene_name, dataset_id):
    """
    Given a gene name, return the list of variants via a region
    query.  These will mostly be the intronic variants.
    """
    coords = gene_to_coords(gene_name)
    variantList = coords_to_variants(coords["chrom"], coords["start"],
                                     coords["stop"], dataset_id)
    return(set(variantList))

def fetch_data_for_one_variant(variant_id, dataset, max_retries=5):
    variant_detail_query = """
        query GnomadVariant($variantId: String!, $datasetId: DatasetId!) {
            variant(variantId: $variantId, dataset: $datasetId) {
                alt
                chrom
                pos
                ref
                variantId                                              
                ... on GnomadVariantDetails {
                        flags
                        sortedTranscriptConsequences {
                            transcript_id
                            hgvsc
                        }
                    }
                    exome {
                        ac
                        an
                        ac_hom
                        faf95 {
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
                       filters
                            populations {
                            id
                            ac
                            an
                            ac_hom
                       }                 
                    }
                }
             }"""
    headers = { "content-type": "application/json" }
    print("Fetching", variant_id)
    variant_detail_variables = {
        "variantId": variant_id,
        "datasetId": dataset,
    }
    retries = 0
    while retries < max_retries:
        try:
            response = requests.post(
                'http://gnomad.broadinstitute.org/api',
                json={
                    "query": variant_detail_query,
                    "variables": variant_detail_variables
                },
                headers=headers)
            parse = json.loads(response.text)
        except ValueError:
            retries += 1
            time.sleep(0.1)
        else:
            return(parse['data']['variant'])
    return None


def variant_set_to_variant_data(variants, dataset):
    variant_details = []
    for this_variant in variants:
        variant_data = fetch_data_for_one_variant(this_variant, dataset)
        if variant_data is not None:
            variant_details.append(variant_data)
            time.sleep(0.01)
    return(variant_details)


def compile_allele_values(df):
    populations = ['AFR', 'AFR_FEMALE', 'AFR_MALE', 'AMR', 'AMR_FEMALE', 'AMR_MALE', 'ASJ', 'ASJ_FEMALE', 'ASJ_MALE', 'EAS',
                   'EAS_JPN', 'EAS_KOR', 'EAS_OEA', 'EAS_FEMALE', 'EAS_MALE', 'FIN', 'FIN_FEMALE', 'FIN_MALE', 'NFE', 'NFE_BGR',
                   'NFE_EST', 'NFE_NWE', 'NFE_ONF', 'NFE_SEU', 'NFE_SWE', 'NFE_FEMALE', 'NFE_MALE', 'OTH', 'OTH_FEMALE', 'OTH_MALE',
                   'SAS', 'SAS_FEMALE', 'SAS_MALE', 'FEMALE', 'MALE']
    for population in populations:
        df['genome_' + population + '_af'] = calculate_frequency(df['genome_' + population + '_ac'], df['genome_' + population + '_an'])
        df['exome_' + population + '_af'] = calculate_frequency(df['exome_' + population + '_ac'], df['exome_' + population + '_an'])
    df['exome_af'] = calculate_frequency(df['exome_ac'], df['exome_an'])
    df['genome_af'] = calculate_frequency(df['genome_ac'], df['genome_an'])
    return df


def calculate_frequency(ac, an):
    freq = pd.to_numeric(ac, errors='coerce').divide(pd.to_numeric(an, errors='coerce'))
    return freq.apply(round_four_sigfigs)


def round_four_sigfigs(num):
    if isnan(num):
        return num
    elif num == 0 or num == 0.0:
        return 0
    else:
        return round(num, -int(floor(log10(abs(num))) - (3)))


def flatten(variant, field, genome_or_exome):
    for f in field:
        if f not in ['populations', 'filters']:
            variant[genome_or_exome + '_' + f] = field[f]
    populations = field['populations']
    for population in populations:
        name = population['id']
        keys = population.keys()
        for key in keys:
            if name != key:
                variant[genome_or_exome + '_' + name + '_' + key] = population[key]
    return variant


def flatten_populations(variants):
    for variant in variants:
        genome = variant['genome']
        exome = variant['exome']
        if genome:
            variant = flatten(variant, genome, 'genome')
        if exome:
            variant = flatten(variant, exome, 'exome')
        del variant['genome']
        del variant['exome']
    return variants


def find_correct_hgvs(variants, transcripts):
    """
    Given the set of transcript IDs that we queried for (one per gene), and
    given the data for one particular variant, return the cDNA HGVS string
    for that variant, corresponding to one of the transcripts we're intereted
    in.  There should be one and only one.  If there is no such transcript
    and HGVS string for this variant, then make a note and throw out the variant.
    """
    variants_in_expected_transcripts = []
    for variant in variants:
        sortedTranscriptConsequences = variant["sortedTranscriptConsequences"]
        for transcript_hgvs in sortedTranscriptConsequences:
            if transcript_hgvs["transcript_id"] in transcripts:
                variant["hgvs"] = transcript_hgvs["hgvsc"].split(':')[1]
                variant["transcript"] = transcript_hgvs["transcript_id"]
                del variant["sortedTranscriptConsequences"]
                variants_in_expected_transcripts.append(variant)
                break
        if "hgvs" not in variant:
            print("Warning: variant %s falls outside the expected transcripts"
                  % variant["variantId"])
    return variants_in_expected_transcripts


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-o', '--output', help='Ouput tsv file result.')
    options = parser.parse_args()
    return options


def main():
    f_out = parse_args().output
    dataset = "gnomad_r2_1_non_cancer"
    brca1_transcript ="ENST00000357654"
    brca2_transcript = "ENST00000544455"
    transcripts = (brca1_transcript, brca2_transcript)

    # organize brca1 request
    brca1_exonic_variants = transcript_to_variants(brca1_transcript, dataset)
    brca1_intronic_variants = gene_to_region_variants("BRCA1", dataset)
    brca1_variants = brca1_intronic_variants | brca1_exonic_variants

    # organize brca2 request
    brca2_exonic_variants = transcript_to_variants(brca2_transcript, dataset)
    brca2_intronic_variants = gene_to_region_variants("BRCA2", dataset)
    brca2_variants = brca2_intronic_variants | brca2_exonic_variants

    # combine requests and get brca1 and brca2 data from gnomAD
    brca12_variants = brca1_variants | brca2_variants
    brca12_variant_data = variant_set_to_variant_data(brca12_variants, dataset)

    # find hgvs, flatten, convert to dataframe, compute allele frequencies, and normalize
    variants_with_hgvs = find_correct_hgvs(brca12_variant_data, transcripts)
    variants_with_flattened_populations = flatten_populations(variants_with_hgvs)
    variants_df = json_normalize(variants_with_flattened_populations)
    variants_df['flags'] = variants_df['flags'].apply(', '.join)
    df_with_allele_values = compile_allele_values(variants_df)
    stringified_df_with_allele_values = df_with_allele_values.replace(np.nan, '-', regex=True).replace('', '-', regex=True)

    # output to .tsv
    stringified_df_with_allele_values.to_csv(f_out, sep='\t', index=False)


if __name__ == "__main__":
    main()