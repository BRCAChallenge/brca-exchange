#!/usr/bin/env python
import requests
import json
import time
import pandas as pd
from math import floor, log10, isnan
import pdb

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

def variant_set_to_variant_data(variants, dataset):
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
    variant_details = []
    count = 0
    for this_variant in variants:
        variant_detail_variables = {
            "variantId": this_variant,
            "datasetId": dataset,
            }
        response = requests.post(
            'http://gnomad.broadinstitute.org/api',
            json={
                "query": variant_detail_query,
                "variables": variant_detail_variables
                },
            headers=headers)
        parse = json.loads(response.text)
        variant_details.append(parse['data']['variant'])
        print('fetched', this_variant)
        count += 1
        if (count > 1):
            break
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

# def flatten_consequences_for_selected_transcript(variants):
#     for variant in variants:
#         try:
#             for f in variant['consequencesForSelectedTranscript']:
#                 variant[f] = variant['consequencesForSelectedTranscript'][f]
#             del variant['consequencesForSelectedTranscript']
#         except:
#             pass
#         del variant['sortedTranscriptConsequences']
#     return variants

def main():
    dataset = "gnomad_r2_1_non_cancer"

    # gather brca1 data from gnomad
    brca1_exonic_variants = transcript_to_variants("ENST00000357654", dataset)
    brca1_intronic_variants = gene_to_region_variants("BRCA1", dataset)
    brca1_variants = brca1_intronic_variants | brca1_exonic_variants

    # gather brca2 data from gnomad
    brca2_exonic_variants = transcript_to_variants("ENST00000351666", dataset)
    brca2_intronic_variants = gene_to_region_variants("BRCA2", dataset)
    brca2_variants = brca2_intronic_variants | brca2_exonic_variants
    
    # combine brca1 and brca2 datasets into a pandas dataframe
    brca12_variants = brca1_variants | brca2_variants
    brca12_variant_data = variant_set_to_variant_data(brca12_variants, dataset)


    # flatten data
    variants_with_flattened_populations = flatten_populations(brca12_variant_data)
    # variants_with_flattened_fields = flatten_consequences_for_selected_transcript(variants_with_flattened_populations)
    
    # convert to dataframe, compute allele frequencies, and normalize
    variants_df = pd.DataFrame.from_dict(variants_with_flattened_populations)
    variants_df['flags'] = variants_df['flags'].apply(', '.join)
    brca12_dataframe = pd.json_normalize(brca12_variant_data)
    df_with_allele_values = compile_allele_values(brca12_dataframe)
    stringified_df_with_allele_values = df_with_allele_values.replace(np.nan, '-', regex=True).replace('', '-', regex=True)
    pdb.set_trace()


if __name__ == "__main__":
    main()