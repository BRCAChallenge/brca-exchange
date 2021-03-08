#!/usr/bin/env python
import requests
import json
import numpy as np
import pandas as pd
import argparse
import logging
import time
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


def transcript_to_variants(transcript_id, dataset, reference):
    """
    Given a transcript, return the list of variants that map to the exons
    of the transcript, and were observed in samples from the indicated
    dataset.
    """
    fmt_graphql = """
    {
      transcript(transcript_id: "%s", reference_genome: %s) {
        variants(dataset: %s) {
        alt
          chrom
          pos
          ref
          variant_id
          transcript_consequence {
            transcript_id
            hgvsc
          }
          genome {
            ac
            an
            homozygote_count
            populations {
              id
              ac
              an
              homozygote_count
           } 
          }
          flags
      lof
      lof_filter
      lof_flags
      hgvs
    } 
      }
    }
    """
    req_variantlist = {
        "query": fmt_graphql % (transcript_id, reference, dataset),
        "variables": {}
    }
    response = fetch(req_variantlist)
    return response["data"]["transcript"]["variants"]


def gene_to_coords(gene_id, reference, max_retries=5):
    """                                                                         
    Given a gene symbol, return the coordinates.                                
    """
    graphql_query = """                                                         
    {                                                                           
        gene(gene_symbol: "%s", reference_genome: %s) {                         
            chrom                                                               
            start                                                               
            stop                                                                
        }                                                                       
    }"""
    graphql_request = {
        "query": graphql_query % (gene_id, reference),
        "variables": {}
    }
    retries = 0
    while retries < max_retries:
        response = fetch(graphql_request)
        gene_data = response["data"]["gene"]
        if gene_data is not None:
            return gene_data
        retries += 1
    logging.info(f'Request for gene coords {gene_id} failed')
    print(f'Request for gene coords {gene_id} failed')
    return None



def coords_to_variants(chrom, start, stop, dataset, reference):
    region_query = """
        {   region(chrom: "%s", start: %s, stop: %s, reference_genome: %s) {
            variants(dataset: %s) {
              variant_id
            }
          }
        }
    """
    r_q = region_query % (chrom, start, stop, reference, dataset)
    response = requests.post(
        'http://gnomad.broadinstitute.org/api',
        json={"query": r_q, "variables": {}},
        headers={"content-type": "application/json"})
    max_retries=5
    retries = 0
    while retries < max_retries:
        try:
            parse = json.loads(response.text)
            variants = parse['data']['region']['variants']
            if variants is not None:
                return variants
            else:
                continue
        except Exception as e:
            logging.info(e)
            print(f"coords to variants error {e}")
            retries += 1
            time.sleep(3)
            continue
    logging.info(f'Request for intronic variants non_cancer failed 5 attempts')


def fetch_data_for_one_variant(variant_id, dataset, max_retries=5):
    variant_query = """
            {   variant(variantId: "%s", dataset: %s) {
                variant_id
                chrom
                pos
                ref
                alt
                genome {
                  ac
                  an
                  af
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
                flags
                transcript_consequences {
                  hgvsc
                  lof
                  lof_flags
                  lof_filter
                  major_consequence
                  transcript_id
                }
              }
            }
      """
    v_q = variant_query % (variant_id, dataset)
    retries = 0
    while retries < max_retries:
        try:
            # https://stackoverflow.com/questions/49064398/requests-exceptions-chunkedencodingerror-connection-broken-incompleteread0
            response = requests.post(
                'http://gnomad.broadinstitute.org/api',
                json={ "query": v_q,
                       "variables": {}},
                headers={"content-type": "application/json"})
            if response.ok:
                time.sleep(0.1)
                parse = json.loads(response.text)
                print(f'success for variant {variant_id}')
                return(parse['data']['variant'])
            else:
                print('rate limit reached, waiting 60 seconds...')
                time.sleep(60)
                continue
        except json.decoder.JSONDecodeError:
            logging.info(f'Request for variant {variant_id} failed, kicking off retry #{retries}')
            retries += 1
            time.sleep(3)
            continue
        except KeyError as e:
            print(f"keyerror for variant {variant_id}, error is {e}")
            retries += 1
            time.sleep(300)
            continue
        except Exception as e:
            print(f"unknown error: {e}, sleeping for a couple minutes...")
            retries += 1
            time.sleep(300)
            continue
    print(f"failed to get data for {variant_id}")
    logging.info(f"failed to get data for {variant_id}")
    return None
    

def gene_to_region_variants(gene, dataset, reference):
    """
    Given a gene name, return the list of variants via a region
    query.  These will mostly be the intronic variants.
    """
    coords = gene_to_coords(gene, reference)
    variantList = coords_to_variants(coords["chrom"], coords["start"],
                                     coords["stop"], dataset, reference)
    return variantList


def getVariants(transcript, gene, dataset, reference):
    exonic_set = set()
    variant_ids = []
    exonic_variants_non_cancer = transcript_to_variants(transcript, dataset, reference)
    for variant in exonic_variants_non_cancer:
        variant_id = variant['variant_id']
        if variant_id not in variant_ids:
            variant = fetch_data_for_one_variant(variant_id, dataset)
            exonic_set.add(json.dumps(variant, sort_keys=True))
            variant_ids.append(variant_id)
    retries = 0
    max_retries = 5
    while retries < max_retries:
        intronic_set = set()
        variant_ids = []
        intronic_variants_non_cancer = gene_to_region_variants(gene, dataset, reference)
        if intronic_variants_non_cancer is None:
            retries += 1
            continue
        for variant in intronic_variants_non_cancer:
            variant_id = variant['variant_id']
            if variant_id not in variant_ids:
                variant = fetch_data_for_one_variant(variant_id, dataset)
                intronic_set.add(json.dumps(variant, sort_keys=True))
                variant_ids.append(variant_id)
        combined_set = exonic_set | intronic_set
        return combined_set
    logging.info(f'Request for intronic variants non_cancer failed 5 attempts')
    pdb.set_trace()


def convertSetToDict(mySet):
    myDict = dict()
    for elt in mySet:
        myElt = eval(elt.replace('null', 'None'))
        myDict[myElt['variant_id']] = myElt
    return myDict


def compile_allele_values(df):
    populations = ['afr', 'afr_XX', 'afr_XY', 'amr', 'amr_XX', 'amr_XY', 'asj', 'asj_XX', 'asj_XY', 'eas', 'eas_XX',
                   'eas_XY', 'fin', 'fin_XX', 'fin_XY', 'nfe', 'nfe_XX', 'nfe_XY', 'oth', 'oth_XX', 'oth_XY',
                   'sas', 'sas_XX', 'sas_XY', 'mid', 'mid_XX', 'mid_XY', 'ami', 'ami_XX', 'ami_XY', 'XX', 'XY']
    for population in populations:
        df['genome_' + population + '_af'] = calculate_frequency(df['genome_' + population + '_ac'], df['genome_' + population + '_an'])
    df['genome_af'] = calculate_frequency(df['genome_ac'], df['genome_an'])
    return df


def round_popmax(df):
    df['genome_popmax'] = pd.to_numeric(df['genome_popmax'], errors='coerce').apply(round_four_sigfigs)
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


def flatten(variant, field):
    for f in field:
        if f not in ['populations', 'filters']:
            if f == 'faf95':
                for p in field[f]:
                    variant['genome_' + p] = field[f][p]
            else:
                variant['genome_' + f] = field[f]
    populations = field['populations']
    for population in populations:
        name = population['id']
        if name.startswith('1kg') or 'hgdp' in name:
            continue
        keys = population.keys()
        for key in keys:
            if name != key:
                variant['genome_' + name + '_' + key] = population[key]
    return variant


def flatten_populations(variants):
    for variant in variants:
        genome = variant['genome']
        if genome:
            variant = flatten(variant, genome)
        del variant['genome']
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
        variant = variants[variant]
        transcriptConsequences = variant["transcript_consequences"]
        for transcript_hgvs in transcriptConsequences:
            if transcript_hgvs["transcript_id"] in transcripts:
                variant["hgvs"] = transcript_hgvs["hgvsc"]
                variant["transcript"] = transcript_hgvs["transcript_id"]
                del variant["transcript_consequences"]
                variants_in_expected_transcripts.append(variant)
                break
        if "hgvs" not in variant:
            print("Warning: variant %s falls outside the expected transcripts"
                  % variant["variant_id"])
    return variants_in_expected_transcripts


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-o', '--output', help='Output file base name')
    parser.add_argument('-l', '--logfile', help='Ouput logfile.')
    options = parser.parse_args()
    return options


def main():
    release = "r3"
    reference = "GRCh38"
    outputFile = parse_args().output + "_" + release + "_" + reference + ".tsv"
    log_file_path = parse_args().logfile
    logging.basicConfig(filename=log_file_path, filemode="w",
                        format=' %(asctime)s %(filename)-15s %(message)s')
    non_cancer_dataset = "gnomad_" + release + "_non_cancer"
    full_dataset = "gnomad_" + release
    brca1_transcript ="ENST00000357654"
    brca1_gene = "BRCA1"
    brca2_transcript = "ENST00000544455"
    brca2_gene = "BRCA2"

    # organize and combine variants
    brca1_variants_non_cancer = getVariants(brca1_transcript, brca1_gene, non_cancer_dataset, reference)
    variants_dict = convertSetToDict(brca1_variants_non_cancer)
    brca2_variants_non_cancer = getVariants(brca2_transcript, brca2_gene, non_cancer_dataset, reference)
    variants_dict.update(convertSetToDict(brca2_variants_non_cancer))

    # find hgvs, flatten, convert to dataframe, compute allele frequencies, and normalize
    variants_with_hgvs = find_correct_hgvs(variants_dict, (brca1_transcript, brca2_transcript))
    variants_with_flattened_populations = flatten_populations(variants_with_hgvs)
    variants_df = pd.json_normalize(variants_with_flattened_populations)
    variants_df['flags'] = variants_df['flags'].apply(', '.join)
    df_with_allele_values = compile_allele_values(variants_df)
    df_with_rounded_popmax = round_popmax(df_with_allele_values)
    stringified_df_with_allele_values = df_with_rounded_popmax.replace(np.nan, '-', regex=True).replace('', '-', regex=True)

    # output to .tsv
    stringified_df_with_allele_values.to_csv(outputFile, sep='\t', index=False)
    pdb.set_trace()

if __name__ == "__main__":
    main()