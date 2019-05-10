import requests
import pandas
import json
import argparse
import logging


def parse_args():
        parser = argparse.ArgumentParser(description='Download gnomad data and convert to .tsv format.')
        parser.add_argument('-o', '--output', type=argparse.FileType('w'),
                            help='Ouput TSV file result.')
        parser.add_argument('-l', '--logfile', default='/tmp/download_gnomad_data.log')
        parser.add_argument('-v', '--verbose', action='count', default=False, help='determines logging')
        options = parser.parse_args()
        return options


def main():
    options = parse_args()
    output = options.output
    logfile = options.logfile

    if options.verbose:
        logging_level = logging.DEBUG
    else:
        logging_level = logging.CRITICAL

    logging.basicConfig(filename=logfile, filemode="w", level=logging_level)

    variants_brca1 = generate_data("BRCA1")
    variants_brca2 = generate_data("BRCA2")

    variants = variants_brca1 + variants_brca2

    normalized_variants_df = normalize_variants(variants)

    normalized_variants_df.to_csv(output, sep='\t', index=False, na_rep='-')


def flatten_populations(variants):
    for variant in variants:
        genome = variant['genome']
        exome = variant['exome']
        if genome:
            for field in genome:
                if field not in ['populations', 'filters']:
                    variant['genome_' + field] = genome[field]
            genome_populations = genome['populations']
            for pop in genome_populations:
                name = pop['id']
                keys = pop.keys()
                for key in keys:
                    if name != key:
                        variant['genome_' + name + '_' + key] = pop[key]
        if exome:
            for field in exome:
                if field not in ['populations', 'filters']:
                    variant['exome_' + field] = exome[field]
            exome_populations = exome['populations']
            for pop in exome_populations:
                name = pop['id']
                keys = pop.keys()
                for key in keys:
                    if name != key:
                        variant['exome_' + name + '_' + key] = pop[key]
        del variant['genome']
        del variant['exome']
    return variants


def normalize_variants(variants):
    variants_with_flattened_populations = flatten_populations(variants)
    variants_df = pandas.DataFrame.from_dict(variants_with_flattened_populations)
    variants_df['flags'] = variants_df['flags'].apply(', '.join)
    return variants_df


def build_query(gene):
    return """{
        gene(gene_name: "%s") {
            _id
            omim_description
            gene_id
            omim_accession
            chrom
            strand
            full_gene_name
            gene_name_upper
            other_names
            canonical_transcript
            start
            stop
            xstop
            xstart
            gene_name
            variants(dataset: gnomad_r2_1_non_cancer) {
                alt
                chrom
                pos
                ref
                variantId
                xpos
                genome {
                    ac
                    ac_hemi
                    ac_hom
                    an
                    af
                    filters
                    populations {
                      id
                      ac
                      an
                      ac_hemi
                      ac_hom
                    }
                }
                exome {
                    ac
                    ac_hemi
                    ac_hom
                    an
                    af
                    filters
                    populations {
                        id
                        ac
                        an
                        ac_hemi
                        ac_hom
                    }
                }
                consequence
                consequence_in_canonical_transcript
                flags
                hgvs
                hgvsc
                hgvsp
                rsid
            }
        }
    }""" % (gene)


def generate_data(gene):
    query = build_query(gene)
    headers = { "content-type": "application/graphql" }
    response = requests.post('https://gnomad.broadinstitute.org/api', data=query, headers=headers)
    parsed_json = json.loads(response.text)
    return parsed_json['data']['gene']['variants']


if __name__ == "__main__":
    main()

