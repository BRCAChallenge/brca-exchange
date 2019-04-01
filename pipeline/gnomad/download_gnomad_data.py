import requests
import pandas
import json
import argparse
import logging
import pdb

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

    variants_with_flattened_populations = flatten_populations(variants)

    df = pandas.DataFrame.from_dict(variants_with_flattened_populations)

    df.to_csv(output, sep='\t', index=False, na_rep='-')


def flatten_populations(variants):
    for variant in variants:
        populations = variant['populations']
        for pop in populations:
            name = pop['id']
            keys = pop.keys()
            for key in keys:
                if name != key:
                    variant[name + '_' + key] = pop[key]
        del variant['populations']
    return variants


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
        variants(dataset:gnomad_r2_1_non_cancer) {
            alt
            chrom
            pos
            ref
            variantId
            xpos
            ac
            ac_hemi
            ac_hom
            an
            af
            consequence
            consequence_in_canonical_transcript
            datasets
            filters
            flags
            hgvs
            hgvsc
            hgvsp
            rsid
            populations {
                id
                ac
                an
                ac_hemi
                ac_hom
            }
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

