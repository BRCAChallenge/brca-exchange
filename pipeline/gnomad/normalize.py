import pandas
import json
import argparse
import logging
import pdb

def parse_args():
        parser = argparse.ArgumentParser(description='Download gnomad data and convert to .tsv format.')
        parser.add_argument('-i', '--input', type=argparse.FileType('r'),
                            help='Input json file.')
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
    variants = options.input

    if options.verbose:
        logging_level = logging.DEBUG
    else:
        logging_level = logging.CRITICAL

    logging.basicConfig(filename=logfile, filemode="w", level=logging_level)

    variants = json.load(variants)

    normalized_variants_df = normalize_variants(variants)

    normalized_variants_df.to_csv(output, sep='\t', index=False, na_rep='-')


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

def flatten_consequences_for_selected_transcript(variants):
    for variant in variants:
        try:
            for f in variant['consequencesForSelectedTranscript']:
                variant[f] = variant['consequencesForSelectedTranscript'][f]
            del variant['consequencesForSelectedTranscript']
        except:
            pdb.set_trace()
            pass
        del variant['sortedTranscriptConsequences']
    return variants      
            
def normalize_variants(variants):
    variants_with_flattened_populations = flatten_populations(variants)
    variants_with_flattened_fields = flatten_consequences_for_selected_transcript(variants_with_flattened_populations)
    variants_df = pandas.DataFrame.from_dict(variants_with_flattened_fields)
    variants_df['flags'] = variants_df['flags'].apply(', '.join)
    return variants_df


if __name__ == "__main__":
    main()


