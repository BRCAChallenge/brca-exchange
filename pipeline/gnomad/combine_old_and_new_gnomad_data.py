import pandas as pd
import json
import argparse
import logging


def parse_args():
        parser = argparse.ArgumentParser(description='Download gnomad data and convert to .tsv format.')
        parser.add_argument('-n', '--new', type=argparse.FileType('r'),
                            help='Input new processd tsv gnomad file using updated download/processing.')
        parser.add_argument('-p', '--previous', type=argparse.FileType('r'),
                            help='Input previous processed tsv gnomad file, eg gnomAD.clean.tsv from 8/19/19.')
        parser.add_argument('-o', '--output', type=argparse.FileType('w'),
                            help='Output concattenated tsv gnomad file.')
        parser.add_argument('-l', '--logfile', default='/tmp/download_gnomad_data.log')
        parser.add_argument('-v', '--verbose', action='count', default=False, help='determines logging')
        options = parser.parse_args()
        return options


def main():
    options = parse_args()
    output = options.output
    logfile = options.logfile
    new = options.new
    previous = options.previous

    if options.verbose:
        logging_level = logging.DEBUG
    else:
        logging_level = logging.CRITICAL

    logging.basicConfig(filename=logfile, filemode="w", level=logging_level)

    new_df = pd.read_csv(new, sep='\t')

    # remove transcript from hgvsc column
    new_df['hgvsc'] = new_df['hgvsc'].str.split(':').str[1]

    previous_df = pd.read_csv(previous, sep='\t')

    # drop obsolete columns and rename variantID
    previous_df_clean = previous_df.drop(columns=['genome_AFR_ac_hemi', 'genome_AMR_ac_hemi', 'genome_ASJ_ac_hemi',
                                                  'genome_EAS_ac_hemi', 'genome_FIN_ac_hemi', 'genome_NFE_ac_hemi',
                                                  'genome_OTH_ac_hemi', 'genome_SAS_ac_hemi', 'exome_AFR_ac_hemi',
                                                  'exome_AMR_ac_hemi', 'exome_ASJ_ac_hemi', 'exome_EAS_ac_hemi',
                                                  'exome_FIN_ac_hemi', 'exome_NFE_ac_hemi', 'exome_OTH_ac_hemi',
                                                  'exome_SAS_ac_hemi', 'consequence']).rename(columns={'variantId': 'variant_id'})

    concattenated_df = pd.concat([new_df, previous_df_clean])
    concattenated_df.to_csv(output, sep='\t', index=False, na_rep='-')


if __name__ == "__main__":
    main()
