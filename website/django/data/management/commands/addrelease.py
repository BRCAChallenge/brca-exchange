from django.core.management.base import BaseCommand, CommandError
from django.conf import settings
from data.models import Variant, DataRelease, ChangeType
from argparse import FileType
import json
import csv
import psycopg2


class Command(BaseCommand):
    help = 'Add a new variant release to the database'

    def add_arguments(self, parser):
        parser.add_argument('variants', type=FileType('r'), help='Variants to be added, in TSV format')
        parser.add_argument('notes', type=FileType('r'), help='Release notes and metadata, in JSON format')
        parser.add_argument('deletions', nargs='?', default=None, type=FileType('r'),
                            help='Deleted variants, in TSV format, same schema sans change_type')

    def update_autocomplete_words(self):
        dbname = settings.DATABASES['default']['NAME']
        user = settings.DATABASES['default']['USER']

        # Connect to an existing database
        conn = psycopg2.connect("dbname=" + dbname + " user=" + user)

        # Open a cursor to perform database operations
        cur = conn.cursor()

        # Drop words table and recreate with latest data
        cur.execute("""
            DROP TABLE IF EXISTS words;
            CREATE TABLE words AS SELECT DISTINCT left(word, 300) as word FROM (
            SELECT regexp_split_to_table(lower("Genomic_Coordinate_hg38"), '[\s|''"]') as word from variant UNION
            SELECT regexp_split_to_table(lower("Genomic_Coordinate_hg37"), '[\s|''"]') as word from variant UNION
            SELECT regexp_split_to_table(lower("Genomic_Coordinate_hg36"), '[\s|''"]') as word from variant UNION
            SELECT regexp_split_to_table(lower("Condition_category_ENIGMA"), '[\s|''"]') as word from variant UNION
            SELECT regexp_split_to_table(lower("Clinical_significance_ENIGMA"), '[\s|''"]') as word from variant UNION
            SELECT regexp_split_to_table(lower("Gene_Symbol"), '[\s|''"]') as word from variant UNION
            SELECT regexp_split_to_table(lower("Reference_Sequence"), '[\s|''"]') as word from variant UNION
            SELECT regexp_split_to_table(lower("HGVS_cDNA"), '[\s|:''"]') as word from variant UNION
            SELECT regexp_split_to_table(lower("BIC_Nomenclature"), '[\s|''"]') as word from variant UNION
            SELECT regexp_split_to_table(lower("HGVS_Protein"), '[\s|''"]') as word from variant
            )
            AS combined_words;

            CREATE INDEX words_idx ON words(word text_pattern_ops);
        """)

        # Make the changes to the database persistent
        conn.commit()

        # Close communication with the database
        cur.close()
        conn.close()

    def handle(self, *args, **options):
        variants_tsv = options['variants']
        notes = json.load(options['notes'])
        deletions_tsv = options['deletions']

        notes['sources'] = ', '.join(notes['sources'])
        release_id = DataRelease.objects.create(**notes).id

        reader = csv.reader(variants_tsv, dialect="excel-tab")
        header = reader.next()


        change_types = {ct['name']: ct['id'] for ct in ChangeType.objects.values()}

        for row in reader:
            # split Source column into booleans
            row_dict = dict(zip(header, row))
            if 'change_type' in row_dict and row_dict['change_type']:
                for source in row_dict['Source'].split(','):
                    row_dict['Variant_in_' + source] = True
                row_dict['Data_Release_id'] = release_id
                row_dict['Change_Type_id'] = change_types[row_dict.pop('change_type')]
                # use cleaned up genomic coordinates
                row_dict['Genomic_Coordinate_hg38'] = row_dict.pop('pyhgvs_Genomic_Coordinate_38')
                row_dict['Genomic_Coordinate_hg37'] = row_dict.pop('pyhgvs_Genomic_Coordinate_37')
                row_dict['Genomic_Coordinate_hg36'] = row_dict.pop('pyhgvs_Genomic_Coordinate_36')
                row_dict['Hg37_Start'] = row_dict.pop('pyhgvs_Hg37_Start')
                row_dict['Hg37_End'] = row_dict.pop('pyhgvs_Hg37_End')
                row_dict['Hg36_Start'] = row_dict.pop('pyhgvs_Hg36_Start')
                row_dict['Hg36_End'] = row_dict.pop('pyhgvs_Hg36_End')
                row_dict['HGVS_cDNA'] = row_dict.pop('pyhgvs_cDNA')
                row_dict['HGVS_Protein'] = row_dict.pop('pyhgvs_Protein')

                Variant.objects.create_variant(row_dict)

        # deleted variants
        if (deletions_tsv):
            reader = csv.reader(deletions_tsv, dialect="excel-tab")
            header = reader.next()
            for row in reader:
                # split Source column into booleans
                row_dict = dict(zip(header, row))
                for source in row_dict['Source'].split(','):
                    row_dict['Variant_in_' + source] = True
                row_dict['Data_Release_id'] = release_id
                # remove change type property, Variant only has Change_Type_id property
                row_dict.pop('change_type')
                row_dict['Change_Type_id'] = change_types['deleted']
                # use cleaned up genomic coordinates
                row_dict['Genomic_Coordinate_hg38'] = row_dict.pop('pyhgvs_Genomic_Coordinate_38')
                row_dict['Genomic_Coordinate_hg37'] = row_dict.pop('pyhgvs_Genomic_Coordinate_37')
                row_dict['Genomic_Coordinate_hg36'] = row_dict.pop('pyhgvs_Genomic_Coordinate_36')
                row_dict['Hg37_Start'] = row_dict.pop('pyhgvs_Hg37_Start')
                row_dict['Hg37_End'] = row_dict.pop('pyhgvs_Hg37_End')
                row_dict['Hg36_Start'] = row_dict.pop('pyhgvs_Hg36_Start')
                row_dict['Hg36_End'] = row_dict.pop('pyhgvs_Hg36_End')
                row_dict['HGVS_cDNA'] = row_dict.pop('pyhgvs_cDNA')
                row_dict['HGVS_Protein'] = row_dict.pop('pyhgvs_Protein')

                Variant.objects.create_variant(row_dict)

        self.update_autocomplete_words()
