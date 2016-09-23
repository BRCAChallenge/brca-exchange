from django.core.management.base import BaseCommand, CommandError
from data.models import Variant, DataRelease
from argparse import FileType
import json
import csv

class Command(BaseCommand):
    help = 'Add a new variant release to the database'

    def add_arguments(self, parser):
        parser.add_argument('variants', type=FileType('r'), help='Variants to be added, in TSV format')
        parser.add_argument('notes', type=FileType('r'), help='Release notes and metadata, in JSON format')

    def handle(self, *args, **options):
        variants_tsv = options['variants']
        notes = json.load(options['notes'])
        release_number = notes['release_number']

        #print "Variant Count:", sum(1 for line in variants_tsv)
        print "Release:", notes['release_number']
        print "Timestamp:", notes['timestamp']
        print "Comment:", notes['comment']

        DataRelease.objects.create(id=release_number, timestamp=notes['timestamp'], comment=notes['comment'])

        reader = csv.reader(variants_tsv, dialect="excel-tab")

        header = reader.next()

        variant_count = 0
        for row in reader:
            # split Source column into booleans
            row_dict = dict(zip(header, row))
            for source in row_dict['Source'].split(','):
                row_dict['Variant_in_' + source] = True
            row_dict['Data_Release_id'] = release_number
            Variant.objects.create_variant(row_dict)
            variant_count += 1;

        print variant_count, "variants added."
