from django.core.management.base import BaseCommand, CommandError
from data.models import Variant, DataRelease, ChangeType
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

        change_types = {ct['name']:ct['id'] for ct in ChangeType.objects.values()}

        for row in reader:
            # split Source column into booleans
            row_dict = dict(zip(header, row))
            for source in row_dict['Source'].split(','):
                row_dict['Variant_in_' + source] = True
            row_dict['Data_Release_id'] = release_number
            previous_version = Variant.objects.filter(Genomic_Coordinate_hg38 = row_dict['Genomic_Coordinate_hg38']).order_by('Data_Release_id').values()
            if previous_version:
                previous_version = previous_version[0]
                if row_dict['Pathogenicity_Default'] != previous_version['Pathogenicity_Default']:
                    row_dict['Change_Type_id'] = change_types['major']
                else:
                    for field in header:
                        if field not in row_dict:
                            print "Field [", field, "] not present in new data"
                        else:
                            new = row_dict[field]
                            old = previous_version[field]
                            if new != old:
                                print "Changed: [", field, "]", old, "->", new
                                row_dict['Change_Type_id'] = change_types['minor']
            else:
                row_dict['Change_Type_id'] = change_types['added']

            Variant.objects.create_variant(row_dict)
