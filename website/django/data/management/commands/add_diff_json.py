from django.core.management.base import BaseCommand, CommandError
# from django.conf import settings
from django.db import connection
from data.models import Variant, VariantDiff 
from argparse import FileType
import json
import psycopg2


class Command(BaseCommand):
    help = 'Add a new variant release to the database'

    def add_arguments(self, parser):
        parser.add_argument('release', type=int, help='release id')
        parser.add_argument('diffJSON', type=FileType('r'), help='JSON diff file')

    def handle(self, *args, **options):
        diff = json.load(options['diffJSON'])
        release_id = options['release']
        for key in diff:
            try:
                variant_instance = Variant.objects.filter(Data_Release_id = release_id).filter(Genomic_Coordinate_hg38 = key).get()
                VariantDiff.objects.create(id = variant_instance, diff = diff[key])
            except Variant.DoesNotExist:
                print "Error adding Variant Diff: Variant", key, "in release", release_id, "not found."
