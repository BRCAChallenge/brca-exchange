from django.core.management.base import BaseCommand, CommandError
# from django.conf import settings
from django.db import connection
from data.models import Variant, DataRelease, Report
from argparse import FileType
import json
import csv
import psycopg2
from django.core.management import call_command


class Command(BaseCommand):
    help = 'Add a new variant release to the database'

    def add_arguments(self, parser):
        parser.add_argument('reports', type=FileType('r'), help='Reports to be added, in TSV format')

    def handle(self, *args, **options):
        reports_tsv = options['repots']

        reader = csv.reader(reports_tsv, dialect="excel-tab")
        header = reader.next()

        for row in reader:
            # split Source column into booleans
            row_dict = dict(zip(header, row))
            row_dict['Data_Release_id'] = 10
            row_dict['Variant'] = 10

            Report.objects.create_report(row_dict)
