import csv
import json

from django.core.management.base import BaseCommand, CommandError
from django.db import connection, transaction
from data.models import Variant, VariantRepresentation
from argparse import FileType


class Command(BaseCommand):
    help = 'Add VR representations from a static file to the database'

    def add_arguments(self, parser):
        parser.add_argument('vr_reps', type=FileType('r'), help='VR representations to add, in TSV format')

    @transaction.atomic
    def handle(self, *args, **options):
        vr_rep_fp = csv.DictReader(options['vr_reps'], dialect='excel-tab')

        for row in vr_rep_fp:
            candidate = VariantRepresentation(
                Variant=Variant.objects.get(Genomic_Coordinate_hg38=row['Variant_id']),
                Description=json.loads(row['VR'])
            )
            candidate.save()
