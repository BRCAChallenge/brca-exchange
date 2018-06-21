from django.core.management.base import BaseCommand, CommandError
from django.db import connection
from data.models import Variant
from django.db import transaction
import csv
import pdb


class Command(BaseCommand):
    help = 'Creates a file that lists what release a variant was first found and whether it was only submitted by ENIGMA in that release.'

    @transaction.atomic
    def handle(self, *args, **options):
        seen = {}
        with open('/tmp/firstVersions.csv', 'wb') as csvfile:
            fieldnames = ['Genomic_Coordinate_hg38', 'Data_Release_id', 'Only_submitted_by_ENIGMA']
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()
            variants = Variant.objects.all().exclude(Data_Release_id='1').order_by("Data_Release_id")
            for variant in variants:
                key = variant.Genomic_Coordinate_hg38
                if key not in seen:
                    seen[key] = 1
                    writer.writerow({'Genomic_Coordinate_hg38': key, 'Data_Release_id': variant.Data_Release_id, 'Only_submitted_by_ENIGMA': variant.Source == "ENIGMA"})

        print "Done!"
