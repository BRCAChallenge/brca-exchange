from django.core.management.base import BaseCommand, CommandError
from django.db import connection
from data.models import Variant, VariantDiff, Report, ReportDiff
from argparse import FileType
import json
import psycopg2
from tqdm import tqdm

class Command(BaseCommand):
    help = 'Add diff information to variants'

    def add_arguments(self, parser):
        parser.add_argument('release', type=int, help='release id')
        parser.add_argument('diffJSON', type=FileType('r'), help='JSON diff file')
        parser.add_argument('reportsDiffJSON', type=FileType('r'), help='Reports JSON diff file')

    def handle(self, *args, **options):
        diff = json.load(options['diffJSON'])
        reports_diff = json.load(options['reportsDiffJSON'])
        release_id = options['release']

        print "Creating variant diffs..."
        for key in tqdm(diff, total=len(diff)):
            try:
                variant_instance = Variant.objects.filter(Data_Release_id=release_id).filter(Genomic_Coordinate_hg38=key).get()
                VariantDiff.objects.create(variant=variant_instance, diff=diff[key])
            except Variant.DoesNotExist:
                tqdm.write("Error adding Variant Diff: Variant %s in release %s not found" % (self._encode(key), release_id))

        print "Creating report diffs..."
        for key in tqdm(reports_diff, total=len(reports_diff)):

            # Only handles ClinVar and LOVD reports for now
            if "SCV" in key:
                # handle clinvar reports
                try:
                    report_instance = Report.objects.filter(Data_Release_id=release_id).filter(SCV_ClinVar=key).get()
                    ReportDiff.objects.create(report=report_instance, report_diff=reports_diff[key])
                except Report.DoesNotExist:
                    tqdm.write("Error adding Report Diff: Report %s in release %s not found" % (self._encode(key), release_id))
            else:
                # handle lovd reports
                try:
                    report_instance = Report.objects.filter(Data_Release_id=release_id).filter(Submission_ID_LOVD=key).get()
                    ReportDiff.objects.create(report=report_instance, report_diff=reports_diff[key])
                except Report.DoesNotExist:
                    tqdm.write("Error adding Report Diff: Report %s in release %s not found" % (self._encode(key), release_id))

    def _encode(self, s):
        return unicode(s).encode('utf-8')
