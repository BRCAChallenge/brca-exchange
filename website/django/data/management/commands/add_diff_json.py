from django.core.management.base import BaseCommand, CommandError
from django.db import connection
from data.models import Variant, VariantDiff, Report, ReportDiff
from argparse import FileType
import json
import psycopg2


class Command(BaseCommand):
    help = 'Add diff information to variants'

    def add_arguments(self, parser):
        parser.add_argument('release', type=int, help='release id')
        parser.add_argument('diffJSON', type=FileType('r'), help='JSON diff file')
        parser.add_argument('reportsDiffJSON', type=FileType('r'), help='Reports JSON diff file')

    def handle(self, *args, **options):
        diff = json.load(options['diffJSON'])
        report_diff = json.load(options['reportsDiffJSON'])
        release_id = options['release']
        for key in diff:
            try:
                variant_instance = Variant.objects.filter(Data_Release_id=release_id).filter(Genomic_Coordinate_hg38=key).get()
                VariantDiff.objects.create(variant=variant_instance, diff=diff[key])
            except Variant.DoesNotExist:
                print "Error adding Variant Diff: Variant", key, "in release", release_id, "not found."
        for key in reports_diff:
            # Only handles ClinVar and LOVD reports for now
            if "scv" in key.lowercase():
                # handle clinvar reports
                try:
                    report_instance = Report.objects.filter(Data_Release_id=release_id).filter(SCV_ClinVar=key).get()
                    ReportDiff.objects.create(report=report_instance, diff=reports_diff[key])
                except Report.DoesNotExist:
                    print "Error adding Report Diff: Report", key, "in release", release_id, "not found."
            else:
                # handle lovd reports
                try:
                    report_instance = Report.objects.filter(Data_Release_id=release_id).filter(DBID_LOVD=key).get()
                    ReportDiff.objects.create(report=report_instance, diff=reports_diff[key])
                except Report.DoesNotExist:
                    print "Error adding Report Diff: Report", key, "in release", release_id, "not found."
