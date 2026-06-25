from django.core.management.base import BaseCommand, CommandError
from django.db import connection, transaction
from data.models import Variant, VariantDiff, Report, ReportDiff
from argparse import FileType
import json
import psycopg2
from tqdm import tqdm

from django.db import connection
from io import StringIO

from data.utilities import Benchmark


class Command(BaseCommand):
    help = 'Add diff information to variants'

    def add_arguments(self, parser):
        parser.add_argument('release', type=int, help='release id')
        parser.add_argument('diffJSON', type=FileType('r'), help='JSON diff file')
        parser.add_argument('reportsDiffJSON', type=FileType('r'), help='Reports JSON diff file')

    @transaction.atomic
    def handle(self, *args, **options):
        diff = options['diffJSON']
        reports_diff = options['reportsDiffJSON']
        release_id = options['release']

        add_diffs(diff, release_id, reports_diff)


def add_diffs(diff_fp, release_id, reports_diff_fp):
    diff = json.load(diff_fp)
    reports_diff = json.load(reports_diff_fp)

    print("Creating variant diffs...")
    with Benchmark("variant diffs"):
        with connection.cursor() as cursor:
            cursor.execute("""create temporary table _var_diffs (key text, diff json)""")

            with cursor.connection.cursor() as psycon:
                # psycon = cursor.connection  # get the underlying psycopg2 handle so we can use copy_from()
                buf = StringIO("".join("%s\t%s\n" % (k, json.dumps(diff[k]).replace('\\', '\\\\')) for k in diff))
                psycon.copy_from(file=buf, table="_var_diffs")

            cursor.execute("""
            insert into data_variantdiff (variant_id, diff)
            select variant.id, _var_diffs.diff from _var_diffs
            inner join variant on _var_diffs.key = variant."Genomic_Coordinate_hg38"
            and variant."Data_Release_id"=%s
            -- on conflict DO NOTHING;
            """, [release_id])

    print("Creating report diffs...")
    with Benchmark("creating report diffs"):
        with connection.cursor() as cursor:
            cursor.execute("""create temporary table _report_diffs (key text, diff json)""")

            with cursor.connection.cursor() as psycon:
                # psycon = cursor.connection  # get the underlying psycopg2 handle so we can use copy_from()
                buf = StringIO("".join("%s\t%s\n" % (k, json.dumps(reports_diff[k]).replace('\\', '\\\\')) for k in reports_diff))
                psycon.copy_from(buf, table="_report_diffs")

            cursor.execute("""
            insert into data_reportdiff (report_id, report_diff)
            select report.id, _report_diffs.diff from _report_diffs
            inner join report on _report_diffs.key = (
              case when _report_diffs.key like 'SCV%%' then report."SCV_ClinVar" else report."Submission_ID_LOVD" end
            )
            and report."Data_Release_id"=%s
            -- on conflict DO NOTHING;
            """, [release_id])


def _encode(s):
    return str(s).encode('utf-8')
