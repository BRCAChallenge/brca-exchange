from django.core.management.base import BaseCommand, CommandError
# from django.conf import settings
from django.db import connection, transaction
from data.models import Variant, DataRelease, ChangeType, Report
from argparse import FileType
import json
import csv
import psycopg2
from django.core.management import call_command
from data.utilities import update_autocomplete_words


OLD_MAF_ESP_FIELD_NAMES = ['Minor_allele_frequency_ESP', 'Minor_allele_frequency_ESP_percent']


class Command(BaseCommand):
    help = 'Add a new variant release to the database'

    def add_arguments(self, parser):
        parser.add_argument('variants', type=FileType('r'), help='Variants to be added, in TSV format')
        parser.add_argument('reports', type=FileType('r'), help='Reports to be added, in TSV format')
        parser.add_argument('notes', type=FileType('r'), help='Release notes and metadata, in JSON format')
        parser.add_argument('deletions', nargs='?', default=None, type=FileType('r'),
                            help='Deleted variants, in TSV format, same schema sans change_type')
        parser.add_argument('diffJSON', help='JSON diff file')

    @transaction.atomic
    def handle(self, *args, **options):
        variants_tsv = options['variants']
        reports_tsv = options['reports']
        notes = json.load(options['notes'])
        deletions_tsv = options['deletions']
        diff_json = options['diffJSON']

        sources = notes['sources']

        notes['sources'] = ', '.join(sources)

        release_id = DataRelease.objects.create(**notes).id

        reports_reader = csv.reader(reports_tsv, dialect="excel-tab")
        reports_header = reports_reader.next()
        reports_dict = self.build_report_dictionary_by_source(reports_reader, reports_header, sources)

        variants_reader = csv.reader(variants_tsv, dialect="excel-tab")
        variants_header = variants_reader.next()

        change_types = {ct['name']: ct['id'] for ct in ChangeType.objects.values()}

        for row in variants_reader:
            # split Source column into booleans
            row_dict = dict(zip(variants_header, row))
            if 'change_type' in row_dict and row_dict['change_type']:
                row_dict = self.update_variant_values_for_insertion(row_dict, release_id, change_types)
                variant = Variant.objects.create_variant(row_dict)
                self.create_and_associate_reports_to_variant(variant, reports_dict, sources, release_id)

        # deleted variants
        if (deletions_tsv):
            deletions_reader = csv.reader(deletions_tsv, dialect="excel-tab")
            deletions_header = deletions_reader.next()
            for row in deletions_reader:
                # split Source column into booleans
                row_dict = dict(zip(deletions_header, row))
                row_dict = self.update_variant_values_for_insertion(row_dict, release_id, change_types, True)
                variant = Variant.objects.create_variant(row_dict)

        update_autocomplete_words()

        # update materialized view of current variants
        with connection.cursor() as cursor:
            cursor.execute("REFRESH MATERIALIZED VIEW currentvariant")

        # calls django/data/management/commands/add_diff_json to add diff to db
        call_command('add_diff_json', str(release_id), diff_json)

    def update_variant_values_for_insertion(self, row_dict, release_id, change_types, set_change_type_to_none=False):
        for source in row_dict['Source'].split(','):
            row_dict['Variant_in_' + source] = True
        row_dict['Data_Release_id'] = release_id
        if set_change_type_to_none is True:
            row_dict.pop('change_type', None)
            row_dict['Change_Type_id'] = change_types['deleted']
        else:
            row_dict['Change_Type_id'] = change_types[row_dict.pop('change_type')]

        # use cleaned up genomic coordinates and other values
        row_dict['Genomic_Coordinate_hg38'] = row_dict.pop('pyhgvs_Genomic_Coordinate_38')
        row_dict['Genomic_Coordinate_hg37'] = row_dict.pop('pyhgvs_Genomic_Coordinate_37')
        row_dict['Genomic_Coordinate_hg36'] = row_dict.pop('pyhgvs_Genomic_Coordinate_36')
        row_dict['Hg37_Start'] = row_dict.pop('pyhgvs_Hg37_Start')
        row_dict['Hg37_End'] = row_dict.pop('pyhgvs_Hg37_End')
        row_dict['Hg36_Start'] = row_dict.pop('pyhgvs_Hg36_Start')
        row_dict['Hg36_End'] = row_dict.pop('pyhgvs_Hg36_End')
        row_dict['HGVS_cDNA'] = row_dict.pop('pyhgvs_cDNA')
        row_dict['HGVS_Protein'] = row_dict.pop('pyhgvs_Protein')

        # Denote percentage in field name, two different fieldnames were used previously so both are handled below
        for oldName in OLD_MAF_ESP_FIELD_NAMES:
            if oldName in row_dict:
                row_dict['Minor_allele_frequency_percent_ESP'] = row_dict.pop(oldName)

        return row_dict

    def build_report_dictionary_by_source(self, reports_reader, reports_header, sources):
        reports_dict = {}
        for row in reports_reader:
            report = dict(zip(reports_header, row))
            source = report['Source']
            bx_id = report['BX_ID_' + source]
            if source not in reports_dict:
                reports_dict[source] = {}
            reports_dict[source][bx_id] = report
        return reports_dict

    def create_and_associate_reports_to_variant(self, variant, reports_dict, sources, release_id):
        for source in sources:
            if source == "Bic":
                source = "BIC"
            elif source == "1000 Genomes":
                source = "1000_Genomes"
            elif source == "ExUV":
                source = "exLOVD"
            bx_id_field = "BX_ID_" + source
            if not self.is_empty(getattr(variant, bx_id_field)):
                bx_ids = getattr(variant, bx_id_field).split(',')
                for bx_id in bx_ids:
                    self.create_and_associate_report_to_variant(bx_id, source, reports_dict, variant, release_id)

    def create_and_associate_report_to_variant(self, bx_id, source, reports_dict, variant, release_id):
        report = reports_dict[source][bx_id]
        report['Data_Release_id'] = release_id
        report['Variant'] = variant

        # Denote percentage in field name, two different fieldnames were used previously so both are handled below
        for oldName in OLD_MAF_ESP_FIELD_NAMES:
            if oldName in report:
                report['Minor_allele_frequency_percent_ESP'] = report.pop(oldName)

        Report.objects.create_report(report)

    def is_empty(self, value):
        return value is None or value is '' or value is '-'
