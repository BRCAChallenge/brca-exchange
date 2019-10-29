from collections import defaultdict

from django.core.management.base import BaseCommand, CommandError
# from django.conf import settings
from django.db import connection, transaction

from brca.site_settings import DATABASES
from data.models import Variant, DataRelease, ChangeType, Report, MupitStructure, InSilicoPriors
from argparse import FileType
import json
import csv
import psycopg2
from data.utilities import update_autocomplete_words
from tqdm import tqdm

from data.management.commands.add_diff_json import add_diffs


def get_num_lines(file_path):
    import mmap
    fp = open(file_path, "r+")
    buf = mmap.mmap(fp.fileno(), 0)
    lines = 0
    while buf.readline():
        lines += 1
    return lines


OLD_MAF_ESP_FIELD_NAMES = ['Minor_allele_frequency_ESP', 'Minor_allele_frequency_ESP_percent']

SKIP_VAR_INSERTION = True


class Command(BaseCommand):
    help = 'Add a new variant release to the database'

    def __init__(self):
        super(Command, self).__init__()
        self.previous_release_id = None  # filled in within handle()

    def add_arguments(self, parser):
        parser.add_argument('variants', type=FileType('r'), help='Variants to be added, in TSV format')
        parser.add_argument('notes', type=FileType('r'), help='Release notes and metadata, in JSON format')
        parser.add_argument('deletions', nargs='?', default=None, type=FileType('r'),
                            help='Deleted variants, in TSV format, same schema sans change_type')
        parser.add_argument('diffJSON', type=FileType('r'), help='JSON diff file')
        parser.add_argument('reports', type=FileType('r'), help='Reports to be added, in TSV format')
        parser.add_argument('removedReports', type=FileType('r'), help='Removed reports to be added, in TSV format')
        parser.add_argument('reportsDiffJSON', type=FileType('r'), help='Reports JSON diff file')

    @transaction.atomic
    def handle(self, *args, **options):
        variants_tsv = options['variants']
        reports_tsv = options['reports']
        removed_reports_tsv = options['removedReports']
        notes = json.load(options['notes'])
        deletions_tsv = options['deletions']
        diff_json = options['diffJSON']
        reports_diff_json = options['reportsDiffJSON']

        sources = notes['sources']

        notes['sources'] = ', '.join(sources)

        # Currently, the release name is just an ascending number starting at one for the first release.
        # To name the release we're adding, find the most recently added release and add 1 to its name.
        release_name = int(DataRelease.objects.all().order_by('-name')[0].name) + 1
        release_id = DataRelease.objects.create(name=release_name, **notes).id

        self.previous_release_id = DataRelease.objects.order_by('-id')[1].id

        print(("Creating new release with ID %d and name %s in db %s" % (release_id, release_name, DATABASES['default']['NAME'])))

        reports_dict = self.build_report_dictionary_by_source(reports_tsv, removed_reports_tsv)
        change_types = {ct['name']: ct['id'] for ct in list(ChangeType.objects.values())}
        mupit_structures = {ms['name']: ms['id'] for ms in list(MupitStructure.objects.values())}

        # TODO: split variants_tsv ahead of time into separate variant and insilicopriors lists
        # TODO: use django-postgres-copy to load these files, although we'll have to run update_variant_values_for_insertion on it beforehand

        self.insert_variants(
            variants_tsv, label="inserting variants", is_deletion=False, release_id=release_id,
            change_types=change_types, mupit_structures=mupit_structures, reports_dict=reports_dict, sources=sources
        )

        # deleted variants
        if deletions_tsv:
            self.insert_variants(
                deletions_tsv, label="applying deletions", is_deletion=True, release_id=release_id,
                change_types=change_types, mupit_structures=mupit_structures, reports_dict=reports_dict, sources=sources
            )

        # calls django/data/management/commands/add_diff_json to add diff to db
        # converted to a method call so that we don't launch a separate process that's not in the current transaction
        # call_command('add_diff_json', str(release_id), diff_json, reports_diff_json)
        add_diffs(diff_json, str(release_id), reports_diff_json)

        # raise Exception("terminating b/c i don't want to commit")

        update_autocomplete_words()

        # update materialized view of current variants
        with connection.cursor() as cursor:
            cursor.execute("REFRESH MATERIALIZED VIEW currentvariant")

        # raise Exception("terminating b/c i don't want to commit")

    def insert_variants(self, tsv_fp, label, is_deletion, release_id, change_types, mupit_structures, reports_dict, sources):
        """
        Given a file pointer to a TSV file containing variants, insert each variant + in-silico priors information
        into the database. If is_deletion is true, mark the variants we're inserting as deletions instead of
        new or changed variant information.
        :param tsv_fp: a handle to a tab-delimited file containing variant data
        :param label: the message to display to the user in the progress bar
        :param is_deletion: whether these variants should be treated as new/changed variants, or deletions
        :param release_id: the number of this release
        :param change_types: a mapping from change_type strings (e.g., "added_information") to ChangeType model IDs
        :param mupit_structures: a mapping from mupit strings (e.g., "4igk") to MupitStructure model IDs
        :param reports_dict: an aggregation of reports, returned by build_report_dictionary_by_source()
        :param sources: a list of sources included in this release, e.g. ["Bic", "ClinVar", "ESP", "ExAC"...]
        :return:
        """
        reader = csv.reader(tsv_fp, dialect="excel-tab")
        header = next(reader)

        # collect list of insilico prior column names so we can split those into a separate table, insilicopriors
        insilicopriors_cols = [x.name for x in InSilicoPriors._meta.get_fields() if x.name not in ['id', 'Variant']]

        # tqdm() displays a progress bar as we read elements from the csvreader, making this all slightly more pleasant
        for row in tqdm(reader, desc=label, total=get_num_lines(tsv_fp.name) - 1):
            # split Source column into booleans
            row_dict = dict(list(zip(header, row)))

            # transfer insilico prior columns from row_dict to insilico_dict
            insilico_dict = dict((col, row_dict.pop(col)) for col in insilicopriors_cols)

            # basically, we process every variant if we're loading deletions,
            # or if it's a non-deletion only variants that have a valid change type field
            if is_deletion or ('change_type' in row_dict and row_dict['change_type']):
                # remap certain variant columns to account for historical differences in the variant description
                row_dict = self.update_variant_values_for_insertion(
                    row_dict, release_id, change_types, mupit_structures,
                    set_ct_and_ms_to_none=is_deletion
                )

                # create the actual variant
                variant = Variant.objects.create_variant(row_dict)

                if not is_deletion:
                    self.create_and_associate_reports_to_variant(variant, reports_dict, sources, release_id, change_types)

                # also create an instance in insilicopriors linked to this variant
                InSilicoPriors.objects.create(Variant=variant, **insilico_dict)

    def update_variant_values_for_insertion(self, row_dict, release_id, change_types, mupit_structures, set_ct_and_ms_to_none=False):
        for source in row_dict['Source'].split(','):
            row_dict['Variant_in_' + source] = True
        row_dict['Data_Release_id'] = release_id
        if set_ct_and_ms_to_none is True:
            row_dict.pop('change_type', None)
            row_dict.pop('mupit_structure', None)
            row_dict['Change_Type_id'] = change_types['deleted']
            row_dict['Mupit_Structure_id'] = None
        else:
            row_dict['Change_Type_id'] = change_types[row_dict.pop('change_type')]
            mupit_structure = row_dict.pop('mupit_structure')
            if self.is_empty(mupit_structure):
                row_dict['Mupit_Structure_id'] = None
            else:
                row_dict['Mupit_Structure_id'] = mupit_structures[mupit_structure]

        # use cleaned up genomic coordinates and other values
        row_dict['Genomic_Coordinate_hg38'] = row_dict.pop('pyhgvs_Genomic_Coordinate_38')
        row_dict['Genomic_Coordinate_hg37'] = row_dict.pop('pyhgvs_Genomic_Coordinate_37')
        row_dict['Hg37_Start'] = row_dict.pop('pyhgvs_Hg37_Start')
        row_dict['Hg37_End'] = row_dict.pop('pyhgvs_Hg37_End')
        row_dict['HGVS_cDNA'] = row_dict.pop('pyhgvs_cDNA')
        row_dict['HGVS_Protein'] = row_dict.pop('pyhgvs_Protein')

        # Denote percentage in field name, two different fieldnames were used previously so both are handled below
        for oldName in OLD_MAF_ESP_FIELD_NAMES:
            if oldName in row_dict:
                row_dict['Minor_allele_frequency_percent_ESP'] = row_dict.pop(oldName)

        return row_dict

    def build_report_dictionary_by_source(self, reports_tsv, removed_reports_tsv):
        reports_reader = csv.reader(reports_tsv, dialect="excel-tab")
        reports_header = next(reports_reader)
        removed_reports_reader = csv.reader(removed_reports_tsv, dialect="excel-tab")
        removed_reports_header = next(removed_reports_reader)

        reports_dict = {
            'reports': defaultdict(dict), 'removed_reports': defaultdict(dict)
        }

        for row in reports_reader:
            report = dict(list(zip(reports_header, row)))
            if self.is_empty(report['change_type']):
                report['change_type'] = 'none'
            source = report['Source']
            bx_id = report['BX_ID_' + source]
            reports_dict['reports'][source][bx_id] = report

        # add removed reports to dict
        for row in removed_reports_reader:
            report = dict(list(zip(reports_header, row)))
            report['change_type'] = 'deleted'
            source = report['Source']
            bx_id = report['BX_ID_' + source]
            reports_dict['removed_reports'][source][bx_id] = report

        return reports_dict

    def create_and_associate_reports_to_variant(self, variant, reports_dict, sources, release_id, change_types):
        # Used to associate removed reports with variants because removed report bx_ids refer
        # to bx_ids from the previous release.\
        previous_version_of_variant_query = Variant.objects.filter(
            Data_Release_id=self.previous_release_id, Genomic_Coordinate_hg38=variant.Genomic_Coordinate_hg38
        )

        for source in sources:
            if source == "Bic":
                source = "BIC"
            elif source == "1000 Genomes":
                source = "1000_Genomes"
            elif source == "ExUV":
                source = "exLOVD"
            elif source == "Findlay BRCA1 Ring Function Scores":
                source = "Findlay_BRCA1_Ring_Function_Scores"
            bx_id_field = "BX_ID_" + source

            if not self.is_empty(getattr(variant, bx_id_field)):
                bx_ids = getattr(variant, bx_id_field).split(',')
                for bx_id in bx_ids:
                    self.create_and_associate_report_to_variant(bx_id, source, reports_dict, variant, release_id, change_types, removed=False)

            '''
            Associate removed reports with variant -- note that bx_ids are release specific,
            so we must check the previous version of a variant (if it exists) for its bx_ids and compare them
            to the removed report bx_ids (removed report bx_ids are from the previous release).
            '''
            if len(previous_version_of_variant_query) > 0:
                previous_version_of_variant = previous_version_of_variant_query[0]
                if not self.is_empty(getattr(previous_version_of_variant, bx_id_field)):
                    bx_ids = getattr(previous_version_of_variant, bx_id_field).split(',')
                    for bx_id in bx_ids:
                        if source in reports_dict['removed_reports'] and bx_id in reports_dict['removed_reports'][source]:
                            self.create_and_associate_report_to_variant(bx_id, source, reports_dict, variant, release_id, change_types, removed=True)

    def create_and_associate_report_to_variant(self, bx_id, source, reports_dict, variant, release_id, change_types, removed=False):
        if removed:
            report = reports_dict['removed_reports'][source][bx_id]
        else:
            report = reports_dict['reports'][source][bx_id]
        report['Data_Release_id'] = release_id
        report['Variant'] = variant
        report['Change_Type_id'] = change_types[report.pop('change_type')]

        # Denote percentage in field name, two different fieldnames were used previously so both are handled below
        for oldName in OLD_MAF_ESP_FIELD_NAMES:
            if oldName in report:
                report['Minor_allele_frequency_percent_ESP'] = report.pop(oldName)

        Report.objects.create_report(report)

    def is_empty(self, value):
        return value is None or value == '' or value == '-'
