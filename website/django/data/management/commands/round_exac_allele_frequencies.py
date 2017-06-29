from django.core.management.base import BaseCommand, CommandError
from django.db import connection
from data.models import Variant, Report
from django.db import transaction
import pdb


EMPTY = '-'

FIELDS_TO_ROUND = ["Allele_frequency_ExAC", "Allele_frequency_AFR_ExAC", "Allele_frequency_AMR_ExAC", "Allele_frequency_EAS_ExAC", "Allele_frequency_FIN_ExAC", "Allele_frequency_NFE_ExAC", "Allele_frequency_OTH_ExAC", "Allele_frequency_SAS_ExAC"]


class Command(BaseCommand):
    help = 'Rounds existing allele frequencies for ExAC fields to correct number of significant figures.'

    def update_autocomplete_words(self):
        # Drop words table and recreate with latest data
        with connection.cursor() as cursor:
            cursor.execute(
                """
                DROP MATERIALIZED VIEW IF EXISTS currentvariant;
                    CREATE MATERIALIZED VIEW currentvariant AS (
                        SELECT * FROM "variant" WHERE (
                            "id" IN ( SELECT DISTINCT ON ("Genomic_Coordinate_hg38") "id" FROM "variant" ORDER BY "Genomic_Coordinate_hg38" ASC, "Data_Release_id" DESC )
                        )
                    );
                """
            )

    def round_sigfigs(num, sig_figs):
        if num != 0:
            return round(num, -int(math.floor(math.log10(abs(num))) - (sig_figs - 1)))
        else:
            return 0  # Can't take the log of 0

    @transaction.atomic
    def handle(self, *args, **options):
        for Obj in [Variant, Report]:
            objs = Obj.objects.all()
            for obj in objs:
                for field in FIELDS_TO_ROUND:
                    if obj[field] != EMPTY:
                        pdb.set_trace()
                        obj[field] = round_sigfigs(float(obj[field], 3))
                obj.save()

        self.update_autocomplete_words()

        print "Done!"
