from django.core.management.base import BaseCommand, CommandError
from django.db import connection
from data.models import Variant, Report
from django.db import transaction
from math import floor, log10


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

    def round_sigfigs(self, num, sig_figs):
        if num != 0:
            return round(num, -int(floor(log10(abs(num))) - (sig_figs - 1)))
        else:
            return 0  # Can't take the log of 0

    @transaction.atomic
    def handle(self, *args, **options):
        for Obj in [Variant, Report]:
            objs = Obj.objects.all()
            for obj in objs:
                for field in FIELDS_TO_ROUND:
                    val = getattr(obj, field)
                    if val is not None and val != EMPTY:
                        try:
                            setattr(obj, field, str(self.round_sigfigs(float(val), 3)))
                        except ValueError:
                            setattr(obj, field, EMPTY)
                obj.save()

        self.update_autocomplete_words()

        print "Done!"
