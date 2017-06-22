from django.core.management.base import BaseCommand, CommandError
from django.db import connection
from data.models import Variant, Report
from django.db import transaction


EMPTY = '-'


class Command(BaseCommand):
    help = 'Split existing allele frequencies into appropriate columns'

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

    @transaction.atomic
    def handle(self, *args, **options):
        for Obj in [Variant, Report]:
            objs = Obj.objects.all()
            for obj in objs:
                if obj.Minor_allele_frequency_ESP and obj.Minor_allele_frequency_ESP != EMPTY:
                    eaAlleleFrequency = EMPTY
                    aaAlleleFrequency = EMPTY
                    alleleFrequency = EMPTY
                    maf = obj.Minor_allele_frequency_ESP.split(',')

                    if len(maf) > 2:
                        alleleFrequency = "%s" % (float(maf[2]) / 100)
                    if len(maf) > 1:
                        aaAlleleFrequency = "%s" % (float(maf[1]) / 100)
                    if len(maf) > 0:
                        eaAlleleFrequency = "%s" % (float(maf[0]) / 100)

                    obj.EA_Allele_Frequency_ESP = eaAlleleFrequency
                    obj.AA_Allele_Frequency_ESP = aaAlleleFrequency
                    obj.Allele_Frequency_ESP = alleleFrequency
                    obj.save()

        self.update_autocomplete_words()

        print "Done!"
