from django.core.management.base import BaseCommand, CommandError
from django.db import connection
from data.models import Variant, Report
from django.db import transaction
from data.utilities import update_materialized_view



EMPTY = '-'


class Command(BaseCommand):
    help = 'Split existing allele frequencies into appropriate columns'


    @transaction.atomic
    def handle(self, *args, **options):
        for Obj in [Variant, Report]:
            objs = Obj.objects.all()
            for obj in objs:
                if obj.Minor_allele_frequency_ESP_percent and obj.Minor_allele_frequency_ESP_percent != EMPTY:
                    eaAlleleFrequency = EMPTY
                    aaAlleleFrequency = EMPTY
                    alleleFrequency = EMPTY
                    maf = obj.Minor_allele_frequency_ESP_percent.split(',')

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

        update_materialized_view()

        print "Done!"
