from django.core.management.base import BaseCommand, CommandError
from django.db import connection
from data.models import Variant, Report, VariantDiff
from django.db import transaction
from math import floor, log10
from data.utilities import update_materialized_view



EMPTY = '-'

FIELDS_TO_ROUND = ["Allele_frequency_ExAC", "Allele_frequency_AFR_ExAC", "Allele_frequency_AMR_ExAC", "Allele_frequency_EAS_ExAC", "Allele_frequency_FIN_ExAC", "Allele_frequency_NFE_ExAC", "Allele_frequency_OTH_ExAC", "Allele_frequency_SAS_ExAC"]


class Command(BaseCommand):
    help = 'Rounds existing allele frequencies for ExAC fields to correct number of significant figures.'

    def round_sigfigs(self, num, sig_figs):
        if num != 0:
            return round(num, -int(floor(log10(abs(num))) - (sig_figs - 1)))
        else:
            return 0  # Can't take the log of 0

    @transaction.atomic
    def handle(self, *args, **options):
        VariantDiffs = VariantDiff.objects.all()
        for diffs in VariantDiffs:
            diff = diffs.diff
            for df in diff:
                if df['field'] in FIELDS_TO_ROUND:
                    val = df['added']
                    if val is not None and val != EMPTY:
                        try:
                            df['added'] = str(self.round_sigfigs(float(val), 3))
                        except ValueError:
                            print val
            diffs.save()
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

        update_materialized_view()

        print "Done!"
