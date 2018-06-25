from django.core.management.base import BaseCommand, CommandError
from django.db import connection
from data.models import Variant, CurrentVariant
from django.db import transaction
import csv


class Command(BaseCommand):
    help = 'Creates a file that lists what release a variant was first found and whether it was only submitted by ENIGMA in that release.'

    @transaction.atomic
    def handle(self, *args, **options):
        seen = []
        current_variant_keys = []
        with open('/tmp/first_versions_of_variants_barring_first_release_ammended.csv', 'wb') as csvfile:
            fieldnames = ['Genomic_Coordinate_hg38', 'Gene_Symbol', 'HGVS_cDNA', 'Data_Release_id', 'Sources', 'ClinVar_clinical_significance_at_time_of_Enigma_Classification', 'ENIGMA_comment_on_clinical_significance', 'Only_submitted_by_ENIGMA']
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()
            variants = Variant.objects.all().order_by("Data_Release_id")
            current_variants = CurrentVariant.objects.all().exclude(Change_Type__name='deleted')
            for current_variant in current_variants:
                current_variant_keys.append(current_variant.Genomic_Coordinate_hg38)
            for variant in variants:
                key = variant.Genomic_Coordinate_hg38
                if key in current_variant_keys and key not in seen:
                    seen.append(key)
                    if "ENIGMA" in variant.Source:
                        Clinical_Significance_ClinVar = variant.Clinical_Significance_ClinVar
                        Comment_on_clinical_significance_ENIGMA = variant.Comment_on_clinical_significance_ENIGMA
                    else:
                        Clinical_Significance_ClinVar = "N/A"
                        Comment_on_clinical_significance_ENIGMA = "N/A"
                    if variant.Data_Release_id > 1:
                        writer.writerow({
                            'Genomic_Coordinate_hg38': key,
                            'Gene_Symbol': variant.Gene_Symbol,
                            'HGVS_cDNA': variant.HGVS_cDNA,
                            'Data_Release_id': variant.Data_Release_id,
                            'Sources': variant.Source,
                            'ClinVar_clinical_significance_at_time_of_Enigma_Classification': Clinical_Significance_ClinVar,
                            'ENIGMA_comment_on_clinical_significance': Comment_on_clinical_significance_ENIGMA,
                            'Only_submitted_by_ENIGMA': variant.Source == "ENIGMA"
                        })

        print "Done!"
