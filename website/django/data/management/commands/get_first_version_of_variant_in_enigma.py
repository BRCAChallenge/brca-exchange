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
        with open('/tmp/first_versions_of_variants_in_enigma_barring_first_release_ammended.csv', 'wb') as csvfile:
            fieldnames = [
                            'Genomic_Coordinate_hg38',
                            'Gene_Symbol',
                            'HGVS_cDNA',
                            'Data_Release_id',
                            'Sources',
                            'Clinical_Significance_ClinVar_prior_to_Enigma_Classification',
                            'Date_Last_Updated_ClinVar_prior_to_Enigma_Classification',
                            'Submitter_ClinVar_prior_to_Enigma_Classification',
                            'SCV_ClinVar_prior_to_Enigma_Classification',
                            'Allele_Origin_ClinVar_prior_to_Enigma_Classification',
                            'Method_ClinVar_prior_to_Enigma_Classification',
                            'ENIGMA_comment_on_clinical_significance',
                            'Only_submitted_by_ENIGMA'
                        ]
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()
            variants = Variant.objects.all().exclude(Data_Release_id=1).order_by("Data_Release_id")
            current_variants = CurrentVariant.objects.all().exclude(Change_Type__name='deleted')
            for current_variant in current_variants:
                current_variant_keys.append(current_variant.Genomic_Coordinate_hg38)
            for variant in variants:
                if "ENIGMA" in variant.Source:
                    key = variant.Genomic_Coordinate_hg38
                    if key in current_variant_keys and key not in seen:
                        seen.append(key)
                        previous_release_id = variant.Data_Release_id - 1
                        try:
                            previous_version = Variant.objects.get(Genomic_Coordinate_hg38=key, Data_Release_id=previous_release_id)
                            Clinical_Significance_ClinVar_prior_to_Enigma_Classification = previous_version.Clinical_Significance_ClinVar
                            Date_Last_Updated_ClinVar_prior_to_Enigma_Classification = previous_version.Date_Last_Updated_ClinVar
                            Submitter_ClinVar_prior_to_Enigma_Classification = previous_version.Submitter_ClinVar
                            SCV_ClinVar_prior_to_Enigma_Classification = previous_version.SCV_ClinVar
                            Allele_Origin_ClinVar_prior_to_Enigma_Classification = previous_version.Allele_Origin_ClinVar
                            Method_ClinVar_prior_to_Enigma_Classification = previous_version.Method_ClinVar
                        except:
                            Clinical_Significance_ClinVar_prior_to_Enigma_Classification = "N/A"
                            Date_Last_Updated_ClinVar_prior_to_Enigma_Classification = "N/A"
                            Submitter_ClinVar_prior_to_Enigma_Classification = "N/A"
                            SCV_ClinVar_prior_to_Enigma_Classification = "N/A"
                            Allele_Origin_ClinVar_prior_to_Enigma_Classification = "N/A"
                            Method_ClinVar_prior_to_Enigma_Classification = "N/A"
                        Comment_on_clinical_significance_ENIGMA = variant.Comment_on_clinical_significance_ENIGMA
                        if variant.Data_Release_id > 1:
                            writer.writerow({
                                'Genomic_Coordinate_hg38': key.encode('utf-8'),
                                'Gene_Symbol': variant.Gene_Symbol.encode('utf-8'),
                                'HGVS_cDNA': variant.HGVS_cDNA.encode('utf-8'),
                                'Data_Release_id': variant.Data_Release_id,
                                'Sources': variant.Source.encode('utf-8'),
                                'Clinical_Significance_ClinVar_prior_to_Enigma_Classification': Clinical_Significance_ClinVar_prior_to_Enigma_Classification.encode('utf-8'),
                                'Date_Last_Updated_ClinVar_prior_to_Enigma_Classification': Date_Last_Updated_ClinVar_prior_to_Enigma_Classification.encode('utf-8'),
                                'Submitter_ClinVar_prior_to_Enigma_Classification': Submitter_ClinVar_prior_to_Enigma_Classification.encode('utf-8'),
                                'SCV_ClinVar_prior_to_Enigma_Classification': SCV_ClinVar_prior_to_Enigma_Classification.encode('utf-8'),
                                'Allele_Origin_ClinVar_prior_to_Enigma_Classification': Allele_Origin_ClinVar_prior_to_Enigma_Classification.encode('utf-8'),
                                'Method_ClinVar_prior_to_Enigma_Classification': Method_ClinVar_prior_to_Enigma_Classification.encode('utf-8'),
                                'ENIGMA_comment_on_clinical_significance': Comment_on_clinical_significance_ENIGMA.encode('utf-8'),
                                'Only_submitted_by_ENIGMA': variant.Source == "ENIGMA"
                            })

        print "Done!"
