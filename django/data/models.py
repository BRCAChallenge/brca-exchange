from django.db import models


class VariantManager(models.Manager):
    def create_variant(self, row):
        return self.create(**row)


class Variant(models.Model):
    Variant_in_ENIGMA = models.BooleanField()
    Variant_in_ClinVar = models.BooleanField()
    Variant_in_1000_Genomes = models.BooleanField()
    Variant_in_ExAC = models.BooleanField()
    Variant_in_LOVD = models.BooleanField()
    Variant_in_BIC = models.BooleanField()
    Gene_symbol = models.TextField()
    Genomic_Coordinate = models.TextField()
    HGVS_genomic = models.TextField()
    Transcript_id = models.TextField()
    HGVS_cDNA = models.TextField()
    HGVS_protein = models.TextField()
    BIC_Nomenclature = models.TextField()
    Abbrev_AA_change = models.TextField()
    URL = models.TextField()
    Condition_ID_type = models.TextField()
    Condition_ID_value = models.TextField()
    Condition_category = models.TextField()
    Clinical_significance = models.TextField()
    Date_last_evaluated = models.TextField()
    Assertion_method = models.TextField()
    Assertion_method_citation = models.TextField()
    Clinical_significance_citations = models.TextField()
    Comment_on_clinical_significance = models.TextField()
    Collection_method = models.TextField()
    Allele_origin = models.TextField()
    ClinVarAccession = models.TextField()
    SAS_Allele_frequency_1000_Genomes = models.TextField()
    EAS_Allele_frequency_1000_Genomes = models.TextField()
    Allele_frequency_1000_Genomes = models.TextField()
    AMR_Allele_frequency_1000Genomes = models.TextField()
    EUR_Allele_frequency_1000Genomes = models.TextField()
    AFR_Allele_frequency_1000Genomes = models.TextField()
    Allele_origin_ClinVar = models.TextField()
    Variant_clinical_significance_ClinVar = models.TextField()
    HGVS_genomic_LOVD = models.TextField()
    Origin_of_variant_LOVD = models.TextField()
    HGVS_protein_LOVD = models.TextField()
    Variant_frequency_LOVD = models.TextField()
    HGVS_cDNA_LOVD = models.TextField()
    Variant_affecting_protein_LOVD = models.TextField()
    Variant_haplotype_LOVD = models.TextField()
    VEP_Gene_ExAC = models.TextField()
    Allele_frequency_ExAC = models.TextField()
    VEP_HGVSc_ExAC = models.TextField()
    VEP_Consequence_ExAC = models.TextField()
    VEP_HGVSp_ExAC = models.TextField()
    Exon_number_exLOVD = models.TextField()
    IARC_class_exLOVD = models.TextField()
    BIC_exLOVD = models.TextField()
    HGVS_cDNA_exLOVD = models.TextField()
    Literature_source_exLOVD = models.TextField()
    HGVS_protein_exLOVD = models.TextField()
    Number_of_family_member_carrying_mutation_BIC = models.TextField()
    HGVS_genomic_BIC = models.TextField()
    Germline_or_Somatic_BIC = models.TextField()
    Mutation_type_BIC = models.TextField()
    BIC_Designation_BIC = models.TextField()
    Literature_citation_BIC = models.TextField()
    Exon_number_BIC = models.TextField()
    Clinical_importance_BIC = models.TextField()
    Clinical_classification_BIC = models.TextField()
    HGVS_protein_BIC = models.TextField()
    HGVS_cDNA_BIC = models.TextField()
    Ethnicity_BIC = models.TextField()
    Patient_nationality_BIC = models.TextField()

    objects = VariantManager()

    class Meta:
        db_table = 'variant'


class Word(models.Model):
    word = models.TextField()

    class Meta:
        db_table = 'words'
