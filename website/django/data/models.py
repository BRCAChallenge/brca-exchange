"""
Normalized database models for BRCA Exchange.

This is a refactored version of the original models.py that normalizes the schema
by splitting large models into smaller, related models grouped by data source.
"""

from django.db import models
from django.db.models import JSONField
from django.contrib.postgres.fields import ArrayField

from postgres_copy import CopyManager


# ------------------------------------------------------------------------
# --- Base models
# ------------------------------------------------------------------------

class DataRelease(models.Model):
    schema = models.TextField()
    archive = models.TextField()
    date = models.DateTimeField()
    notes = models.TextField()
    sources = models.TextField()
    md5sum = models.TextField()
    created = models.DateTimeField(auto_now_add=True, null=True)
    name = models.PositiveIntegerField()

    class Meta:
        db_table = "data_release"
        ordering = ['date']


class ChangeType(models.Model):
    name = models.TextField()

    class Meta:
        db_table = "data_changetype"


class MupitStructure(models.Model):
    name = models.TextField()

    class Meta:
        db_table = "data_mupitstructure"


# ------------------------------------------------------------------------
# --- Core Variant Model (normalized)
# ------------------------------------------------------------------------

class VariantManager(CopyManager):
    def create_variant(self, row):
        return self.create(**row)


class Variant(models.Model):
    """
    Core variant model containing only essential variant information.
    Source-specific data has been moved to related models.
    """
    # Data source info
    Source = models.TextField(db_index=True)
    Variant_in_ENIGMA = models.BooleanField(default=False, db_index=True)
    Variant_in_ClinVar = models.BooleanField(default=False, db_index=True)
    Variant_in_1000_Genomes = models.BooleanField(default=False, db_index=True)
    Variant_in_LOVD = models.BooleanField(default=False, db_index=True)
    Variant_in_Other = models.BooleanField(default=False, db_index=True)
    Variant_in_GnomAD = models.BooleanField(default=False, db_index=True)

    # Variant nomenclature data
    Gene_Symbol = models.TextField(db_index=True)
    Reference_Sequence = models.TextField()
    HGVS_cDNA = models.TextField(db_index=True)
    BIC_Nomenclature = models.TextField()
    HGVS_Protein = models.TextField()
    HGVS_RNA = models.TextField()
    Protein_Change = models.TextField()
    VRS_Digest = models.TextField(null=True, db_index=True)
    VRS_Data = models.JSONField(null=True, blank=True, help_text="GA4GH Variant Representation Standard")
    CA_ID = models.TextField(null=True, db_index=True)

    # Genomic coordinates
    Chr = models.TextField(db_index=True)
    Pos = models.TextField()
    Ref = models.TextField()
    Alt = models.TextField()
    Genomic_Coordinate_hg38 = models.TextField(db_index=True)
    Hg38_Start = models.BigIntegerField(default=1, db_index=True)
    Hg38_End = models.BigIntegerField(default=1)
    Genomic_Coordinate_hg37 = models.TextField()
    Hg37_Start = models.BigIntegerField(default=1, db_index=True)
    Hg37_End = models.BigIntegerField(default=1)
    Genomic_HGVS_38 = models.TextField(null=True)
    Genomic_HGVS_37 = models.TextField(null=True)

    # Allele frequencies (summary)
    Allele_Frequency = models.TextField()
    Allele_Frequency_FAF = models.TextField()

    # Other core fields
    Source_URL = models.TextField()
    Synonyms = models.TextField()
    Pathogenicity_expert = models.TextField(db_index=True)
    Pathogenicity_all = models.TextField()

    # BX IDs for linking
    BX_ID_ENIGMA = models.TextField(default='', db_index=True)
    BX_ID_LOVD = models.TextField(default='', db_index=True)
    BX_ID_ClinVar = models.TextField(default='', db_index=True)
    BX_ID_exLOVD = models.TextField(default='', db_index=True)
    BX_ID_GnomAD = models.TextField(null=True, db_index=True)
    BX_ID_Other = models.TextField(default='', db_index=True)

    # External IDs

    # Data Versioning
    Data_Release = models.ForeignKey(DataRelease, on_delete=models.CASCADE)
    Change_Type = models.ForeignKey(ChangeType, on_delete=models.CASCADE)
    Mupit_Structure = models.ForeignKey(MupitStructure, null=True, on_delete=models.SET_NULL)

    objects = VariantManager()

    class Meta:
        db_table = 'variant'
        indexes = [
            models.Index(fields=['Gene_Symbol', 'Hg38_Start']),
            models.Index(fields=['Chr', 'Hg38_Start']),
        ]


# ------------------------------------------------------------------------
# --- ENIGMA Data
# ------------------------------------------------------------------------

class VariantENIGMA(models.Model):
    """ENIGMA consortium data for a variant."""
    variant = models.OneToOneField(Variant, primary_key=True, on_delete=models.CASCADE, related_name='enigma_data')

    URL_ENIGMA = models.TextField()
    Condition_ID_type_ENIGMA = models.TextField()
    Condition_ID_value_ENIGMA = models.TextField()
    Condition_category_ENIGMA = models.TextField()
    Clinical_significance_ENIGMA = models.TextField(db_index=True)
    Date_last_evaluated_ENIGMA = models.TextField()
    Assertion_method_ENIGMA = models.TextField()
    Assertion_method_citation_ENIGMA = models.TextField()
    Clinical_significance_citations_ENIGMA = models.TextField()
    Comment_on_clinical_significance_ENIGMA = models.TextField()
    Collection_method_ENIGMA = models.TextField()
    Allele_origin_ENIGMA = models.TextField()
    ClinVarAccession_ENIGMA = models.TextField()

    class Meta:
        db_table = 'variant_enigma'


# ------------------------------------------------------------------------
# --- ClinVar Data
# ------------------------------------------------------------------------

class VariantClinVar(models.Model):
    """ClinVar data for a variant."""
    variant = models.OneToOneField(Variant, primary_key=True, on_delete=models.CASCADE, related_name='clinvar_data')

    Clinical_Significance_ClinVar = models.TextField(db_index=True)
    Date_Last_Updated_ClinVar = models.TextField()
    DateSignificanceLastEvaluated_ClinVar = models.TextField(default='-')
    Submitter_ClinVar = models.TextField()
    SCV_ClinVar = models.TextField(db_index=True)
    SCV_Version_ClinVar = models.TextField(default='-')
    Allele_Origin_ClinVar = models.TextField()
    Method_ClinVar = models.TextField()

    class Meta:
        db_table = 'variant_clinvar'


# ------------------------------------------------------------------------
# --- LOVD Data
# ------------------------------------------------------------------------

class VariantLOVD(models.Model):
    """LOVD (Leiden Open Variation Database) data for a variant."""
    variant = models.OneToOneField(Variant, primary_key=True, on_delete=models.CASCADE, related_name='lovd_data')

    Functional_analysis_technique_LOVD = models.TextField()
    Functional_analysis_result_LOVD = models.TextField()
    Variant_frequency_LOVD = models.TextField()
    Variant_haplotype_LOVD = models.TextField()
    HGVS_cDNA_LOVD = models.TextField()
    HGVS_protein_LOVD = models.TextField()
    Individuals_LOVD = models.TextField()
    Variant_effect_LOVD = models.TextField()
    Genetic_origin_LOVD = models.TextField()
    RNA_LOVD = models.TextField()
    Submitters_LOVD = models.TextField()
    Created_date_LOVD = models.TextField(default="-")
    Edited_date_LOVD = models.TextField(default="-")
    Remarks_LOVD = models.TextField(null=True)
    Classification_LOVD = models.TextField(null=True)
    DBID_LOVD = models.TextField(default="-")

    class Meta:
        db_table = 'variant_lovd'


# ------------------------------------------------------------------------
# --- exLOVD Data
# ------------------------------------------------------------------------

class VariantExLOVD(models.Model):
    """exLOVD data for a variant."""
    variant = models.OneToOneField(Variant, primary_key=True, on_delete=models.CASCADE, related_name='exlovd_data')

    IARC_class_exLOVD = models.TextField()
    Sum_family_LR_exLOVD = models.TextField()
    Combined_prior_probablility_exLOVD = models.TextField()
    Literature_source_exLOVD = models.TextField()
    Co_occurrence_LR_exLOVD = models.TextField()
    Pathology_LR_exLOVD = models.TextField(null=True)
    Case_control_LR_exLOVD = models.TextField(null=True)
    Posterior_probability_exLOVD = models.TextField()
    Missense_analysis_prior_probability_exLOVD = models.TextField()
    Segregation_LR_exLOVD = models.TextField()

    class Meta:
        db_table = 'variant_exlovd'


# ------------------------------------------------------------------------
# --- GnomAD Data (v2, v3 or v4, depending on the gnomAD_version field)
# ------------------------------------------------------------------------

class VariantGnomAD(models.Model):
    """
    GnomAD (Genome Aggregation Database) data for a variant.
    Due to the large number of population-specific fields, consider
    further normalization if needed.
    """
    variant = models.OneToOneField(Variant, primary_key=True, on_delete=models.CASCADE, related_name='gnomad_data')

    gnomAD_version = models.TextField(null=True)
    HGVS_cDNA_GnomAD = models.TextField(null=True)
    HGVS_GnomAD = models.TextField(null=True)
    HGVS_protein_GnomAD = models.TextField(null=True)
    Flags_GnomAD = models.TextField(null=True)
    Consequence_GnomAD = models.TextField(null=True)
    Variant_id_GnomAD = models.TextField(null=True, db_index=True)
    faf95_popmax_genome_GnomAD = models.TextField(null=True)
    faf95_popmax_population_genome_GnomAD = models.TextField(null=True)
    faf95_popmax_exome_GnomAD = models.TextField(null=True)
    faf95_popmax_population_exome_GnomAD = models.TextField(null=True)
    faf95_popmax_joint_GnomAD = models.TextField(null=True)
    faf95_popmax_population_joint_GnomAD = models.TextField(null=True)

    # Overall genome frequencies
    Allele_count_genome_GnomAD = models.TextField(null=True)
    Allele_number_genome_GnomAD = models.TextField(null=True)
    Allele_frequency_genome_GnomAD = models.TextField(null=True)

    # Overall exome frequencies
    Allele_number_exome_GnomAD = models.TextField(null=True)
    Allele_count_exome_GnomAD = models.TextField(null=True)
    Allele_frequency_exome_GnomAD = models.TextField(null=True)

    # Overall joint frequencies
    Allele_number_joint_GnomAD = models.TextField(null=True)
    Allele_count_joint_GnomAD = models.TextField(null=True)
    Allele_frequency_joint_GnomAD = models.TextField(null=True)

    # Population-specific data stored as JSON for flexibility
    # This approach reduces the number of individual fields while maintaining data integrity
    genome_populations = models.JSONField(null=True, blank=True, help_text="Genome frequency data by population")
    exome_populations = models.JSONField(null=True, blank=True, help_text="Exome frequency data by population")
    joint_populations = models.JSONField(null=True, blank=True, help_text="Joint frequency data by population")

    class Meta:
        db_table = 'variant_gnomad'



# ------------------------------------------------------------------------
# --- Data from other sources
# ------------------------------------------------------------------------

class VariantOther(models.Model):
    """Functional assay data for a variant (Findlay, ENIGMA BRCA12, etc.)."""
    variant = models.OneToOneField(Variant, primary_key=True, on_delete=models.CASCADE, related_name='functional_assays')

    # The data type, such as Functional Assay, Multifactorial, etc
    data_type = models.TextField(null=True) 
    # Different forms of variant data - stored as JSON due to complexity
    variant_data = models.JSONField(null=True, blank=True, help_text="ENIGMA BRCA1/2 functional assay results")

    class Meta:
        db_table = 'variant_other'


# ------------------------------------------------------------------------
# --- In Silico Predictions
# ------------------------------------------------------------------------

class VariantInSilicoPredictions(models.Model):
    """In silico prediction scores for a variant."""
    variant = models.OneToOneField(Variant, primary_key=True, on_delete=models.CASCADE, related_name='insilico_predictions')

    # SpliceAI
    DS_AG_spliceAI = models.TextField(null=True)
    DS_AL_spliceAI = models.TextField(null=True)
    DS_DG_spliceAI = models.TextField(null=True)
    DS_DL_spliceAI = models.TextField(null=True)
    DP_AG_spliceAI = models.TextField(null=True)
    DP_AL_spliceAI = models.TextField(null=True)
    DP_DG_spliceAI = models.TextField(null=True)
    DP_DL_spliceAI = models.TextField(null=True)
    result_spliceai = models.TextField(null=True)

    # BayesDel
    BayesDel_nsfp33a_noAF = models.TextField(null=True)

    # Provisional Evidence
    Provisional_Evidence_Code_Popfreq = models.TextField(default='')
    Provisional_Evidence_Description_Popfreq = models.TextField(default='')
    Provisional_Evidence_Code_Bioinfo = models.TextField(default='')
    Provisional_Evidence_Description_Bioinfo = models.TextField(default='')

    class Meta:
        db_table = 'variant_insilico_predictions'


# ------------------------------------------------------------------------
# --- Diff Models
# ------------------------------------------------------------------------

class VariantDiff(models.Model):
    variant = models.OneToOneField(Variant, primary_key=True, on_delete=models.CASCADE)
    # Django 3.1+ JSONField natively supports PostgreSQL JSON type
    diff = JSONField()

    class Meta:
        db_table = 'data_variantdiff'


class ReportDiff(models.Model):
    report = models.OneToOneField('Report', primary_key=True, on_delete=models.CASCADE)
    # Django 3.1+ JSONField natively supports PostgreSQL JSON type
    report_diff = JSONField()

    class Meta:
        db_table = 'data_reportdiff'


# ------------------------------------------------------------------------
# --- Report Model (similarly normalized)
# ------------------------------------------------------------------------

class ReportManager(models.Manager):
    def create_report(self, row):
        return self.create(**row)


class Report(models.Model):
    """
    Core report model. Source-specific fields moved to related models.
    """
    Variant = models.ForeignKey(Variant, on_delete=models.CASCADE)
    Source = models.TextField(db_index=True)

    # Core genomic information
    Gene_symbol_ENIGMA = models.TextField()
    Genomic_Coordinate = models.TextField()
    Chr = models.TextField()
    Pos = models.TextField()
    Ref = models.TextField()
    Alt = models.TextField()
    Reference_sequence_ENIGMA = models.TextField()
    HGVS_cDNA_ENIGMA = models.TextField()
    BIC_Nomenclature_ENIGMA = models.TextField()
    Abbrev_AA_change_ENIGMA = models.TextField()

    objects = ReportManager()

    class Meta:
        db_table = 'data_report'


class ReportENIGMA(models.Model):
    """ENIGMA-specific report data."""
    report = models.OneToOneField(Report, primary_key=True, on_delete=models.CASCADE, related_name='enigma_report')

    URL_ENIGMA = models.TextField()
    Condition_ID_type_ENIGMA = models.TextField()
    Condition_ID_value_ENIGMA = models.TextField()
    Condition_category_ENIGMA = models.TextField()
    Clinical_significance_ENIGMA = models.TextField()
    Date_last_evaluated_ENIGMA = models.TextField()
    Assertion_method_ENIGMA = models.TextField()
    Assertion_method_citation_ENIGMA = models.TextField()
    Clinical_significance_citations_ENIGMA = models.TextField()
    Comment_on_clinical_significance_ENIGMA = models.TextField()
    Collection_method_ENIGMA = models.TextField()
    Allele_origin_ENIGMA = models.TextField()
    ClinVarAccession_ENIGMA = models.TextField()
    HGVS_protein_ENIGMA = models.TextField()
    BX_ID_ENIGMA = models.TextField()

    class Meta:
        db_table = 'report_enigma'


class ReportClinVar(models.Model):
    """ClinVar-specific report data."""
    report = models.OneToOneField(Report, primary_key=True, on_delete=models.CASCADE, related_name='clinvar_report')

    Clinical_Significance_ClinVar = models.TextField()
    Date_Last_Updated_ClinVar = models.TextField()
    DateSignificanceLastEvaluated_ClinVar = models.TextField(default='-')
    BX_ID_ClinVar = models.TextField()
    HGVS_ClinVar = models.TextField()
    Submitter_ClinVar = models.TextField()
    Protein_ClinVar = models.TextField()
    SCV_ClinVar = models.TextField(db_index=True)
    SCV_Version_ClinVar = models.TextField(default='-')
    Allele_Origin_ClinVar = models.TextField()
    Method_ClinVar = models.TextField()
    Description_ClinVar = models.TextField(default="-")
    Summary_Evidence_ClinVar = models.TextField(default="-")
    Review_Status_ClinVar = models.TextField(default="-")
    Condition_Type_ClinVar = models.TextField(default="-")
    Condition_Value_ClinVar = models.TextField(default="-")
    Condition_DB_ID_ClinVar = models.TextField(default="-")
    Synonyms_ClinVar = models.TextField(default="-")

    class Meta:
        db_table = 'report_clinvar'


class ReportLOVD(models.Model):
    """LOVD-specific report data."""
    report = models.OneToOneField(Report, primary_key=True, on_delete=models.CASCADE, related_name='lovd_report')

    BX_ID_LOVD = models.TextField()
    Variant_frequency_LOVD = models.TextField()
    HGVS_cDNA_LOVD = models.TextField()
    HGVS_protein_LOVD = models.TextField()
    Individuals_LOVD = models.TextField()
    Variant_effect_LOVD = models.TextField()
    Genetic_origin_LOVD = models.TextField()
    RNA_LOVD = models.TextField()
    Submitters_LOVD = models.TextField()
    Functional_analysis_technique_LOVD = models.TextField(default='-')
    Functional_analysis_result_LOVD = models.TextField(default='-')
    Created_date_LOVD = models.TextField(default="-")
    Edited_date_LOVD = models.TextField(default="-")
    DBID_LOVD = models.TextField(default="-")
    Remarks_LOVD = models.TextField(null=True)
    Classification_LOVD = models.TextField(null=True)
    Submission_ID_LOVD = models.TextField(default='-', db_index=True)

    class Meta:
        db_table = 'report_lovd'


# Additional Report models for other data sources can be added similarly...


# ------------------------------------------------------------------------
# --- Other Models
# ------------------------------------------------------------------------

class Paper(models.Model):
    Variant_in_gnomAD = models.BooleanField(default=False)
    Variant_in_ExAC = models.BooleanField(default=False)
    PMID = models.TextField(db_index=True)
    Title = models.TextField()
    Author = models.TextField()
    Year = models.TextField()

    class Meta:
        db_table = "data_paper"


class VariantPaper(models.Model):
    Variant = models.ForeignKey(Variant, on_delete=models.CASCADE)
    Paper = models.ForeignKey(Paper, on_delete=models.CASCADE)

    class Meta:
        db_table = "data_variantpaper"
        unique_together = ('Variant', 'Paper')


class CurrentVariant(models.Model):
    """
    Represents the current version of a variant.
    This model can also be normalized similar to Variant if needed.
    """
    # For brevity, keeping as-is. Can be normalized similarly to Variant.
    pass


class InSilicoPriors(models.Model):
    """In silico prior probabilities."""
    pass


class VariantRepresentation(models.Model):
    """Variant representation in different formats."""
    pass
