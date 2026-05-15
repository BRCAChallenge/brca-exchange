from django.db import models
from django.db.models import JSONField

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




# ------------------------------------------------------------------------
# --- Variant model
# ------------------------------------------------------------------------

class VariantManager(CopyManager):
    def create_variant(self, row):
        return self.create(**row)


class Variant(models.Model):
    VRS_Digest = models.TextField(primary_key=True)

    # Variant nomenclature
    Gene_Symbol = models.TextField(db_index=True)
    Reference_Sequence = models.TextField()
    HGVS_cDNA = models.TextField(db_index=True)
    BIC_Nomenclature = models.TextField()
    HGVS_Protein = models.TextField()
    Protein_Change = models.TextField()
    CA_ID = models.TextField(null=True, db_index=True)
    VRS = models.JSONField(null=True, blank=True)

    Synonyms = models.TextField()
    title = models.TextField(default='-')
    ensembl_cdna = models.TextField(default='-')
    ensembl_protein = models.TextField(default='-')

    objects = VariantManager()

    class Meta:
        db_table = 'variant'
        managed = False


class Genomic_Coordinates(models.Model):
    VRS_Digest = models.ForeignKey(Variant, on_delete=models.CASCADE, related_name='genomic_coordinates')
    assembly = models.TextField()
    hgvs = models.TextField()
    genome_browser_url = models.TextField()
    chr = models.TextField()
    pos = models.TextField()
    ref = models.TextField()
    alt = models.TextField()

    class Meta:
        db_table = 'genomic_coordinates'
        managed = False


# ------------------------------------------------------------------------
# --- Variant source sub-models
# ------------------------------------------------------------------------

class Variant_in_ClinVar(models.Model):
    """ClinVar data for a variant."""
    VRS_Digest = models.OneToOneField(Variant, primary_key=True, on_delete=models.CASCADE, related_name='clinvar_data')

    Source_URL = models.TextField()

    class Meta:
        db_table = 'variant_clinvar'
        managed = False


class Variant_in_LOVD(models.Model):
    """LOVD (Leiden Open Variation Database) data for a variant."""
    VRS_Digest = models.OneToOneField(Variant, primary_key=True, on_delete=models.CASCADE, related_name='lovd_data')

    Source_URL = models.TextField()
    Variant_haplotype = models.TextField()

    class Meta:
        db_table = 'variant_lovd'
        managed = False


class Variant_in_ExLOVD(models.Model):
    """exLOVD expert-curated BRCA1/2 data for a variant."""
    VRS_Digest = models.OneToOneField(Variant, primary_key=True, on_delete=models.CASCADE, related_name='exlovd_data')

    Source_URL = models.TextField(default='-')
    Exon = models.TextField(default='-')
    DNA_Change = models.TextField(default='-')
    BIC_DNA_Change = models.TextField(default='-')
    Protein_Change = models.TextField(default='-')
    Posterior_P = models.TextField(default='-')
    IARC_Class = models.TextField(default='-')
    DBID = models.TextField(default='-')
    Missense_Analysis_Prior_P = models.TextField(default='-')
    Splicing_Prior_P = models.TextField(default='-')
    Combined_Prior_P = models.TextField(default='-')
    Segregation_LR = models.TextField(default='-')
    Pathology_LR = models.TextField(null=True)
    Co_Occurrence_LR = models.TextField(default='-')
    Case_Control_LR = models.TextField(null=True)
    Product_Of_LRs = models.TextField(default='-')
    Comments = models.TextField(default='-')

    class Meta:
        db_table = 'variant_exlovd'
        managed = False


class Variant_in_GnomAD(models.Model):
    """Per-variant gnomAD anchor row (one per variant across all versions)."""
    VRS_Digest = models.OneToOneField(Variant, primary_key=True, on_delete=models.CASCADE, related_name='gnomad_data')

    Source_URL = models.TextField(default='-')

    class Meta:
        db_table = 'variant_gnomad'
        managed = False


class Report_in_GnomAD(models.Model):
    """GnomAD frequency data — one row per variant per version (v2/v3/v4)."""
    VRS_Digest = models.ForeignKey(Variant_in_GnomAD, on_delete=models.CASCADE, related_name='gnomad_reports')

    version = models.TextField()
    Variant_id = models.TextField(default='-', db_index=True)
    Flags = models.TextField(default='-')
    Allele_count = models.TextField(default='-')
    Allele_number = models.TextField(default='-')
    Allele_frequency = models.TextField(default='-')
    faf95_popmax = models.TextField(default='-')
    faf95_popmax_population = models.TextField(default='-')
    populations = models.JSONField(null=True, blank=True)

    class Meta:
        db_table = 'report_gnomad'
        managed = False
        unique_together = ('VRS_Digest', 'version')


class Variant_in_Other(models.Model):
    """Data from other sources (functional assays, multifactorial, etc.)."""
    VRS_Digest = models.OneToOneField(Variant, primary_key=True, on_delete=models.CASCADE, related_name='other_data')

    data_type = models.TextField(null=True)
    variant_data = models.JSONField(null=True, blank=True)

    class Meta:
        db_table = 'variant_other'
        managed = False



class Variant_in_ENIGMA(models.Model):
    """ENIGMA-specific report data."""
    VRS_Digest = models.ForeignKey(Variant, on_delete=models.CASCADE, related_name='enigma_reports')
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
    Pathogenicity = models.TextField()

    class Meta:
        db_table = 'variant_enigma'
        managed = False


class Report_in_ClinVar(models.Model):
    """ClinVar-specific report data."""
    VRS_Digest = models.ForeignKey(Variant_in_ClinVar, on_delete=models.CASCADE, related_name='clinvar_reports')

    Clinical_Significance = models.TextField()
    Date_Last_Updated = models.TextField()
    DateSignificanceLastEvaluated = models.TextField(default='-')
    Submitter = models.TextField()
    SCV = models.TextField(db_index=True)
    SCV_Version = models.TextField(default='-')
    Allele_Origin = models.TextField()
    Method = models.TextField()
    Description = models.TextField(default="-")
    Summary_Evidence = models.TextField(default="-")
    Review_Status = models.TextField(default="-")
    Condition_Type = models.TextField(default="-")
    Condition_Value = models.TextField(default="-")
    Condition_DB_ID = models.TextField(default="-")

    class Meta:
        db_table = 'report_clinvar'
        managed = False


class Report_in_LOVD(models.Model):
    """LOVD-specific report data."""
    VRS_Digest = models.ForeignKey(Variant_in_LOVD, on_delete=models.CASCADE, related_name='lovd_reports')

    Variant_frequency = models.TextField()
    Individuals = models.TextField()
    Variant_effect = models.TextField()
    Genetic_origin = models.TextField()
    Submitters = models.TextField()
    Functional_analysis_technique = models.TextField(default='-')
    Functional_analysis_result = models.TextField(default='-')
    Created_date = models.TextField(default="-")
    Edited_date = models.TextField(default="-")
    DBID = models.TextField(default="-")
    Remarks = models.TextField(null=True)
    Classification = models.TextField(null=True)
    Submission_ID = models.TextField(default='-', db_index=True)

    class Meta:
        db_table = 'report_lovd'
        managed = False


# ------------------------------------------------------------------------
# --- Other models
# ------------------------------------------------------------------------

class Paper(models.Model):
    PMID = models.TextField(db_index=True)
    Title = models.TextField()
    Author = models.TextField()
    Year = models.TextField()

    class Meta:
        db_table = "data_paper"


class Variant_in_Paper(models.Model):
    VRS_Digest = models.ForeignKey(Variant, on_delete=models.CASCADE)
    Paper = models.ForeignKey(Paper, on_delete=models.CASCADE)

    class Meta:
        db_table = "data_variantpaper"
        unique_together = ('VRS_Digest', 'Paper')



class InSilicoPriors(models.Model):
    """In silico prior probabilities."""

    class Meta:
        db_table = 'data_insilicopriors'
        managed = False


class VariantRepresentation(models.Model):
    """Variant representation in different formats."""

    class Meta:
        db_table = 'data_variantrepresentation'
        managed = False


class AnalysisVEP(models.Model):
    """VEP annotation results for a variant."""
    VRS_Digest    = models.OneToOneField(Variant, primary_key=True, on_delete=models.CASCADE,
                                         related_name='vep_analysis')
    variant_class = models.TextField(null=True)

    class Meta:
        db_table = 'analysis_vep'
        managed = False


class AnalysisBayesDel(models.Model):
    """BayesDel scores for a variant."""
    VRS_Digest               = models.OneToOneField(Variant, primary_key=True, on_delete=models.CASCADE,
                                                     related_name='bayesdel_analysis')
    BayesDel_nsfp33a_noAF   = models.TextField(null=True)

    class Meta:
        db_table = 'analysis_bayesdel'
        managed = False


class AnalysisSpliceAI(models.Model):
    """SpliceAI scores for a variant."""
    VRS_Digest = models.OneToOneField(Variant, primary_key=True, on_delete=models.CASCADE,
                                      related_name='spliceai_analysis')
    DS_AG  = models.TextField(null=True)
    DS_AL  = models.TextField(null=True)
    DS_DG  = models.TextField(null=True)
    DS_DL  = models.TextField(null=True)
    DP_AG  = models.TextField(null=True)
    DP_AL  = models.TextField(null=True)
    DP_DG  = models.TextField(null=True)
    DP_DL  = models.TextField(null=True)
    result = models.TextField(null=True)

    class Meta:
        db_table = 'analysis_spliceai'
        managed = False


class AnalysisProvisionalEvidenceCodes(models.Model):
    """Provisional evidence codes for a variant."""
    VRS_Digest          = models.OneToOneField(Variant, primary_key=True, on_delete=models.CASCADE,
                                               related_name='provisional_evidence_codes')
    popfreq_code        = models.TextField(null=True)
    popfreq_description = models.TextField(null=True)

    class Meta:
        db_table = 'analysis_provisional_evidence_codes'
        managed = False


class AnalysisPriors(models.Model):
    """In silico prior probability fields shown in the UI."""
    VRS_Digest = models.OneToOneField(Variant, primary_key=True, on_delete=models.CASCADE,
                                      related_name='priors_analysis')

    varLoc                        = models.TextField(null=True)

    proteinPrior                  = models.TextField(null=True)

    refDonorPrior                 = models.TextField(null=True)
    refRefDonorMES                = models.TextField(null=True)
    refRefDonorZ                  = models.TextField(null=True)
    altRefDonorMES                = models.TextField(null=True)
    altRefDonorZ                  = models.TextField(null=True)
    refRefDonorSeq                = models.TextField(null=True)
    altRefDonorSeq                = models.TextField(null=True)

    deNovoDonorPrior              = models.TextField(null=True)
    refDeNovoDonorMES             = models.TextField(null=True)
    refDeNovoDonorZ               = models.TextField(null=True)
    altDeNovoDonorMES             = models.TextField(null=True)
    altDeNovoDonorZ               = models.TextField(null=True)
    refDeNovoDonorSeq             = models.TextField(null=True)
    altDeNovoDonorSeq             = models.TextField(null=True)
    deNovoDonorGenomicSplicePos   = models.TextField(null=True)
    deNovoDonorTranscriptSplicePos = models.TextField(null=True)
    closestDonorGenomicSplicePos  = models.TextField(null=True)
    closestDonorTranscriptSplicePos = models.TextField(null=True)
    closestDonorRefMES            = models.TextField(null=True)
    closestDonorRefZ              = models.TextField(null=True)
    closestDonorRefSeq            = models.TextField(null=True)
    closestDonorAltMES            = models.TextField(null=True)
    closestDonorAltZ              = models.TextField(null=True)
    closestDonorAltSeq            = models.TextField(null=True)

    refAccPrior                   = models.TextField(null=True)
    refRefAccMES                  = models.TextField(null=True)
    refRefAccZ                    = models.TextField(null=True)
    altRefAccMES                  = models.TextField(null=True)
    altRefAccZ                    = models.TextField(null=True)
    refRefAccSeq                  = models.TextField(null=True)
    altRefAccSeq                  = models.TextField(null=True)

    class Meta:
        db_table = 'analysis_priors'
        managed = False


class EnigmaDomain(models.Model):
    """ENIGMA Consortium functional domains of potential clinical importance."""
    gene     = models.TextField()
    name     = models.TextField()
    chrom    = models.TextField()
    assembly = models.TextField(default='GRCh38')
    start    = models.IntegerField()
    end      = models.IntegerField()

    class Meta:
        db_table = 'enigma_domain'
        managed = False
