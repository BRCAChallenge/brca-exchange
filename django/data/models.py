from django.db import models
from django.db.models import JSONField
from django.contrib.postgres.fields import ArrayField


class LegacyJSONField(JSONField):
    """
    Referenced by migrations 0009_variantdiff.py and 0027_reportdiff.py, whose
    VariantDiff/ReportDiff models were dropped from the schema without a
    corresponding migration. Kept only so the migration graph still loads;
    not used by any current model.
    """
    def from_db_value(self, value, expression, connection):
        # psycopg3 already deserializes JSONB to Python objects, so handle both
        # str and already-deserialized values.
        if value is None:
            return value
        if isinstance(value, str):
            return super().from_db_value(value, expression, connection)
        return value


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

    # Static, stored release-over-release stats. The old schema computed these
    # live from a per-release Variant history with a Change_Type FK, neither of
    # which exist in the finalized schema (Variant is current-state-only), so
    # these are plain columns instead. variants_classified is a placeholder
    # (always 0) until classification-change tracking is designed.
    variants_added = models.PositiveIntegerField(default=0)
    variants_modified = models.PositiveIntegerField(default=0)
    variants_deleted = models.PositiveIntegerField(default=0)
    variants_classified = models.PositiveIntegerField(default=0)

    class Meta:
        db_table = "data_release"
        ordering = ['date']




# ------------------------------------------------------------------------
# --- Variant model
# ------------------------------------------------------------------------

class VariantManager(models.Manager):
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


class Genomic_Coordinates(models.Model):
    VRS_Digest = models.ForeignKey(Variant, on_delete=models.CASCADE, related_name='genomic_coordinates', db_column='VRS_Digest')
    assembly = models.TextField()
    hgvs = models.TextField()
    chr = models.TextField()
    pos = models.TextField()
    end_pos = models.TextField()
    ref = models.TextField()
    alt = models.TextField()

    class Meta:
        db_table = 'variant_genomic_coordinates'
        unique_together = (('VRS_Digest', 'assembly'),)


# ------------------------------------------------------------------------
# --- Variant source sub-models
# ------------------------------------------------------------------------

class Variant_in_ClinVar(models.Model):
    """ClinVar data for a variant."""
    VRS_Digest = models.OneToOneField(Variant, primary_key=True, on_delete=models.CASCADE, related_name='clinvar_data', db_column='VRS_Digest')

    Source_URL = models.TextField()

    # VCV-level aggregate germline classification
    Variation_ID = models.TextField(default='-')
    VCV_Accession = models.TextField(default='-')
    VCV_Version = models.TextField(default='-')
    Clinical_Significance = models.TextField(default='-')
    Review_Status = models.TextField(default='-')
    Date_Last_Evaluated = models.TextField(default='-')
    Number_Of_Submissions = models.TextField(default='-')
    Number_Of_Submitters = models.TextField(default='-')
    Most_Recent_Submission = models.TextField(default='-')
    Explanation = models.TextField(default='-')
    Conditions = models.TextField(default='-')
    Condition_DB_IDs = models.TextField(default='-')
    Date_Last_Updated = models.TextField(default='-')

    class Meta:
        db_table = 'variant_clinvar'


class Variant_in_LOVD(models.Model):
    """LOVD (Leiden Open Variation Database) data for a variant."""
    VRS_Digest = models.OneToOneField(Variant, primary_key=True, on_delete=models.CASCADE, related_name='lovd_data', db_column='VRS_Digest')

    Source_URL = models.TextField()
    Variant_haplotype = models.TextField()

    class Meta:
        db_table = 'variant_lovd'


class Variant_in_ExLOVD(models.Model):
    """exLOVD expert-curated BRCA1/2 data for a variant."""
    VRS_Digest = models.OneToOneField(Variant, primary_key=True, on_delete=models.CASCADE, related_name='exlovd_data', db_column='VRS_Digest')

    Source_URL = models.TextField(null=True)
    Exon = models.TextField(null=True)
    DNA_Change = models.TextField(null=True)
    BIC_DNA_Change = models.TextField(null=True)
    Protein_Change = models.TextField(null=True)
    DBID = models.TextField(null=True)
    Posterior_P = models.TextField(default='-', null=True)
    IARC_Class = models.TextField(default='-', null=True)
    Missense_Analysis_Prior_P = models.TextField(default='-', null=True)
    Combined_Prior_P = models.TextField(default='-', null=True)
    Segregation_LR = models.TextField(default='-', null=True)
    Splicing_Prior_P = models.TextField(null=True)
    Pathology_LR = models.TextField(null=True)
    Co_Occurrence_LR = models.TextField(default='-', null=True)
    Case_Control_LR = models.TextField(null=True)
    Product_Of_LRs = models.TextField(null=True)
    Comments = models.TextField(default='-', null=True)

    class Meta:
        db_table = 'variant_exlovd'


class Variant_in_GnomAD(models.Model):
    """Per-variant gnomAD anchor row (one per variant across all versions)."""
    VRS_Digest = models.OneToOneField(Variant, primary_key=True, on_delete=models.CASCADE, related_name='gnomad_data', db_column='VRS_Digest')

    Source_URL = models.TextField(default='-', null=True)

    class Meta:
        db_table = 'variant_gnomad'


class Report_in_GnomAD(models.Model):
    """GnomAD frequency data — one row per variant per version and data type."""
    VRS_Digest = models.ForeignKey(Variant_in_GnomAD, on_delete=models.CASCADE, related_name='gnomad_reports', db_column='VRS_Digest')

    version = models.TextField(null=True)
    data_type = models.TextField(default='-')    # 'joint', 'genome', or 'exome'
    Variant_id = models.TextField(default='-', null=True, db_index=True)
    Flags = models.TextField(default='-', null=True)
    coverage = models.TextField(default='-')
    Allele_count = models.TextField(default='-', null=True)
    Allele_number = models.TextField(default='-', null=True)
    Allele_frequency = models.TextField(default='-', null=True)
    faf95_popmax = models.TextField(default='-', null=True)
    faf95_popmax_population = models.TextField(default='-', null=True)
    populations = models.JSONField(null=True, blank=True)

    class Meta:
        db_table = 'report_gnomad'
        unique_together = ('VRS_Digest', 'version', 'data_type')


class Variant_in_Other(models.Model):
    """Data from other sources (functional assays, multifactorial, etc.)."""
    VRS_Digest = models.OneToOneField(Variant, primary_key=True, on_delete=models.CASCADE, related_name='other_data', db_column='VRS_Digest')

    data_type = models.TextField(null=True)
    variant_data = models.JSONField(null=True, blank=True)

    class Meta:
        db_table = 'variant_other'



class Variant_in_ENIGMA(models.Model):
    """ENIGMA-specific report data."""
    VRS_Digest = models.ForeignKey(Variant, on_delete=models.CASCADE, related_name='enigma_reports', db_column='VRS_Digest')
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


class Report_in_ClinVar(models.Model):
    """ClinVar-specific report data."""
    VRS_Digest = models.ForeignKey(Variant_in_ClinVar, on_delete=models.CASCADE, related_name='clinvar_reports', db_column='VRS_Digest')

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


class Report_in_LOVD(models.Model):
    """LOVD-specific report data."""
    VRS_Digest = models.ForeignKey(Variant_in_LOVD, on_delete=models.CASCADE, related_name='lovd_reports', db_column='VRS_Digest')

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


# ------------------------------------------------------------------------
# --- Other models
# ------------------------------------------------------------------------

class Paper(models.Model):
    PMID = models.TextField(db_index=True)
    Title = models.TextField()
    Author = models.TextField()
    Year = models.TextField()
    Journal = models.TextField()
    DOI = models.TextField()

    class Meta:
        db_table = "data_paper"


class Variant_in_Paper(models.Model):
    VRS_Digest = models.ForeignKey(Variant, on_delete=models.CASCADE, db_column='VRS_Digest')
    Paper = models.ForeignKey(Paper, on_delete=models.CASCADE)
    mentions = ArrayField(models.TextField())
    variant_mentioned_as = ArrayField(models.TextField())

    class Meta:
        db_table = "data_variantpaper"
        unique_together = ('VRS_Digest', 'Paper')



# VEP consequence_terms that indicate a variant introduces a premature termination codon.
# Splice-altering terms are deliberately excluded: VEP's consequence call for those doesn't
# predict the resulting spliced transcript, so there's no reliable basis for a PTC call here.
PTC_CONSEQUENCE_TERMS = {'stop_gained', 'frameshift_variant'}


class AnalysisVEPQuerySet(models.QuerySet):
    def introduces_ptc(self):
        q = models.Q()
        for term in PTC_CONSEQUENCE_TERMS:
            q |= models.Q(consequences__regex=rf'(^|,){term}(,|$)')
        return self.filter(q)


class AnalysisVEP(models.Model):
    """VEP annotation results for a variant."""
    objects       = AnalysisVEPQuerySet.as_manager()
    VRS_Digest    = models.OneToOneField(Variant, primary_key=True, on_delete=models.CASCADE,
                                         related_name='vep_analysis', db_column='VRS_Digest')
    variant_class = models.TextField(null=True)
    variant_type  = models.TextField(null=True)
    # VEP consequence_terms for the Variant.Reference_Sequence RefSeq transcript, comma-separated
    consequences  = models.TextField(null=True)
    # VEP HGVSp notation for the Variant.Reference_Sequence RefSeq transcript, e.g.
    # "NP_009225.1:p.Gln356GlufsTer9"
    hgvsp         = models.TextField(null=True)
    # GRCh38 1-based genomic position of the first base of the premature stop codon
    # described by hgvsp; only set when consequences intersects PTC_CONSEQUENCE_TERMS
    # and the position could be resolved. Chromosome is implied by Gene_Symbol
    # (BRCA1 -> chr17, BRCA2 -> chr13).
    ptc_genomic_pos = models.IntegerField(null=True)
    VA_Spec       = models.JSONField(null=True, blank=True)

    class Meta:
        db_table = 'analysis_vep'

    @property
    def introduces_ptc(self):
        terms = set((self.consequences or '').split(','))
        return bool(terms & PTC_CONSEQUENCE_TERMS)


class AnalysisBayesDel(models.Model):
    """BayesDel scores for a variant."""
    VRS_Digest               = models.OneToOneField(Variant, primary_key=True, on_delete=models.CASCADE,
                                                     related_name='bayesdel_analysis', db_column='VRS_Digest')
    BayesDel_nsfp33a_noAF   = models.TextField(null=True)
    VA_Spec                 = models.JSONField(null=True, blank=True)

    class Meta:
        db_table = 'analysis_bayesdel'


class AnalysisSpliceAI(models.Model):
    """SpliceAI scores for a variant."""
    VRS_Digest = models.OneToOneField(Variant, primary_key=True, on_delete=models.CASCADE,
                                      related_name='spliceai_analysis', db_column='VRS_Digest')
    DS_AG  = models.TextField(null=True)
    DS_AL  = models.TextField(null=True)
    DS_DG  = models.TextField(null=True)
    DS_DL  = models.TextField(null=True)
    DP_AG  = models.TextField(null=True)
    DP_AL  = models.TextField(null=True)
    DP_DG  = models.TextField(null=True)
    DP_DL  = models.TextField(null=True)
    result = models.TextField(null=True)
    VA_Spec = models.JSONField(null=True, blank=True)

    class Meta:
        db_table = 'analysis_spliceai'


class AnalysisProvisionalEvidenceCodes(models.Model):
    """Provisional evidence codes for a variant. Multiple rows per variant,
    one per method_name."""
    VRS_Digest          = models.ForeignKey(Variant, on_delete=models.CASCADE,
                                             related_name='provisional_evidence_codes', db_column='VRS_Digest')
    popfreq_code        = models.TextField(null=True)
    popfreq_description = models.TextField(null=True)
    method_name         = models.TextField()
    gnomad_version      = models.TextField(null=True)
    gnomad_data_type    = models.TextField(null=True)
    VA_Spec             = models.JSONField(null=True, blank=True)

    class Meta:
        db_table = 'analysis_provisional_evidence_codes'
        unique_together = (('VRS_Digest', 'method_name'),)


class AnalysisPriors(models.Model):
    """In silico prior probability fields shown in the UI."""
    VRS_Digest = models.OneToOneField(Variant, primary_key=True, on_delete=models.CASCADE,
                                      related_name='priors_analysis', db_column='VRS_Digest')

    varLoc                        = models.TextField(null=True)
    applicablePrior               = models.TextField(null=True)

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

    VA_Spec                       = models.JSONField(null=True, blank=True)

    class Meta:
        db_table = 'analysis_priors'


class EnigmaDomain(models.Model):
    """ENIGMA Consortium functional domains of potential clinical importance."""
    gene     = models.TextField()
    name     = models.TextField()
    chrom    = models.TextField()
    assembly = models.TextField(default='GRCh38')
    start    = models.IntegerField()
    end      = models.IntegerField()

    class Meta:
        db_table = 'data_enigma_domain'
