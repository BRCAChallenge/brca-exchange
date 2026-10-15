"""
Maps the old flat schema's column names (still used by the frontend's
column= query param, e.g. website/js/VariantTable.js) onto the finalized,
normalized schema (django/data/models.py).

Three kinds of mapping:
  - a plain string: a values()-safe field path, used as-is.
  - ('subquery', build_fn): build_fn(outer_ref) returns a Subquery expression
    to annotate onto the queryset (for one-to-many relations flattened to a
    single "winning" value per variant).
  - Columns not in COLUMN_MAP have no home in the finalized schema (either
    genuinely dropped, or not yet resolved) and should be rejected.
"""

from django.db.models import F, OuterRef, Subquery
from django.db.models.fields.json import KeyTextTransform

from .models import Report_in_ClinVar, Report_in_LOVD, Report_in_GnomAD, Genomic_Coordinates


def _genomic_coordinate_subquery(assembly, field='hgvs'):
    def build(outer_ref):
        return Subquery(
            Genomic_Coordinates.objects
            .filter(VRS_Digest=OuterRef(outer_ref), assembly=assembly)
            .values(field)[:1]
        )
    return ('subquery', build)


def _clinvar_subquery(field):
    def build(outer_ref):
        return Subquery(
            Report_in_ClinVar.objects
            .filter(VRS_Digest__VRS_Digest=OuterRef(outer_ref))
            .order_by('-Date_Last_Updated')
            .values(field)[:1]
        )
    return ('subquery', build)


def _lovd_subquery(field):
    def build(outer_ref):
        return Subquery(
            Report_in_LOVD.objects
            .filter(VRS_Digest__VRS_Digest=OuterRef(outer_ref))
            .order_by('-Edited_date')
            .values(field)[:1]
        )
    return ('subquery', build)


def _gnomad_subquery(version, data_type, field=None, population=None):
    """field: a top-level Report_in_GnomAD field. population: a lowercase
    population code, extracted from the populations JSONField's "af" key."""
    def build(outer_ref):
        qs = Report_in_GnomAD.objects.filter(
            VRS_Digest__VRS_Digest=OuterRef(outer_ref), version=version, data_type=data_type,
        )
        if population:
            qs = qs.annotate(_val=KeyTextTransform('af', KeyTextTransform(population, 'populations')))
            qs = qs.values('_val')
        else:
            qs = qs.values(field)
        return Subquery(qs[:1])
    return ('subquery', build)


_GNOMAD_POPULATIONS = ['AFR', 'AMR', 'ASJ', 'EAS', 'FIN', 'NFE', 'OTH', 'SAS', 'MID', 'AMI']
_GNOMAD_V41_POPULATIONS = ['AFR', 'AMI', 'AMR', 'ASJ', 'EAS', 'FIN', 'MID', 'NFE', 'REMAINING', 'SAS']

def _build_gnomad_columns():
    cols = {}
    # current _GnomAD family -> v2/exome (confirmed always empty right now; that's expected)
    cols['Allele_frequency_exome_GnomAD'] = _gnomad_subquery('v2', 'exome', field='Allele_frequency')
    cols['faf95_popmax_exome_GnomAD'] = _gnomad_subquery('v2', 'exome', field='faf95_popmax')
    for pop in _GNOMAD_POPULATIONS[:8]:  # v2 doesn't have MID/AMI
        cols[f'Allele_frequency_exome_{pop}_GnomAD'] = _gnomad_subquery('v2', 'exome', population=pop.lower())
        cols[f'Allele_frequency_genome_{pop}_GnomAD'] = _gnomad_subquery('v2', 'exome', population=pop.lower())
    cols['Allele_frequency_genome_GnomAD'] = _gnomad_subquery('v2', 'exome', field='Allele_frequency')

    # _GnomADv3 family -> v3/genome
    cols['Allele_frequency_genome_GnomADv3'] = _gnomad_subquery('v3', 'genome', field='Allele_frequency')
    cols['faf95_popmax_genome_GnomADv3'] = _gnomad_subquery('v3', 'genome', field='faf95_popmax')
    for pop in _GNOMAD_POPULATIONS:
        cols[f'Allele_frequency_genome_{pop}_GnomADv3'] = _gnomad_subquery('v3', 'genome', population=pop.lower())
    return cols

    # _GnomADv41 families -> v4.1/joint and v4.1/exome
    for dt in ('joint', 'exome'):
        cols[f'Allele_frequency_{dt}_GnomADv41'] = _gnomad_subquery('v4.1', dt, field='Allele_frequency')
        cols[f'faf95_popmax_{dt}_GnomADv41'] = _gnomad_subquery('v4.1', dt, field='faf95_popmax')
        cols[f'faf95_popmax_population_{dt}_GnomADv41'] = _gnomad_subquery('v4.1', dt, field='faf95_popmax_population')
        for pop in _GNOMAD_V41_POPULATIONS:
            if dt == 'exome' and pop == 'AMI':  # v4.1 exomes have no Amish
                continue
            cols[f'Allele_frequency_{dt}_{pop}_GnomADv41'] = _gnomad_subquery('v4.1', dt, population=pop.lower())
    return cols


# Old suffix -> (related_name on Variant, model field name). Model field
# defaults to the same name as the frontend's (suffix-stripped) prop unless
# given explicitly in the RENAMES dict below.
_SINGLE_VALUED = {
    '_exLOVD': 'exlovd_data',
    '_spliceAI': 'spliceai_analysis',
    '_spliceai': 'spliceai_analysis',
}

_RENAMES = {
    # exLOVD
    'Posterior_probability_exLOVD': 'exlovd_data__Posterior_P',
    'Combined_prior_probablility_exLOVD': 'exlovd_data__Combined_Prior_P',
    'Missense_analysis_prior_probability_exLOVD': 'exlovd_data__Missense_Analysis_Prior_P',
    'Literature_source_exLOVD': 'exlovd_data__Comments',
    'Sum_family_LR_exLOVD': 'exlovd_data__Segregation_LR',
    'Case_control_LR_exLOVD': 'exlovd_data__Case_Control_LR',
    'Co_occurrence_LR_exLOVD': 'exlovd_data__Co_Occurrence_LR',
    'Segregation_LR_exLOVD': 'exlovd_data__Segregation_LR',
    'Pathology_LR_exLOVD': 'exlovd_data__Pathology_LR',
    # base field renames
    'VR_ID': 'VRS_Digest',
    'HGVS_Protein_ID': 'HGVS_Protein',
    # ENIGMA (Variant_in_ENIGMA, currently 1:1 in practice)
    'Allele_origin_ENIGMA': 'enigma_reports__Allele_origin',
    'Assertion_method_ENIGMA': 'enigma_reports__Assertion_method',
    'Assertion_method_citation_ENIGMA': 'enigma_reports__Assertion_method_citation',
    'ClinVarAccession_ENIGMA': 'enigma_reports__ClinVarAccession',
    'Clinical_significance_ENIGMA': 'enigma_reports__Clinical_significance',
    'Clinical_significance_citations_ENIGMA': 'enigma_reports__Clinical_significance_citations',
    'Collection_method_ENIGMA': 'enigma_reports__Collection_method',
    'Comment_on_clinical_significance_ENIGMA': 'enigma_reports__Comment_on_clinical_significance',
    'Condition_ID_type_ENIGMA': 'enigma_reports__Condition_ID_type',
    'Condition_ID_value_ENIGMA': 'enigma_reports__Condition_ID_value',
    'Condition_category_ENIGMA': 'enigma_reports__Condition_category',
    'Date_last_evaluated_ENIGMA': 'enigma_reports__Date_last_evaluated',
    'Pathogenicity_expert': 'enigma_reports__Pathogenicity',
    # LOVD anchor table (not Report_in_LOVD)
    'Variant_haplotype_LOVD': 'lovd_data__Variant_haplotype',
    # SpliceAI / BayesDel / Priors
    'result_spliceai': 'spliceai_analysis__result',
    'BayesDel_nsfp33a_noAF': 'bayesdel_analysis__BayesDel_nsfp33a_noAF',
    'applicablePrior': 'priors_analysis__applicablePrior',
    'proteinPrior': 'priors_analysis__proteinPrior',
    'refDonorPrior': 'priors_analysis__refDonorPrior',
    'deNovoDonorPrior': 'priors_analysis__deNovoDonorPrior',
    'refAccPrior': 'priors_analysis__refAccPrior',
}

for suf in ['DP_AG', 'DP_AL', 'DP_DG', 'DP_DL', 'DS_AG', 'DS_AL', 'DS_DG', 'DS_DL']:
    _RENAMES[f'{suf}_spliceAI'] = f'spliceai_analysis__{suf}'

# Genomic coordinates: one row per (variant, assembly), needs the same fixed-
# filter Subquery pattern as GnomAD, not a plain direct traversal.
_RENAMES['Genomic_Coordinate_hg37'] = _genomic_coordinate_subquery('GRCh37')
_RENAMES['Genomic_Coordinate_hg38'] = _genomic_coordinate_subquery('GRCh38')

# GRCh38 ref/alt/position, used by the splice-site diagram (Splicing.js).
# Pos and Hg38_Start are the same value under two old names (the old
# frontend used both inconsistently) - end_pos is currently null for all
# locally-restored data (not yet backfilled by the pipeline), so the
# diagram's overlap detection will be degraded until that lands, but it
# won't crash.
_RENAMES['Ref'] = _genomic_coordinate_subquery('GRCh38', field='ref')
_RENAMES['Alt'] = _genomic_coordinate_subquery('GRCh38', field='alt')
_RENAMES['Pos'] = _genomic_coordinate_subquery('GRCh38', field='pos')
_RENAMES['Hg38_Start'] = _genomic_coordinate_subquery('GRCh38', field='pos')
_RENAMES['Hg38_End'] = _genomic_coordinate_subquery('GRCh38', field='end_pos')

# GRCh37 chromosome/position, used by the "Beacons" external link.
_RENAMES['Chr'] = _genomic_coordinate_subquery('GRCh37', field='chr')
_RENAMES['Hg37_Start'] = _genomic_coordinate_subquery('GRCh37', field='pos')

# ENIGMA functional-assay results: single JSONField on Variant_in_Other,
# data_type='functional_assay_results'. A few keys were renamed from the old
# column name; everything else is the same name with the suffix stripped.
_FUNCTIONAL_ASSAY_KEY_RENAMES = {
    'HGVS_Nucleotide': 'HGVS_Nucleotide_Variant',
    'Functional_Enrichment_Findlay': 'Function_Score_Findlay',
    'Cell_Survival_Biwas': 'HAT_(Cell_Survival)_Biwas',
    'Drug_Sensitivity_Biwas': 'DS_(Drug_Sensitivity)_Biwas',
    'HAT_DS_Score_Biwas': 'HAT+DS_Score_Biwas',
}
_FUNCTIONAL_ASSAY_COLUMNS = [
    'CBDCA_fClass_Ikegami', 'Cell_Survival_Biwas', 'Chromosomal_Variant', 'Cisplatin_Bouwman2',
    'Cisplatin_Mesman', 'Complementation_Mesman', 'Control_Group_Petitalot', 'DRGFP_Bouwman2',
    'Drug_Sensitivity_Biwas', 'Functional_Enrichment_Findlay', 'HAT_DS_Score_Biwas', 'HDR_Mesman',
    'HDR_Richardson', 'HGVS_Nucleotide', 'Niraparif_fClass_Ikegami', 'Olaparib_Bouwman2',
    'Olaparib_fClass_Ikegami', 'RNA_Class_Findlay', 'RNA_Score_Findlay', 'Result_Biwas',
    'Result_Bouwman1', 'Result_Bouwman2', 'Result_Fernandes', 'Result_Findlay', 'Result_Ikegami',
    'Result_Mesman', 'Result_Petitalot', 'Result_Richardson', 'Result_Starita',
    'Rucaparib_fClass_Ikegami', 'Selection_Bouwman1',
]
for base in _FUNCTIONAL_ASSAY_COLUMNS:
    json_key = _FUNCTIONAL_ASSAY_KEY_RENAMES.get(base, base)
    _RENAMES[f'{base}_ENIGMA_BRCA12_Functional_Assays'] = f'other_data__variant_data__{json_key}'


# Old column name -> Report_in_ClinVar/Report_in_LOVD field name. Used both
# to build the most-recent-wins Subquery mapping below (for index()'s
# flattened table view) and directly by variant_reports() (for the detail
# page's per-report list, unflattened - see views.py).
CLINVAR_REPORT_FIELDS = {
    'Allele_Origin_ClinVar': 'Allele_Origin',
    'Clinical_Significance_ClinVar': 'Clinical_Significance',
    'DateSignificanceLastEvaluated_ClinVar': 'DateSignificanceLastEvaluated',
    'Date_Last_Updated_ClinVar': 'Date_Last_Updated',
    'Description_ClinVar': 'Description',
    'Method_ClinVar': 'Method',
    'Review_Status_ClinVar': 'Review_Status',
    'SCV_ClinVar': 'SCV',
    'Submitter_ClinVar': 'Submitter',
    'Summary_Evidence_ClinVar': 'Summary_Evidence',
}

LOVD_REPORT_FIELDS = {
    'Classification_LOVD': 'Classification',
    'Created_date_LOVD': 'Created_date',
    'DBID_LOVD': 'DBID',
    'Edited_date_LOVD': 'Edited_date',
    'Functional_analysis_result_LOVD': 'Functional_analysis_result',
    'Functional_analysis_technique_LOVD': 'Functional_analysis_technique',
    'Genetic_origin_LOVD': 'Genetic_origin',
    'Individuals_LOVD': 'Individuals',
    'Remarks_LOVD': 'Remarks',
    'Submitters_LOVD': 'Submitters',
    'Variant_frequency_LOVD': 'Variant_frequency',
}

# Variant_haplotype_LOVD lives on the LOVD anchor table (Variant_in_LOVD),
# not per-report, but the frontend's reportBinding.cols still expects it on
# each LOVD report item.
LOVD_ANCHOR_REPORT_FIELDS = {'Variant_haplotype_LOVD': 'Variant_haplotype'}


def _build_map():
    m = {}

    base_fields = [
        'Gene_Symbol', 'Reference_Sequence', 'HGVS_cDNA', 'BIC_Nomenclature', 'HGVS_Protein',
        'Protein_Change', 'CA_ID', 'Synonyms',
    ]
    for f in base_fields:
        m[f] = f

    # ClinVar/LOVD (Report_in_ClinVar/Report_in_LOVD, one-to-many, most-recently-updated wins)
    for old, field in CLINVAR_REPORT_FIELDS.items():
        m[old] = _clinvar_subquery(field)
    for old, field in LOVD_REPORT_FIELDS.items():
        m[old] = _lovd_subquery(field)

    m.update(_RENAMES)
    m.update(_build_gnomad_columns())
    return m


COLUMN_MAP = _build_map()

# Confirmed genuinely dropped from the finalized schema, and columns whose
# mapping is still an open question (aggregate/computed fields with no
# single-field home) - both rejected explicitly rather than silently
# returning null, so a stale frontend request is loud, not quietly wrong.
DROPPED_COLUMNS = {
    'Ethnicity_BIC', 'Germline_or_Somatic_BIC', 'Literature_citation_BIC', 'Mutation_type_BIC',
    'Number_of_family_member_carrying_mutation_BIC', 'Patient_nationality_BIC',
    'URL_ENIGMA', 'RNA_LOVD', 'HGVS_RNA', 'Mupit_Structure', 'Max_Allele_Frequency',
}
PENDING_COLUMNS = {
    'Pathogenicity_all', 'Allele_Frequency', 'Allele_Frequency_Charts_Exome_GnomAD',
    'Allele_Frequency_Charts_Genome_GnomADv3', 'Source', 'Source_URL',
}

# Every resolvable column except ClinVar/LOVD, which the detail page renders
# from variant_reports() instead (see views.py) - those are genuinely
# one-to-many and the flat variant() response only has room for one
# (most-recently-updated) value per column, which isn't the right answer here.
DETAIL_COLUMNS = sorted(
    set(COLUMN_MAP) - set(CLINVAR_REPORT_FIELDS) - set(LOVD_REPORT_FIELDS)
)


def resolve_columns(queryset, column_names):
    """Returns (annotated_queryset, values_field_names) for the given old
    column names, or raises KeyError with the unresolvable names."""
    unresolvable = [c for c in column_names if c not in COLUMN_MAP]
    if unresolvable:
        raise KeyError(unresolvable)

    values_fields = []
    for name in column_names:
        spec = COLUMN_MAP[name]
        if isinstance(spec, str) and spec == name:
            # old column name already matches the real field name (base
            # fields) - no aliasing needed.
            values_fields.append(name)
            continue

        # apply_filters and index()'s column= selection can both resolve the
        # same column - annotating it twice would raise, so skip if it's
        # already on the queryset. Always annotate under the *old* column
        # name so JSON output keys match what the frontend expects, even for
        # a plain renamed/related path.
        if name not in queryset.query.annotations:
            if isinstance(spec, str):
                queryset = queryset.annotate(**{name: F(spec)})
            else:
                _, build = spec
                queryset = queryset.annotate(**{name: build('VRS_Digest')})
        values_fields.append(name)
    return queryset, values_fields
