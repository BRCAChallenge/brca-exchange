'use strict';

var React = require('react');
var _ = require('underscore');

function renderCell(val) {
    return <span>{val}</span>;
}

const filterColumns = [
    {name: 'Gene', prop: 'Gene_Symbol', values: ['BRCA1', 'BRCA2']},
    {name: 'Pathogenicity', prop: 'Pathogenicity_expert', values: ['Pathogenic', 'Benign / Little Clinical Significance', 'Not Yet Reviewed']}
];

const expertModeGroups = [
    {groupTitle: 'Variant Nomenclature', internalGroupName: 'Variant Nomenclature', innerCols: [
        {title: 'Gene', prop: 'Gene_Symbol', render: gene => <i>{gene}</i>},
        {title: 'HGVS Nucleotide', prop: 'HGVS_cDNA', render: nucleotide => nucleotide.split(':')[1]},
        {title: 'Transcript Identifier', prop: 'Reference_Sequence'},
        {title: 'HGVS RNA', prop: 'HGVS_RNA'},
        {title: 'HGVS Protein', prop: 'HGVS_Protein', render: protein => protein.split(':')[1]},
        // Protein Identfifier is pulled from HGVS_Protein, this is handled in VariantDetail (index.js)
        {title: 'Protein Identifier', prop: 'HGVS_Protein_ID'},
        {title: 'Protein Abbrev', prop: 'Protein_Change'},
        {title: 'BIC Designation', prop: 'BIC_Nomenclature'},
        {title: 'Genomic Nomenclature (GRCh38)', prop: 'Genomic_Coordinate_hg38'},
        {title: 'Genomic Nomenclature (GRCh37)', prop: 'Genomic_Coordinate_hg37'}
    ]},

    {groupTitle: 'Clinical Significance (ENIGMA)', internalGroupName: 'Significance (ENIGMA)', innerCols: [
        {title: 'Clinical Significance', prop: 'Pathogenicity_expert'},
        {title: 'IARC Class', prop: 'Clinical_significance_ENIGMA'},
        {title: 'Comment on Clinical Significance', prop: 'Comment_on_clinical_significance_ENIGMA'},
        {title: 'Clinical Significance Citations', prop: 'Clinical_significance_citations_ENIGMA'},
        {title: 'Supporting Evidence URL(s)', prop: 'URL_ENIGMA'},
        {title: 'Date Last Evaluated', prop: 'Date_last_evaluated_ENIGMA'},
        {title: 'Assertion Method', prop: 'Assertion_method_ENIGMA'},
        {title: 'Assertion Method Citation', prop: 'Assertion_method_citation_ENIGMA'},
        {title: 'Allele Origin', prop: 'Allele_origin_ENIGMA'},
        {title: 'ClinVar Accession', prop: 'ClinVarAccession_ENIGMA'}
    ]},
];

const researchModeGroups = [
    {groupTitle: 'Variant Nomenclature', internalGroupName: 'Variant Nomenclature', innerCols: [
        {title: 'Gene Symbol', prop: 'Gene_Symbol', render: gene => <i>{gene}</i>, core: true},
        {title: 'Reference cDNA Sequence', prop: 'Reference_Sequence', core: true},
        {title: 'HGVS Nucleotide', prop: 'HGVS_cDNA', render: nucleotide => nucleotide.split(':')[1], core: true},
        {title: 'HGVS Protein', prop: 'HGVS_Protein', render: protein => protein.split(':')[1], core: true},
        {title: 'Protein Amino Acid Change', prop: 'Protein_Change', core: true},
        {title: 'BIC Designation', prop: 'BIC_Nomenclature', core: true},
        {title: 'Genome (GRCh38)', prop: 'Genomic_Coordinate_hg38', core: true},
        {title: 'Genome (GRCh37)', prop: 'Genomic_Coordinate_hg37'},
        {title: 'Genome (GRCh36)', prop: 'Genomic_Coordinate_hg36'},
        {title: 'RNA (LOVD)', prop: 'RNA_LOVD'}
    ]},

    {groupTitle: 'Clinical Significance (ENIGMA)', internalGroupName: 'Significance (ENIGMA)', innerCols: [
        {title: 'Clinical Significance', prop: 'Clinical_significance_ENIGMA', core: true},
        {title: 'Comment on Clinical Significance', prop: 'Comment_on_clinical_significance_ENIGMA', core: true},
        {title: 'Assertion Method', prop: 'Assertion_method_ENIGMA', core: true},
        {title: 'Date last evaluated', prop: 'Date_last_evaluated_ENIGMA', core: true},
        {title: 'Collection Method', prop: 'Collection_method_ENIGMA', core: true},
        {title: 'Clinical Significance Citation', prop: 'Clinical_significance_citations_ENIGMA', core: true},
        {title: 'Allele Origin', prop: 'Allele_origin_ENIGMA', core: true},
    ]},

    {groupTitle: 'Clinical Significance (ClinVar)', internalGroupName: 'Significance (ClinVar)', innerCols: [
        {title: 'Clinical Significance', prop: 'Clinical_Significance_ClinVar', core: true},
        {title: 'Submitter', prop: 'Submitter_ClinVar', core: true},
        {title: 'Analysis Method', prop: 'Method_ClinVar', core: true},
        {title: 'Date last updated', prop: 'Date_Last_Updated_ClinVar', core: true},
        {title: 'SCV Accession', prop: 'SCV_ClinVar', core: true},
        {title: 'Allele Origin', prop: 'Allele_Origin_ClinVar', core: true},
    ]},

    {groupTitle: 'Clinical Significance (LOVD)', internalGroupName: 'Significance (LOVD)', innerCols: [
        {title: 'Variant Frequency', prop: 'Variant_frequency_LOVD'},
        {title: 'Variant Haplotype', prop: 'Variant_haplotype_LOVD'},
        {title: 'Submitters', prop: 'Submitters_LOVD'},
        {title: 'Genetic Origin', prop: 'Genetic_origin_LOVD'},
        {title: 'Individuals', prop: 'Individuals_LOVD'},
        {title: 'Variant Effect', prop: 'Variant_effect_LOVD'},
    ]},

    {groupTitle: 'Clinical Significance (BIC)', internalGroupName: 'Significance (BIC)', innerCols: [
        {title: 'Clinical Significance', prop: 'Clinical_classification_BIC', core: true},
        {title: 'Clinical Importance', prop: 'Clinical_importance_BIC', core: true},
        {title: 'Patient Nationality', prop: 'Patient_nationality_BIC'},
        {title: 'Ethnicity', prop: 'Ethnicity_BIC'},
        {title: 'Family members carrying this variant', prop: 'Number_of_family_member_carrying_mutation_BIC'},
        {title: 'Literature Reference', prop: 'Literature_citation_BIC', core: true},
        {title: 'Allele Origin', prop: 'Germline_or_Somatic_BIC'},
    ]},

    {groupTitle: 'Allele Frequency Reference Sets', internalGroupName: 'Allele Frequency Reference Sets', innerCols: [
        {title: 'Maximum Allele Frequency (1000 Genomes and ESP)', prop: 'Max_Allele_Frequency', core: true},
        {title: 'Allele Frequency (1000 Genomes)', prop: 'Allele_frequency_1000_Genomes', core: true},
        {title: 'AFR Allele Frequency (1000 Genomes)', prop: 'AFR_Allele_frequency_1000_Genomes', core: true},
        {title: 'AMR Allele Frequency (1000 Genomes)', prop: 'AMR_Allele_frequency_1000_Genomes', core: true},
        {title: 'EAS Allele Frequency (1000 Genomes)', prop: 'EAS_Allele_frequency_1000_Genomes', core: true},
        {title: 'EUR Allele Frequency (1000 Genomes)', prop: 'EUR_Allele_frequency_1000_Genomes', core: true},
        {title: 'SAS Allele Frequency (1000 Genomes)', prop: 'SAS_Allele_frequency_1000_Genomes', core: true},
        {title: 'EA Allele Frequency (ESP)', prop: 'EA_Allele_Frequency_ESP', core: true},
        {title: 'AA Allele Frequency (ESP)', prop: 'AA_Allele_Frequency_ESP', core: true},
        {title: 'Allele Frequency (ESP)', prop: 'Allele_Frequency_ESP', core: true},
        {title: 'Allele frequency (ExAC minus TCGA)', prop: 'Allele_frequency_ExAC', core: true},
        {title: 'AFR Allele frequency (ExAC minus TCGA)', prop: 'Allele_frequency_AFR_ExAC', core: true},
        {title: 'AMR Allele frequency (ExAC minus TCGA)', prop: 'Allele_frequency_AMR_ExAC', core: true},
        {title: 'EAS Allele frequency (ExAC minus TCGA)', prop: 'Allele_frequency_EAS_ExAC', core: true},
        {title: 'FIN Allele frequency (ExAC minus TCGA)', prop: 'Allele_frequency_FIN_ExAC', core: true},
        {title: 'NFE Allele frequency (ExAC minus TCGA)', prop: 'Allele_frequency_NFE_ExAC', core: true},
        {title: 'OTH Allele frequency (ExAC minus TCGA)', prop: 'Allele_frequency_OTH_ExAC', core: true},
        {title: 'SAS Allele frequency (ExAC minus TCGA)', prop: 'Allele_frequency_SAS_ExAC', core: true},
    ]},

    {groupTitle: 'Allele Counts (ExAC minus TCGA)', internalGroupName: 'Allele Counts (ExAC minus TCGA)', innerCols: [
        {title: 'AFR Allele count', prop: 'Allele_count_AFR_ExAC', core: true},
        {title: 'AMR Allele count', prop: 'Allele_count_AMR_ExAC', core: true},
        {title: 'EAS Allele count', prop: 'Allele_count_EAS_ExAC', core: true},
        {title: 'FIN Allele count', prop: 'Allele_count_FIN_ExAC', core: true},
        {title: 'NFE Allele count', prop: 'Allele_count_NFE_ExAC', core: true},
        {title: 'OTH Allele count', prop: 'Allele_count_OTH_ExAC', core: true},
        {title: 'SAS Allele count', prop: 'Allele_count_SAS_ExAC', core: true},
    ]},

    {groupTitle: 'Multifactorial Likelihood Analysis', internalGroupName: 'Multifactorial Likelihood Analysis', innerCols: [
        {title: 'Posterior probability of pathogenicity', prop: 'Posterior_probability_exLOVD', core: true},
        {title: 'Prior probability of pathogenicity', prop: 'Combined_prior_probablility_exLOVD', core: true},
        {title: 'Missense analysis pathogenicity prior', prop: 'Missense_analysis_prior_probability_exLOVD', core: true},
        {title: 'Co-occurrence likelihood', prop: 'Co_occurrence_LR_exLOVD', core: true},
        {title: 'Segregation Likelihood Ratio', prop: 'Segregation_LR_exLOVD', core: true},
        {title: 'Summary Family History Likelihood Ratio', prop: 'Sum_family_LR_exLOVD', core: true},
        {title: 'Literature Reference', prop: 'Literature_source_exLOVD', core: true}
    ]},
];

// subColumns populate the column selection checkboxes.
// They should match the variant detail groupings.
const subColumns = _.map(researchModeGroups, function (group) {
    return {
        subColTitle: group.groupTitle,
        subColList: _.map(group.innerCols, function (innerCol) {
            return {
                title: innerCol.title,
                prop: innerCol.prop,
                render: renderCell
            };
        })
    };
});

const columns = [
    {title: 'Gene', prop: 'Gene_Symbol', render: gene => <i>{gene}</i>},
    {title: 'HGVS Nucleotide', prop: 'HGVS_cDNA', render: nucleotide => nucleotide.split(':')[1]},
    {title: 'Transcript Identifier', prop: 'Reference_Sequence'},
    {title: 'HGVS RNA', prop: 'HGVS_RNA'},
    {title: 'HGVS Protein', prop: 'HGVS_Protein', render: protein => protein.split(':')[1]},
    // Protein Identfifier is pulled from HGVS_Protein, this is handled in VariantDetail (index.js)
    {title: 'Protein Identifier', prop: 'HGVS_Protein_ID'},
    {title: 'Protein Abbrev', prop: 'Protein_Change'},
    {title: 'BIC Designation', prop: 'BIC_Nomenclature'},
    {title: 'Genomic Nomenclature (GRCh38)', prop: 'Genomic_Coordinate_hg38'},
    {title: 'Genomic Nomenclature (GRCh37)', prop: 'Genomic_Coordinate_hg37'},
    {title: 'Clinical Significance', prop: 'Pathogenicity_expert'},
    {title: 'IARC Class', prop: 'Clinical_significance_ENIGMA'},
    {title: 'Comment on Clinical Significance', prop: 'Comment_on_clinical_significance_ENIGMA'},
    {title: 'Clinical Significance Citations', prop: 'Clinical_significance_citations_ENIGMA'},
    {title: 'Supporting Evidence URL(s)', prop: 'URL_ENIGMA'},
    {title: 'Date Last Evaluated', prop: 'Date_last_evaluated_ENIGMA'},
    {title: 'Assertion Method', prop: 'Assertion_method_ENIGMA'},
    {title: 'Assertion Method Citation', prop: 'Assertion_method_citation_ENIGMA'},
    {title: 'Allele Origin', prop: 'Allele_origin_ENIGMA'},
    {title: 'ClinVar Accession', prop: 'ClinVarAccession_ENIGMA'}
];

const researchModeColumns = [
    {title: 'Gene Symbol', prop: 'Gene_Symbol', render: gene => <i>{gene}</i>},
    {title: 'Genome (GRCh36)', prop: 'Genomic_Coordinate_hg36'},
    {title: 'Genome (GRCh37)', prop: 'Genomic_Coordinate_hg37'},
    {title: 'Genome (GRCh38)', prop: 'Genomic_Coordinate_hg38'},
    {title: 'Mutation category (BIC)', prop: 'Mutation_type_BIC'},
    {title: 'SIFT score', prop: 'Sift_Score'},
    {title: 'BIC Variant Identifier', prop: 'BIC_Nomenclature'},
    {title: 'Nucleotide', prop: 'HGVS_cDNA'},
    {title: 'Protein', prop: 'HGVS_Protein'},
    {title: 'SCV Accession (ClinVar)', prop: 'SCV_ClinVar'},
    {title: 'Source(s)', prop: 'Source'},
    {title: 'Source URL(s)', prop: 'Source_URL'},
    {title: 'Synonyms', prop: 'Synonyms'},
    {title: 'Protein Amino Acid Change', prop: 'Protein_Change'},
    {title: 'Reference cDNA Sequence', prop: 'Reference_Sequence'},
    {title: 'RNA (LOVD)', prop: 'RNA_LOVD'},
    {title: 'Submitters (LOVD)', prop: 'Submitters_LOVD'},
    {title: 'Genetic Origin (LOVD)', prop: 'Genetic_origin_LOVD'},
    {title: 'Individuals (LOVD)', prop: 'Individuals_LOVD'},
    {title: 'Variant Effect (LOVD)', prop: 'Variant_effect_LOVD'},
    {title: 'Allele Origin (ClinVar)', prop: 'Allele_Origin_ClinVar'},
    {title: 'Allele Origin (ENIGMA)', prop: 'Allele_origin_ENIGMA'},
    {title: 'Ethnicity (BIC)', prop: 'Ethnicity_BIC'},
    {title: 'Allele Origin (BIC)', prop: 'Germline_or_Somatic_BIC'},
    {title: 'Patient Nationality (BIC)', prop: 'Patient_nationality_BIC'},
    {title: 'Variant Haplotype (LOVD)', prop: 'Variant_haplotype_LOVD'},
    {title: 'Family members carrying this variant (BIC)', prop: 'Number_of_family_member_carrying_mutation_BIC'},
    {title: 'Co-occurrence likelihood (exLOVD)', prop: 'Co_occurrence_LR_exLOVD'},
    {title: 'Prior probability of pathogenicity (exLOVD)', prop: 'Combined_prior_probablility_exLOVD'},
    {
        title: 'Missense analysis probability of pathogenicity (exLOVD)',
        prop: 'Missense_analysis_prior_probability_exLOVD'
    },
    {title: 'Probability of pathogenicity (exLOVD)', prop: 'Posterior_probability_exLOVD'},
    {title: 'Segregation Likelihood Ratio (exLOVD)', prop: 'Segregation_LR_exLOVD'},
    {title: 'Summary Family History Likelihood Ratio (exLOVD)', prop: 'Sum_family_LR_exLOVD'},
    {title: 'Assertion Method (ENIGMA)', prop: 'Assertion_method_citation_ENIGMA'},
    {title: 'Clinical Significance Citation (ENIGMA)', prop: 'Clinical_significance_citations_ENIGMA'},
    {title: 'Literature Reference (BIC)', prop: 'Literature_citation_BIC'},
    {title: 'Literature Reference (exLOVD)', prop: 'Literature_source_exLOVD'},
    {title: 'Pathogenicity', prop: 'Pathogenicity_all'},
    {title: 'Assertion Method (ENIGMA)', prop: 'Assertion_method_ENIGMA'},
    {title: 'Clinical Significance (BIC)', prop: 'Clinical_classification_BIC'},
    {title: 'Clinical Importance (BIC)', prop: 'Clinical_importance_BIC'},
    {title: 'Clinical Significance (ClinVar)', prop: 'Clinical_Significance_ClinVar'},
    {title: 'Clinical Significance (ENIGMA)', prop: 'Clinical_significance_ENIGMA'},
    {title: 'Collection Method (ENIGMA)', prop: 'Collection_method_ENIGMA'},
    {title: 'Comment on Clinical Significance (ENIGMA)', prop: 'Comment_on_clinical_significance_ENIGMA'},
    {title: 'Date last evaluated (ENIGMA)', prop: 'Date_last_evaluated_ENIGMA'},
    {title: 'Date last updated (ClinVar)', prop: 'Date_Last_Updated_ClinVar'},
    {title: 'Has Discordant Evidence', prop: 'Discordant'},
    {title: 'Functional Analysis Result (LOVD)', prop: 'Functional_analysis_result_LOVD'},
    {title: 'Functional Analysis Method (LOVD)', prop: 'Functional_analysis_technique_LOVD'},
    {title: 'Analysis Method (ClinVar)', prop: 'Method_ClinVar'},
    {title: 'ClinVar Accession', prop: 'ClinVarAccession_ENIGMA'},
    {title: 'Condition Category (ENIGMA)', prop: 'Condition_category_ENIGMA'},
    {title: 'Condition ID Type (ENIGMA)', prop: 'Condition_ID_type_ENIGMA'},
    {title: 'Condition ID Value (ENIGMA)', prop: 'Condition_ID_value_ENIGMA'},
    {title: 'Submitter (ClinVar)', prop: 'Submitter_ClinVar'},
    {title: 'URL (ENIGMA)', prop: 'URL_ENIGMA'},
    {title: 'AFR Allele Frequency (1000 Genomes)', prop: 'AFR_Allele_frequency_1000_Genomes'},
    {title: 'Allele Frequency', prop: 'Allele_Frequency'},
    {title: 'Allele Frequency (1000 Genomes)', prop: 'Allele_frequency_1000_Genomes'},
    {title: 'AMR Allele Frequency (1000 Genomes)', prop: 'AMR_Allele_frequency_1000_Genomes'},
    {title: 'EAS Allele Frequency (1000 Genomes)', prop: 'EAS_Allele_frequency_1000_Genomes'},
    {title: 'EUR Allele Frequency (1000 Genomes)', prop: 'EUR_Allele_frequency_1000_Genomes'},
    {title: 'Maximum Allele Frequency', prop: 'Max_Allele_Frequency'},
    {title: 'EA Allele Frequency (ESP)', prop: 'EA_Allele_Frequency_ESP'},
    {title: 'AA Allele Frequency (ESP)', prop: 'AA_Allele_Frequency_ESP'},
    {title: 'Allele Frequency (ESP)', prop: 'Allele_Frequency_ESP'},
    {title: 'Allele frequency (ExAC minus TCGA)', prop: 'Allele_frequency_ExAC'},
    {title: 'AFR Allele count (ExAC minus TCGA)', prop: 'Allele_count_AFR_ExAC'},
    {title: 'AFR Allele number (ExAC minus TCGA)', prop: 'Allele_number_AFR_ExAC'},
    {title: 'AFR Allele frequency (ExAC minus TCGA)', prop: 'Allele_frequency_AFR_ExAC'},
    {title: 'AFR Homozygous count (ExAC minus TCGA)', prop: 'Homozygous_count_AFR_ExAC'},
    {title: 'AMR Allele count (ExAC minus TCGA)', prop: 'Allele_count_AMR_ExAC'},
    {title: 'AMR Allele number (ExAC minus TCGA)', prop: 'Allele_number_AMR_ExAC'},
    {title: 'AMR Allele frequency (ExAC minus TCGA)', prop: 'Allele_frequency_AMR_ExAC'},
    {title: 'AMR Homozygous count (ExAC minus TCGA)', prop: 'Homozygous_count_AMR_ExAC'},
    {title: 'EAS Allele count (ExAC minus TCGA)', prop: 'Allele_count_EAS_ExAC'},
    {title: 'EAS Allele number (ExAC minus TCGA)', prop: 'Allele_number_EAS_ExAC'},
    {title: 'EAS Allele frequency (ExAC minus TCGA)', prop: 'Allele_frequency_EAS_ExAC'},
    {title: 'EAS Homozygous count (ExAC minus TCGA)', prop: 'Homozygous_count_EAS_ExAC'},
    {title: 'FIN Allele count (ExAC minus TCGA)', prop: 'Allele_count_FIN_ExAC'},
    {title: 'FIN Allele number (ExAC minus TCGA)', prop: 'Allele_number_FIN_ExAC'},
    {title: 'FIN Allele frequency (ExAC minus TCGA)', prop: 'Allele_frequency_FIN_ExAC'},
    {title: 'FIN Homozygous count (ExAC minus TCGA)', prop: 'Homozygous_count_FIN_ExAC'},
    {title: 'NFE Allele count (ExAC minus TCGA)', prop: 'Allele_count_NFE_ExAC'},
    {title: 'NFE Allele number (ExAC minus TCGA)', prop: 'Allele_number_NFE_ExAC'},
    {title: 'NFE Allele frequency (ExAC minus TCGA)', prop: 'Allele_frequency_NFE_ExAC'},
    {title: 'NFE Homozygous count (ExAC minus TCGA)', prop: 'Homozygous_count_NFE_ExAC'},
    {title: 'OTH Allele count (ExAC minus TCGA)', prop: 'Allele_count_OTH_ExAC'},
    {title: 'OTH Allele number (ExAC minus TCGA)', prop: 'Allele_number_OTH_ExAC'},
    {title: 'OTH Allele frequency (ExAC minus TCGA)', prop: 'Allele_frequency_OTH_ExAC'},
    {title: 'OTH Homozygous count (ExAC minus TCGA)', prop: 'Homozygous_count_OTH_ExAC'},
    {title: 'SAS Allele count (ExAC minus TCGA)', prop: 'Allele_count_SAS_ExAC'},
    {title: 'SAS Allele number (ExAC minus TCGA)', prop: 'Allele_number_SAS_ExAC'},
    {title: 'SAS Allele frequency (ExAC minus TCGA)', prop: 'Allele_frequency_SAS_ExAC'},
    {title: 'SAS Homozygous count (ExAC minus TCGA)', prop: 'Homozygous_count_SAS_ExAC'},
    {title: 'SAS Allele Frequency (1000 Genomes)', prop: 'SAS_Allele_frequency_1000_Genomes'},
    {title: 'Variant Frequency (LOVD)', prop: 'Variant_frequency_LOVD'}
];

module.exports = ({
    filterColumns,
    expertModeGroups,
    researchModeGroups,
    subColumns,
    columns,
    researchModeColumns
});
