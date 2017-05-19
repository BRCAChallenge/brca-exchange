// A table for variants.
//
// The intent here was to split the generic table code
// in DataTable from the variant domain knowledge, which
// would be here. That division has broken down due to
// peculiarities of react-data-components DataMixin, and
// time pressure. Knowledge about variants is in both files.
// This needs to be revisited.

/*global module: false, require: false, window: false */
'use strict';

var React = require('react');
var PureRenderMixin = require('./PureRenderMixin');
var DataTable = require('./DataTable');
var _ = require('underscore');
var {Col, Panel, Button, Input} = require('react-bootstrap');
var ColumnCheckbox = require('./ColumnCheckbox');
var {getDefaultExpertColumns, getDefaultResearchColumns, getAllSources} = require('./VariantTableDefaults');
var {State} = require('react-router');

require('react-data-components-bd2k/css/table-twbs.css');

function buildHeader(onClick, title) {
    return (
        <span>
            {title}
            <span onClick={ev => {ev.stopPropagation(); onClick(title); }}
                  className='help glyphicon glyphicon-question-sign superscript'/>
        </span>
    );
}

function renderCell(val) {
    return <span>{val}</span>;
}

var filterColumns = [
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
        {title: 'African Allele Frequency (1000 Genomes)', prop: 'AFR_Allele_frequency_1000_Genomes', core: true},
        {title: 'AMR Allele Frequency (1000 Genomes)', prop: 'AMR_Allele_frequency_1000_Genomes', core: true},
        {title: 'EAS Allele Frequency (1000 Genomes)', prop: 'EAS_Allele_frequency_1000_Genomes', core: true},
        {title: 'EUR Allele Frequency (1000 Genomes)', prop: 'EUR_Allele_frequency_1000_Genomes', core: true},
        {title: 'South Asian Allele Frequency (1000 Genomes)', prop: 'SAS_Allele_frequency_1000_Genomes', core: true},
        {title: 'Allele Frequency (ExAC minus TCGA)', prop: 'Allele_frequency_ExAC', core: true},
        {title: '% Allele Frequencies: EA|AA|All (ESP)', prop: 'Minor_allele_frequency_ESP', core: true},
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

var columns = [
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

var researchModeColumns = [

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

    {title: 'African Allele Frequency (1000 Genomes)', prop: 'AFR_Allele_frequency_1000_Genomes'},
    {title: 'Allele Frequency', prop: 'Allele_Frequency'},
    {title: 'Allele Frequency (1000 Genomes)', prop: 'Allele_frequency_1000_Genomes'},
    {title: 'Allele Frequency (ExAC)', prop: 'Allele_frequency_ExAC'},
    {title: 'AMR Allele Frequency (1000 Genomes)', prop: 'AMR_Allele_frequency_1000_Genomes'},
    {title: 'EAS Allele Frequency (1000 Genomes)', prop: 'EAS_Allele_frequency_1000_Genomes'},
    {title: 'EUR Allele Frequency (1000 Genomes)', prop: 'EUR_Allele_frequency_1000_Genomes'},
    {title: 'Maximum Allele Frequency', prop: 'Max_Allele_Frequency'},
    {title: 'Allele Frequencies: EA|AA|All (ESP)', prop: 'Minor_allele_frequency_ESP'},
    {title: 'South Asian Allele Frequency (1000 Genomes)', prop: 'SAS_Allele_frequency_1000_Genomes'},
    {title: 'Variant Frequency (LOVD)', prop: 'Variant_frequency_LOVD'}

];

var subColumns = [
    {
        subColTitle: "Variant Nomenclature",
        subColList: [
            {title: 'BIC Variant Identifier', prop: 'BIC_Nomenclature', render: renderCell},
            {title: 'Protein', prop: 'HGVS_Protein', render: renderCell},
            {title: 'SCV Accession (ClinVar)', prop: 'SCV_ClinVar', render: renderCell},
            {title: 'HGVS Nucleotide', prop: 'HGVS_cDNA', render: renderCell},
            {title: 'Protein Amino Acid Change', prop: 'Protein_Change', render: renderCell},
            {title: 'Reference cDNA Sequence', prop: 'Reference_Sequence', render: renderCell}
        ]
    },
    {
        subColTitle: "Origin",
        subColList: [
            {title: 'Allele Origin (ClinVar)', prop: 'Allele_Origin_ClinVar', render: renderCell},
            {title: 'Allele Origin (ENIGMA)', prop: 'Allele_origin_ENIGMA', render: renderCell},
            {title: 'Ethnicity (BIC)', prop: 'Ethnicity_BIC', render: renderCell},
            {title: 'Allele Origin (BIC)', prop: 'Germline_or_Somatic_BIC', render: renderCell},
            {title: 'Patient Nationality (BIC)', prop: 'Patient_nationality_BIC', render: renderCell},
            {title: 'Variant Haplotype (LOVD)', prop: 'Variant_haplotype_LOVD', render: renderCell}
        ]
    },

    {
        subColTitle: "Frequency",
        subColList: [
            {
                title: 'African Allele Frequency (1000 Genomes)',
                prop: 'AFR_Allele_frequency_1000_Genomes',
                render: renderCell
            },
            {title: 'Allele Frequency', prop: 'Allele_Frequency', render: renderCell},
            {title: 'Allele Frequency (1000 Genomes)', prop: 'Allele_frequency_1000_Genomes', render: renderCell},
            {title: 'Allele Frequency (ExAC)', prop: 'Allele_frequency_ExAC', render: renderCell},
            {
                title: 'AMR Allele Frequency (1000 Genomes)',
                prop: 'AMR_Allele_frequency_1000_Genomes',
                render: renderCell
            },
            {
                title: 'EAS Allele Frequency (1000 Genomes)',
                prop: 'EAS_Allele_frequency_1000_Genomes',
                render: renderCell
            },
            {
                title: 'EUR Allele Frequency (1000 Genomes)',
                prop: 'EUR_Allele_frequency_1000_Genomes',
                render: renderCell
            },
            {title: 'Maximum Allele Frequency', prop: 'Max_Allele_Frequency', render: renderCell},
            {title: 'Allele Frequencies: EA|AA|All (ESP)', prop: 'Minor_allele_frequency_ESP', render: renderCell},
            {
                title: 'South Asian Allele Frequency (1000 Genomes)',
                prop: 'SAS_Allele_frequency_1000_Genomes',
                render: renderCell
            },
            {title: 'Variant Frequency (LOVD)', prop: 'Variant_frequency_LOVD', render: renderCell}
        ]
    },

    {
        subColTitle: "Genomic",
        subColList: [
            {title: 'Gene Symbol', prop: 'Gene_Symbol', render: renderCell},
            {title: 'Genome (GRCh38)', prop: 'Genomic_Coordinate_hg38', render: renderCell},
            {title: 'Genome (GRCh36)', prop: 'Genomic_Coordinate_hg36', render: renderCell},
            {title: 'Genome (GRCh37)', prop: 'Genomic_Coordinate_hg37', render: renderCell}
        ]
    },
    {
        subColTitle: "Bioinformatic Annotation",
        subColList: [
            {title: 'Mutation category (BIC)', prop: 'Mutation_type_BIC', render: renderCell},
            {title: 'SIFT score', prop: 'Sift_Score', render: renderCell}
        ]
    },
    {
        subColTitle: "Probability",
        subColList: [
            {title: 'Co-occurrence likelihood (exLOVD)', prop: 'Co_occurrence_LR_exLOVD', render: renderCell},
            {
                title: 'Prior probability of pathogenicity (exLOVD)',
                prop: 'Combined_prior_probablility_exLOVD',
                render: renderCell
            },
            {
                title: 'Missense analysis probability of pathogenicity (exLOVD)',
                prop: 'Missense_analysis_prior_probability_exLOVD',
                render: renderCell
            },
            {title: 'Probability of pathogenicity (exLOVD)', prop: 'Posterior_probability_exLOVD', render: renderCell},
            {title: 'Segregation Likelihood Ratio (exLOVD)', prop: 'Segregation_LR_exLOVD', render: renderCell},
            {
                title: 'Summary Family History Likelihood Ratio (exLOVD)',
                prop: 'Sum_family_LR_exLOVD',
                render: renderCell
            }
        ]
    },
    {
        subColTitle: "Significance",
        subColList: [
            {title: 'Pathogenicity', prop: 'Pathogenicity_all', render: renderCell},
            {title: 'Assertion Method (ENIGMA)', prop: 'Assertion_method_ENIGMA', render: renderCell},
            {title: 'Clinical Significance (BIC)', prop: 'Clinical_classification_BIC', render: renderCell},
            {title: 'Clinical Importance (BIC)', prop: 'Clinical_importance_BIC', render: renderCell},
            {title: 'Clinical Significance (ClinVar)', prop: 'Clinical_Significance_ClinVar', render: renderCell},
            {title: 'Clinical Significance (ENIGMA)', prop: 'Clinical_significance_ENIGMA', render: renderCell},
            {title: 'Collection Method (ENIGMA)', prop: 'Collection_method_ENIGMA', render: renderCell},
            {
                title: 'Comment on Clinical Significance (ENIGMA)',
                prop: 'Comment_on_clinical_significance_ENIGMA',
                render: renderCell
            },
            {title: 'Date last evaluated (ENIGMA)', prop: 'Date_last_evaluated_ENIGMA', render: renderCell},
            {title: 'Date last updated (ClinVar)', prop: 'Date_Last_Updated_ClinVar', render: renderCell},
            {title: 'Has Discordant Evidence', prop: 'Discordant', render: renderCell},
            {title: 'Functional Analysis Result (LOVD)', prop: 'Functional_analysis_result_LOVD', render: renderCell},
            {
                title: 'Functional Analysis Method (LOVD)',
                prop: 'Functional_analysis_technique_LOVD',
                render: renderCell
            },
            {title: 'Analysis Method (ClinVar)', prop: 'Method_ClinVar', render: renderCell}
        ]
    },
    {
        subColTitle: "Pedigree",
        subColList: [
            {
                title: 'Family members carrying this variant (BIC)',
                prop: 'Number_of_family_member_carrying_mutation_BIC',
                render: renderCell
            }
        ]
    },
    {
        subColTitle: "Publications",
        subColList: [
            {title: 'Assertion Method (ENIGMA)', prop: 'Assertion_method_citation_ENIGMA', render: renderCell},
            {
                title: 'Clinical Significance Citation (ENIGMA)',
                prop: 'Clinical_significance_citations_ENIGMA',
                render: renderCell
            },
            {title: 'Literature Reference (BIC)', prop: 'Literature_citation_BIC', render: renderCell},
            {title: 'Literature Reference (exLOVD)', prop: 'Literature_source_exLOVD', render: renderCell}
        ]
    },
    {
        subColTitle: "Source",
        subColList: [
            {title: 'ClinVar Accession', prop: 'ClinVarAccession_ENIGMA', render: renderCell},
            {title: 'Condition Category (ENIGMA)', prop: 'Condition_category_ENIGMA', render: renderCell},
            {title: 'Condition ID Type (ENIGMA)', prop: 'Condition_ID_type_ENIGMA', render: renderCell},
            {title: 'Condition ID Value (ENIGMA)', prop: 'Condition_ID_value_ENIGMA', render: renderCell},
            {title: 'Submitter (ClinVar)', prop: 'Submitter_ClinVar', render: renderCell},
            {title: 'URL (ENIGMA)', prop: 'URL_ENIGMA', render: renderCell},
            {title: 'Source(s)', prop: 'Source', render: renderCell},
            {title: 'Source URL(s)', prop: 'Source_URL', render: renderCell}
        ]
    },
];

/*eslint-enable camelcase */

// Work-around to allow the user to select text in the table. The browser does not distinguish between
// click and drag: if mouseup and mousedown occur on the same element, a click event is fired even if
// the events occur at very different locations. That makes it hard to select text. This workaround
// defeats the click event if text has been selected.
//
// XXX getSelection().isCollapsed is not available on all platforms. On those platforms we
// will always return false (no selection), so the row click will fire. This makes it hard
// for the user to select text in the table. A better solution would be to add a polyfill for
// getSelection and isCollapsed. There are a few available, though they are much larger than
// what we require:
// https://github.com/Modernizr/Modernizr/wiki/HTML5-Cross-Browser-Polyfills#dom-range-and-selection
// We might want to write a minimal isCollapsed that will use whichever DOM method is available. We could
// also add feature detection, and modify the UI if the feature is not available. E.g. we could style the
// text areas to look clickable instead of selectable..
var hasSelection = () => !(window.getSelection && window.getSelection().isCollapsed);

var Table = React.createClass({
    mixins: [PureRenderMixin],
    render: function () {
        // Expert portal always shows all sources and default columns
        var {data, onHeaderClick, onRowClick, hiddenSources, mode, columnSelection, sourceSelection, ...opts} = this.props;
        if (mode === "default") {
            // Always show default columns/sources in expert mode.
            columnSelection = _.object(_.map(columns,
                            c => _.contains(getDefaultExpertColumns(), c.prop) ? [c.prop, true] : [c.prop, false])
                        );
            sourceSelection = getAllSources();
        }
        return (
            <DataTable
                ref='table'
                className='row-clickable data-table'
                {...opts}
                columnSelection={columnSelection}
                sourceSelection={sourceSelection}
                buildRowOptions={r => ({className: r['Change_Type_id'] === 2 ? 'warning data-table-row' : 'data-table-row', title: 'click for details', onClick: () => hasSelection() ? null : onRowClick(r)})}
                buildHeader={title => buildHeader(onHeaderClick, title)}
                onRowClick={onRowClick}
                onHeaderClick={onHeaderClick}
                filterColumns={filterColumns}
                initialData={data}
                initialPageLength={20}
                initialSortBy={{prop: 'Gene_Symbol', order: 'descending'}}
                pageLengthOptions={[ 20, 50, 100 ]}
                mode={mode}/>
        );
    }
});

var ResearchVariantTableSupplier = function (Component) {
    var ResearchVariantTableComponent = React.createClass({
        mixins: [State, PureRenderMixin],

        getInitialState: function () {
            /*
            Selections take the following order of priority:
                1. Query params in URL
                2. Local Storage
                3. Default settings

            To accomplish this, first set all column selections to true and all filters off.
            Then, update selections according to query params if present. If no query params
            are present, update according to local storage if present. If neither query params
            nor local storage specify changes, use default settings.
            */

            // Start with all columns set to true and data showing from all sources.
            var selectedColumns = selectedColumns = _.object(_.map(this.getColumns(),
                                    c => [c.prop, true])
                                );
            var selectedSources = getAllSources();

            // Get query params.
            const urlParams = this.getQuery();
            const useQueryParams = urlParams.hasOwnProperty("hide") || urlParams.hasOwnProperty("hideSources");

            if (useQueryParams) {
                // If query params are present, use them for settings.
                if (urlParams.hasOwnProperty("hide")) {
                    const columnsToHide = urlParams.hide;
                    for (let i = 0; i < columnsToHide.length; i++) {
                        selectedColumns[columnsToHide[i]] = false;
                    }
                }
                if (urlParams.hasOwnProperty("hideSources")) {
                    const sourcesToHide = urlParams.hideSources;
                    for (let i = 0; i < sourcesToHide.length; i++) {
                        selectedSources[sourcesToHide[i]] = 0;
                    }
                }
            } else {
                // If no query params are present, check local storage.
                const lsSelectedColumns = JSON.parse(localStorage.getItem('columnSelection'));
                if (lsSelectedColumns !== null && lsSelectedColumns !== undefined) {
                    selectedColumns = lsSelectedColumns;
                } else {
                    // If no query params and no local storage, use default settings.
                    selectedColumns = this.getDefaultColumnSelections();
                }
                const lsSelectedSources = JSON.parse(localStorage.getItem('sourceSelection'));
                if (lsSelectedSources !== null && lsSelectedSources !== undefined) {
                    selectedSources = lsSelectedSources;
                }
            }

            return {
                sourceSelection: selectedSources,
                columnSelection: selectedColumns
            };
        },
        toggleColumns: function (prop) {
            let {columnSelection} = this.state,
                val = columnSelection[prop],
                cs = {...columnSelection, [prop]: !val};
            localStorage.setItem('columnSelection', JSON.stringify(cs));
            this.setState({columnSelection: cs});
        },
        setSource: function (prop, event) {
            // this function uses 1, 0 and -1 to accommodate excluding sources as well as not-including them
            // currently only uses 1 and 0 because exclusion is not being used
            let {sourceSelection} = this.state;
            let value = event.target.checked ? 1 : 0;
            let ss = {...sourceSelection, [prop]: value};
            localStorage.setItem('sourceSelection', JSON.stringify(ss));
            this.setState({sourceSelection: ss});
        },
        filterFormCols: function (subColList, columnSelection) {
            return _.map(subColList, ({title, prop}) =>
                <ColumnCheckbox onChange={() => this.toggleColumns(prop)} key={prop} label={prop} title={title}
                                initialCheck={columnSelection}/>);
        },
        getColumnSelectors() {
            var filterFormSubCols = _.map(subColumns, ({subColTitle, subColList}) =>
                <Col sm={6} md={4} key={subColTitle}>
                    <Panel header={subColTitle}>
                        {this.filterFormCols(subColList, this.state.columnSelection)}
                    </Panel>
                </Col>
            );
            return (<label className='control-label'>
                <Panel header="Column Selection">
                    {filterFormSubCols}
                </Panel>
            </label>);
        },
        getFilters: function() {
            var sourceCheckboxes = _.map(this.state.sourceSelection, (value, name) =>
                <Col sm={6} md={3} key={name}>
                    <Input type="checkbox"
                        onChange={v => this.setSource(name, v)}
                        label={name.substring(11).replace(/_/g, " ")} // eg "Variant_in_1000_Genomes" => "1000 Genomes"
                        checked={value > 0}/>
                </Col>
            );
            return (<label className='control-label'>
                <Panel className="top-buffer" header="Source Selection">
                    {sourceCheckboxes}
                </Panel>
            </label>);
        },
        getDownloadButton: function (callback) {
            return <Button className="btn-sm rgt-buffer" download="variants.tsv" href={callback()}>Download</Button>;
        },
        getLollipopButton: function (callback, isOpen) {
            return (<Button id="lollipop-chart-toggle" className="btn-sm rgt-buffer"
                           onClick={callback}>{(isOpen ? 'Hide' : 'Show' ) + ' Lollipop Chart'}</Button>);
        },
        getColumns: function () {
            return researchModeColumns;
        },
        getDefaultColumnSelections: function() {
            return _.object(_.map(researchModeColumns,
                c => _.contains(getDefaultResearchColumns(), c.prop) ? [c.prop, true] : [c.prop, false])
            );
        },
        getDefaultSourceSelections: function() {
            return getAllSources();
        },
        researchVariantTableRestoreDefaults: function(callback) {
            const columnSelection = this.getDefaultColumnSelections();
            const sourceSelection = this.getDefaultSourceSelections();
            this.setState({columnSelection: columnSelection,
                           sourceSelection: sourceSelection},
                           function() {
                                this.props.restoreDefaults(callback);
                           });
        },
        render: function () {
            return (
                <Component
                    {...this.props}
                    researchVariantTableRestoreDefaults={this.researchVariantTableRestoreDefaults}
                    columns={this.getColumns()}
                    columnSelectors={this.getColumnSelectors()}
                    filters={this.getFilters()}
                    sourceSelection={this.state.sourceSelection}
                    columnSelection={this.state.columnSelection}
                    downloadButton={this.getDownloadButton}
                    lollipopButton={this.getLollipopButton}/>
            );
        }
    });
    return ResearchVariantTableComponent;
};

var VariantTableSupplier = function (Component) {
    var VariantTableComponent = React.createClass({
        mixins: [PureRenderMixin],
        getColumns: function () {
            return columns;
        },
        expertVariantTableRestoreDefaults: function(callback) {
            this.props.restoreDefaults(callback);
        },
        render: function () {
            let expertColumns = _.object(_.map(this.getColumns(),
                c => _.contains(getDefaultExpertColumns(), c.prop) ? [c.prop, true] : [c.prop, false])
            );
            // Expert portal always shows all sources
            let sourceSelection = getAllSources();
            return (
                <Component
                    {...this.props}
                    columns={this.getColumns()}
                    columnSelection={expertColumns}
                    sourceSelection={sourceSelection}
                    expertVariantTableRestoreDefaults={this.expertVariantTableRestoreDefaults}
                    downloadButton={()=> null}
                    lollipopButton={()=> null}/>
            );
        }
    });
    return VariantTableComponent;
};


module.exports = {
	VariantTable: VariantTableSupplier(Table),
	ResearchVariantTable: ResearchVariantTableSupplier(Table),
	researchModeColumns: researchModeColumns,
	columns: columns,
    researchModeGroups: researchModeGroups,
    expertModeGroups: expertModeGroups
};
