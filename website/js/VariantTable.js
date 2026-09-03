'use strict';

var {computeReviewStatusScore} = require("./components/VariantSubmitter");

var React = require('react');
import DataTable from './DataTable';
var _ = require('underscore');

// React-Bootstrap v2 replacements
import { Card, Col, Row, Collapse, Form, Button } from 'react-bootstrap';

import ColumnCheckbox from './ColumnCheckbox';
var {getDefaultExpertColumns, getDefaultResearchColumns, getAllSources} = require('./VariantTableDefaults');
var alleleFrequencyCharts = require('./AlleleFrequencyCharts');

// TODO: Re-enable or replace after upgrading react-data-components-brcaex
// require('react-data-components-brcaex/css/table-twbs.css');

function buildHeader(onClick, title) {
    return (
        <span>
            {title}
            <span
                onClick={ev => {ev.stopPropagation(); onClick(title); }}
                className='help fa fa-question-circle'
                style={{cursor: 'pointer',
			fontSize: '0.8em',
			color: '#888',
			position: 'relative',
			top: '-8px',
			left: '4px'
		}}
            />
        </span>
    );
}

function renderCell(val) {
    return <span>{val}</span>;
}

const filterColumns = [
    {name: 'Gene', prop: 'Gene_Symbol', values: ['BRCA1', 'BRCA2']},
    {name: 'Pathogenicity', prop: 'Pathogenicity_expert', values: ['Pathogenic', 'Likely Pathogenic', 'Benign / Little Clinical Significance', 'Likely Benign', 'Not Yet Reviewed']}
];

const expertModeGroups = [
    {groupTitle: 'Alternate Names', internalGroupName: 'Variant Nomenclature', subtitle: 'Below are some of the aliases for this gene change.', innerCols: [
        {title: 'Gene', prop: 'Gene_Symbol', render: gene => <i>{gene}</i>},
        {title: 'Coding Change (c.)', prop: 'HGVS_cDNA', render: nucleotide => nucleotide.split(':')[1]},
        // {title: 'Transcript Identifier', prop: 'Reference_Sequence'},
        // {title: 'HGVS RNA', prop: 'HGVS_RNA'},
        {title: 'Protein Change (p.)', prop: 'HGVS_Protein', render: protein => protein.split(':')[1]},
        // Protein Identfifier is pulled from HGVS_Protein, this is handled in VariantDetail (index.js)
        //{title: 'Protein Identifier', prop: 'HGVS_Protein_ID'},
        {title: 'Abbreviated Protein Change', prop: 'Protein_Change'}, // this is manually renamed to 'Abbreviated AA Change' in the front-end
        // {title: 'BIC Designation', prop: 'BIC_Nomenclature'},
        // {title: 'ClinGen Allele Registry', prop: 'CA_ID'},
        {title: 'Genomic Change (g.)', prop: 'Genomic_Coordinate_hg38'}
        // {title: 'Genomic Nomenclature (GRCh37)', prop: 'Genomic_Coordinate_hg37'}
    ]},

    {groupTitle: 'Overall Assessment of Clinical Significance', internalGroupName: 'Significance (ENIGMA)', subtitle: 'This assessment was made either by the ClinGen Variant Curation Expert Panel, the internationally-recognized expert body for genetic variation in this gene, or by a consensus of the laboratories that share their assessments with ClinVar.', innerCols: [
        {title: 'Clinical Significance', prop: 'Pathogenicity_expert'},
        // {title: 'IARC Class', prop: 'Clinical_significance_ENIGMA'},
        {title: 'Comment on Clinical Significance', prop: 'Comment_on_clinical_significance_ENIGMA'},
        // {title: 'Clinical Significance Citations', prop: 'Clinical_significance_citations_ENIGMA'},
        // {title: 'Supporting Evidence URL(s)', prop: 'URL_ENIGMA'},
        {title: 'Date Last Evaluated', prop: 'Date_last_evaluated_ENIGMA'},
        {title: 'Evaluation Method', prop: 'Assertion_method_ENIGMA'}
        // {title: 'Assertion Method Citation', prop: 'Assertion_method_citation_ENIGMA'},
        // {title: 'Allele Origin', prop: 'Allele_origin_ENIGMA'},
        // {title: 'ClinVar Accession', prop: 'ClinVarAccession_ENIGMA'}
    ]},

    {groupTitle: 'Clinical Significance (ClinVar)', internalGroupName: 'Significance (ClinVar)', subtitle: 'Clinical significance assessments shared with ClinVar.', reportSource: 'ClinVar',
        reportBinding: {
            sortBy: (a, b) => {
                // sort first by stars, then by date in reverse chronological order
                const starDiff = (
                    computeReviewStatusScore(b['Review_Status_ClinVar']) -
                    computeReviewStatusScore(a['Review_Status_ClinVar'])
                );

                const datetimeDiff = (
                    new Date(b['Date_Last_Updated_ClinVar']).getTime() -
                    new Date(a['Date_Last_Updated_ClinVar']).getTime()
                );

                return starDiff || datetimeDiff;
            },
            submitter: {title: 'Submitter', prop: 'Submitter_ClinVar'},
            helpKey: 'clinical-significance-clinvar',
            cols: [
                {title: 'Submitter', prop: 'Submitter_ClinVar'},
                {title: 'Clinical Significance', prop: 'Clinical_Significance_ClinVar'},
                {title: 'Date Significance Last Evaluated', prop: 'DateSignificanceLastEvaluated_ClinVar'},
                {title: 'Date Submission Last Updated', prop: 'Date_Last_Updated_ClinVar'},
                {title: 'Submission Type', prop: 'Method_ClinVar'},
                {title: 'SCV Accession', prop: 'SCV_ClinVar'},
                {title: 'Summary Evidence', prop: 'Summary_Evidence_ClinVar', dummy: true},
                {title: 'Supporting Observations', prop: 'Description_ClinVar', dummy: true},
                {title: 'Review Status', prop: 'Review_Status_ClinVar', dummy: true},
            ]
        }
    },
];

const researchModeGroups = [
    {groupTitle: 'Variant Nomenclature', internalGroupName: 'Variant Nomenclature', innerCols: [
        {title: 'Gene Symbol', prop: 'Gene_Symbol', render: gene => <i>{gene}</i>},
        {title: 'Reference cDNA Sequence', prop: 'Reference_Sequence', core: true},
        {title: 'HGVS Nucleotide', prop: 'HGVS_cDNA', render: nucleotide => nucleotide.split(':')[1], core: true},
        {title: 'HGVS Protein', prop: 'HGVS_Protein', render: protein => protein.split(':')[1], core: true},
        {title: 'Protein Amino Acid Change', prop: 'Protein_Change', core: true},
        {title: 'BIC Designation', prop: 'BIC_Nomenclature', core: true},
        {title: 'ClinGen Allele Registry', prop: 'CA_ID'},
        {title: 'Genome (GRCh38)', prop: 'Genomic_Coordinate_hg38', core: true},
        {title: 'Genome (GRCh37)', prop: 'Genomic_Coordinate_hg37'},
        {title: 'RNA (LOVD)', prop: 'RNA_LOVD'},
        {title: 'Beacons', core: true},
        {title: 'GA4GH VRS Identifier (hg38)', prop: 'VR_ID'},
        {title: 'Synonyms', prop: 'Synonyms', toggleable: true}
    ]},

    {groupTitle: 'Clinical Significance (ENIGMA)', internalGroupName: 'Significance (ENIGMA)', innerCols: [
        {title: 'Clinical Significance', prop: 'Clinical_significance_ENIGMA', core: true},
        {title: 'Comment on Clinical Significance', prop: 'Comment_on_clinical_significance_ENIGMA', core: true},
        {title: 'Assertion Method', prop: 'Assertion_method_ENIGMA', core: true},
        {title: 'Date last evaluated', prop: 'Date_last_evaluated_ENIGMA', core: true},
        {title: 'Collection Method', prop: 'Collection_method_ENIGMA', core: true},
        {title: 'Clinical Significance Citation', prop: 'Clinical_significance_citations_ENIGMA', core: true},
        {title: 'Allele Origin', prop: 'Allele_origin_ENIGMA', core: true},
        // Note: Displayed condition value is actually a combo of Condition_ID_type_ENIGMA and
        // Condition_ID_value_ENIGMA, but only one is used as the prop to ensure proper handling downstream
    ]},

    // if a group contains reportSource, it'll be replaced with a SourceReportsTile component.
    // the reportBinding field must also be specified, and describes how individual reports are bound to nested subtiles
    // (note that specific sources are hardcoded to have different appearances in components/VariantSubmitter.js)
    // (note #2: specifying 'helpKey' on a 'cols' entry will cause the label to render as a help link w/'helpKey' as its target)
    {groupTitle: 'Clinical Significance (ClinVar)', internalGroupName: 'Significance (ClinVar)', reportSource: 'ClinVar',
        reportBinding: {
            sortBy: (a, b) => {
                // sort first by stars, then by date in reverse chronological order
                // (we receive strings for the date, so we have to first interpret them as dates)
                // FIXME: is there a more reliable way to parse these strings as dates?

                const starDiff = (
                    computeReviewStatusScore(b['Review_Status_ClinVar']) -
                    computeReviewStatusScore(a['Review_Status_ClinVar'])
                );

                const datetimeDiff = (
                    new Date(b['Date_Last_Updated_ClinVar']).getTime() -
                    new Date(a['Date_Last_Updated_ClinVar']).getTime()
                );

                return starDiff || datetimeDiff;
            },
            submitter: {title: 'Submitter', prop: 'Submitter_ClinVar'},
            helpKey: 'clinical-significance-clinvar',
            cols: [
                // displayed here in full since the display in the header is potentially truncated
                {title: 'Submitter', prop: 'Submitter_ClinVar'},

                // added these fields (redundantly?) b/c they're not noticeable in the header
                {title: 'Clinical Significance', prop: 'Clinical_Significance_ClinVar'},
                {title: 'Date Significance Last Evaluated', prop:
                'DateSignificanceLastEvaluated_ClinVar'},
                {title: 'Date Submission Last Updated', prop:
                'Date_Last_Updated_ClinVar'},

                {title: 'Submission Type', prop: 'Method_ClinVar'},
                {title: 'SCV Accession', prop: 'SCV_ClinVar'},
                {title: 'Summary Evidence', prop: 'Summary_Evidence_ClinVar', dummy: true},
                {title: 'Supporting Observations', prop: 'Description_ClinVar', dummy: true},
                {title: 'Review Status', prop: 'Review_Status_ClinVar', dummy: true},
                // Note: Displayed condition value is actually a combo of Condition_Value and Condition_DB_ID,
                // but only one is used as the prop to ensure proper handling downstream
            ]
        }
    },
    {groupTitle: 'Clinical Significance (LOVD)', internalGroupName: 'Significance (LOVD)', reportSource: 'LOVD',
        reportBinding: {
            submitter: {title: 'Submitter(s)', prop: 'Submitters_LOVD'},
            helpKey: 'clinical-significance-lovd',
            cols: [
                // displayed here in full since the display in the header is potentially truncated
                {title: 'Submitter(s)', prop: 'Submitters_LOVD'},
                {title: 'Clinical Classification', prop: 'Classification_LOVD'},

                {title: 'Variant Data Type', prop: 'Genetic_origin_LOVD'},
                {title: 'Variant Frequency', prop: 'Variant_frequency_LOVD'},
                {title: 'Individuals', prop: 'Individuals_LOVD'},
                {title: 'Variant ID', prop: 'DBID_LOVD'},
                {title: 'Variant Haplotype', prop: 'Variant_haplotype_LOVD'},
                {title: 'Created Date', prop: 'Created_date_LOVD'},
                {title: 'Edited Date', prop: 'Edited_date_LOVD'},
                {title: 'Variant Remarks', prop: 'Remarks_LOVD'},
            ]
        }
    },

    {groupTitle: 'Computational Predictions', internalGroupName: 'Computational Predictions', innerGroups: [
        {
            source: "BayesDel",
            data: [
                {title: 'Result BayesDel', prop: 'BayesDel_nsfp33a_noAF', core: true},
            ]
        },
        {
            source: "SpliceAI",
            data: [
                {title: 'Result SpliceAI', prop: 'result_spliceai', core: true},
                {title: 'Delta Score Acceptor Gain SpliceAI', prop: 'DS_AG_spliceAI', core: true},
                {title: 'Delta Score Acceptor Loss SpliceAI', prop: 'DS_AL_spliceAI', core: true},
                {title: 'Delta Score Donor Gain SpliceAI', prop: 'DS_DG_spliceAI', core: true},
                {title: 'Delta Score Donor Loss SpliceAI', prop: 'DS_DL_spliceAI', core: true},
                {title: 'Delta Position Acceptor Gain SpliceAI', prop: 'DP_AG_spliceAI', core: true},
                {title: 'Delta Position Acceptor Loss SpliceAI', prop: 'DP_AL_spliceAI', core: true},
                {title: 'Delta Position Donor Gain SpliceAI', prop: 'DP_DG_spliceAI', core: true},
                {title: 'Delta Position Donor Loss SpliceAI', prop: 'DP_DL_spliceAI', core: true}
            ]
        },
    ]},

    {groupTitle: 'Allele Frequency Reference Sets', internalGroupName: 'Allele Frequency Reference Sets', alleleFrequencies: true,
        innerGroups: [
            {
                source: "GnomADv3 Genomes",
                chart: [
                    {title: 'Allele Frequency Charts (gnomAD V3.1 Genomes)', prop: 'Allele_Frequency_Charts_Genome_GnomADv3', replace: alleleFrequencyCharts, tableKey: false, dummy: true},
                ],
                data: [
                    {title: 'Total Allele frequency (gnomAD V3.1 Genomes)', prop: 'Allele_frequency_genome_GnomADv3', core: true},
                    {title: 'Popmax Filtering Allele Frequency (95% confidence) (gnomAD V3.1 Genomes)', prop: 'faf95_popmax_genome_GnomADv3', core: true},
                    {title: 'Allele frequency African (gnomAD V3.1 Genomes)', prop: 'Allele_frequency_genome_AFR_GnomADv3', core: true},
                    {title: 'Allele frequency Latino (gnomAD V3.1 Genomes)', prop: 'Allele_frequency_genome_AMR_GnomADv3', core: true},
                    {title: 'Allele frequency Ashkenazi Jewish (gnomAD V3.1 Genomes)', prop: 'Allele_frequency_genome_ASJ_GnomADv3', core: true},
                    {title: 'Allele frequency East Asian (gnomAD V3.1 Genomes)', prop: 'Allele_frequency_genome_EAS_GnomADv3', core: true},
                    {title: 'Allele frequency Finnish (gnomAD V3.1 Genomes)', prop: 'Allele_frequency_genome_FIN_GnomADv3', core: true},
                    {title: 'Allele frequency Non-Finnish European (gnomAD V3.1 Genomes)', prop: 'Allele_frequency_genome_NFE_GnomADv3', core: true},
                    {title: 'Allele frequency Other (gnomAD V3.1 Genomes)', prop: 'Allele_frequency_genome_OTH_GnomADv3', core: true},
                    {title: 'Allele frequency South Asian (gnomAD V3.1 Genomes)', prop: 'Allele_frequency_genome_SAS_GnomADv3', core: true},
                    {title: 'Allele frequency Middle Eastern (gnomAD V3.1 Genomes)', prop: 'Allele_frequency_genome_MID_GnomADv3', core: true},
                    {title: 'Allele frequency Amish (gnomAD V3.1 Genomes)', prop: 'Allele_frequency_genome_AMI_GnomADv3', core: true}
                ]
            },
            {
                source: "GnomAD Exomes",
                chart: [
                    {title: 'Allele Frequency Charts (gnomAD V2.1 Exomes)', prop: 'Allele_Frequency_Charts_Exome_GnomAD', replace: alleleFrequencyCharts, tableKey: false, dummy: true},
                ],
                data: [
                    {title: 'Total Allele frequency (gnomAD V2.1 Exomes)', prop: 'Allele_frequency_exome_GnomAD', core: true},
                    {title: 'Popmax Filtering Allele Frequency (95% confidence) (gnomAD V2.1 Exomes)', prop: 'faf95_popmax_exome_GnomAD', core: true},
                    {title: 'Allele frequency African (gnomAD V2.1 Exomes)', prop: 'Allele_frequency_exome_AFR_GnomAD', core: true},
                    {title: 'Allele frequency Latino (gnomAD V2.1 Exomes)', prop: 'Allele_frequency_exome_AMR_GnomAD', core: true},
                    {title: 'Allele frequency Ashkenazi Jewish (gnomAD V2.1 Exomes)', prop: 'Allele_frequency_exome_ASJ_GnomAD', core: true},
                    {title: 'Allele frequency East Asian (gnomAD V2.1 Exomes)', prop: 'Allele_frequency_exome_EAS_GnomAD', core: true},
                    {title: 'Allele frequency Finnish (gnomAD V2.1 Exomes)', prop: 'Allele_frequency_exome_FIN_GnomAD', core: true},
                    {title: 'Allele frequency Non-Finnish European (gnomAD V2.1 Exomes)', prop: 'Allele_frequency_exome_NFE_GnomAD', core: true},
                    {title: 'Allele frequency Other (gnomAD V2.1 Exomes)', prop: 'Allele_frequency_exome_OTH_GnomAD', core: true},
                    {title: 'Allele frequency South Asian (gnomAD V2.1 Exomes)', prop: 'Allele_frequency_exome_SAS_GnomAD', core: true},
                ]

            },
        ]
    },

    {groupTitle: 'ACMG Variant Evidence Codes, Provisional Assignment',
        innerGroups: [
            {
                source: "Population Frequency",
                data: [
                    {title: "Provisionally Assigned", prop: "Provisional_Evidence_Code_Popfreq", core: true},
                    {title: "Description", prop: "Provisional_Evidence_Description_Popfreq", core: true}
                ]
            },
            {
                source: "Computational Prediction",
                data: [
                    {title: "Provisionally Assigned", prop: "Provisional_Evidence_Code_Bioinfo", core: true},
                    {title: "Description", prop: "Provisional_Evidence_Description_Bioinfo", core: true}
                ]
            }]
    },


    {groupTitle: 'Multifactorial Likelihood Analysis', internalGroupName: 'Multifactorial Likelihood Analysis', innerCols: [
        {title: 'Posterior probability of pathogenicity', prop: 'Posterior_probability_exLOVD', core: true},
        {title: 'Prior probability of pathogenicity', prop: 'Combined_prior_probablility_exLOVD', core: true},
        {title: 'Missense analysis pathogenicity prior', prop: 'Missense_analysis_prior_probability_exLOVD', core: true},
        {title: 'Co-occurrence likelihood', prop: 'Co_occurrence_LR_exLOVD', core: true},
        {title: 'Segregation Likelihood Ratio', prop: 'Segregation_LR_exLOVD', core: true},
        {title: 'Summary Family History Likelihood Ratio', prop: 'Sum_family_LR_exLOVD', core: true},
        {title: 'Pathology likelihood', prop: 'Pathology_LR_exLOVD', core: true},
        {title: 'Case-Control likelihood', prop: 'Case_control_LR_exLOVD', core: true},
        {title: 'Literature Reference', prop: 'Literature_source_exLOVD', core: true}
    ]},

    {groupTitle: 'In Silico Prediction (prior to considering other evidence)', internalGroupName: 'InSilicoPrediction', inSilicoPred: true,
        hideFromColumnSelection: true,
        innerCols: []
    },

    {groupTitle: 'Functional Assay Results', internalGroupName: 'FunctionalAssayResults',
        innerGroups: [
            {
                source: "ENIGMA BRCA12 Functional Assays",
                data: [
                    {title: 'HGVS Nucleotide', prop: 'HGVS_Nucleotide_ENIGMA_BRCA12_Functional_Assays', core: true},
                    {title: 'Chromosomal Variant', prop: 'Chromosomal_Variant_ENIGMA_BRCA12_Functional_Assays', core: true},
                    {title: 'Report Findlay', prop: 'Result_Findlay_ENIGMA_BRCA12_Functional_Assays', core: true},
                    {title: 'Functional Enrichment Findlay', prop: 'Functional_Enrichment_Findlay_ENIGMA_BRCA12_Functional_Assays', core: true},
                    {title: 'RNA Score Findlay', prop: 'RNA_Score_Findlay_ENIGMA_BRCA12_Functional_Assays', core: true},
                    {title: 'RNA Class Findlay', prop: 'RNA_Class_Findlay_ENIGMA_BRCA12_Functional_Assays', core: true},
                    {title: 'Report Starita', prop: 'Result_Starita_ENIGMA_BRCA12_Functional_Assays', core: true},
                    {title: 'Control Group Petitalot', prop: 'Control_Group_Petitalot_ENIGMA_BRCA12_Functional_Assays', core: true},
                    {title: 'Report Petitalot', prop: 'Result_Petitalot_ENIGMA_BRCA12_Functional_Assays', core: true},
                    {title: 'Selection Bouwman1', prop: 'Selection_Bouwman1_ENIGMA_BRCA12_Functional_Assays', core: true},
                    {title: 'Report Bouwman1', prop: 'Result_Bouwman1_ENIGMA_BRCA12_Functional_Assays', core: true},
                    {title: 'Cisplatin Bouwman2', prop: 'Cisplatin_Bouwman2_ENIGMA_BRCA12_Functional_Assays', core: true},
                    {title: 'Olaparib Bouwman2', prop: 'Olaparib_Bouwman2_ENIGMA_BRCA12_Functional_Assays', core: true},
                    {title: 'DRGFP Bouwman2', prop: 'DRGFP_Bouwman2_ENIGMA_BRCA12_Functional_Assays', core: true},
                    {title: 'Report Bouwman2', prop: 'Result_Bouwman2_ENIGMA_BRCA12_Functional_Assays', core: true},
                    {title: 'Report Fernandes', prop: 'Result_Fernandes_ENIGMA_BRCA12_Functional_Assays', core: true},
                    {title: 'Complementation Mesman', prop: 'Complementation_Mesman_ENIGMA_BRCA12_Functional_Assays', core: true},
                    {title: 'HDR Mesman', prop: 'HDR_Mesman_ENIGMA_BRCA12_Functional_Assays', core: true},
                    {title: 'Cisplatin Mesman', prop: 'Cisplatin_Mesman_ENIGMA_BRCA12_Functional_Assays', core: true},
                    {title: 'Report Mesman', prop: 'Result_Mesman_ENIGMA_BRCA12_Functional_Assays', core: true},
                    {title: 'HDR Richardson', prop: 'HDR_Richardson_ENIGMA_BRCA12_Functional_Assays', core: true},
                    {title: 'Report Richardson', prop: 'Result_Richardson_ENIGMA_BRCA12_Functional_Assays', core: true},
                    {title: 'Report Ikegami', prop: 'Result_Ikegami_ENIGMA_BRCA12_Functional_Assays', core: true},
                    {title: 'Olaparib fClass Ikegami', prop: 'Olaparib_fClass_Ikegami_ENIGMA_BRCA12_Functional_Assays', core: true},
                    {title: 'Niraparif fClass Ikegami', prop: 'Niraparif_fClass_Ikegami_ENIGMA_BRCA12_Functional_Assays', core: true},
                    {title: 'Rucaparib fClass Ikegami', prop: 'Rucaparib_fClass_Ikegami_ENIGMA_BRCA12_Functional_Assays', core: true},
                    {title: 'CBDCA fClass Ikegami', prop: 'CBDCA_fClass_Ikegami_ENIGMA_BRCA12_Functional_Assays', core: true},
                    {title: 'Cell Survival Biwas', prop: 'Cell_Survival_Biwas_ENIGMA_BRCA12_Functional_Assays', core: true},
                    {title: 'Drug Sensitivity Biwas', prop: 'Drug_Sensitivity_Biwas_ENIGMA_BRCA12_Functional_Assays', core: true},
                    {title: 'HAT+DS Score Biwas', prop: 'HAT_DS_Score_Biwas_ENIGMA_BRCA12_Functional_Assays', core: true},
                    {title: 'Report Biwas', prop: 'Result_Biwas_ENIGMA_BRCA12_Functional_Assays', core: true}
                ],
                submitters: [{
                    Biwas: [
                        {title: 'Cell Survival', prop: 'Cell_Survival_Biwas_ENIGMA_BRCA12_Functional_Assays', core: true},
                        {title: 'Drug Sensitivity', prop: 'Drug_Sensitivity_Biwas_ENIGMA_BRCA12_Functional_Assays', core: true},
                        {title: 'HAT + DS Score', prop: 'HAT_DS_Score_Biwas_ENIGMA_BRCA12_Functional_Assays', core: true},
                        {title: 'Report', prop: 'Result_Biwas_ENIGMA_BRCA12_Functional_Assays', core: true},
                    ],
                    Bouwman1: [
                        {title: 'Selection', prop: 'Selection_Bouwman1_ENIGMA_BRCA12_Functional_Assays', core: true},
                        {title: 'Report', prop: 'Result_Bouwman1_ENIGMA_BRCA12_Functional_Assays', core: true},
                    ],
                    Bouwman2: [
                        {title: 'Cisplatin', prop: 'Cisplatin_Bouwman2_ENIGMA_BRCA12_Functional_Assays', core: true},
                        {title: 'Olaparib', prop: 'Olaparib_Bouwman2_ENIGMA_BRCA12_Functional_Assays', core: true},
                        {title: 'DRGFP', prop: 'DRGFP_Bouwman2_ENIGMA_BRCA12_Functional_Assays', core: true},
                        {title: 'Report', prop: 'Result_Bouwman2_ENIGMA_BRCA12_Functional_Assays', core: true},
                    ],
                    Fernandes: [
                        {title: 'Report', prop: 'Result_Fernandes_ENIGMA_BRCA12_Functional_Assays', core: true},
                    ],
                    Findlay: [
                        {title: 'Report', prop: 'Result_Findlay_ENIGMA_BRCA12_Functional_Assays', core: true},
                        {title: 'Functional Enrichment Score', prop: 'Functional_Enrichment_Findlay_ENIGMA_BRCA12_Functional_Assays', core: true},
                        {title: 'RNA Score', prop: 'RNA_Score_Findlay_ENIGMA_BRCA12_Functional_Assays', core: true},
                        {title: 'RNA Class', prop: 'RNA_Class_Findlay_ENIGMA_BRCA12_Functional_Assays', core: true}
                    ],
                    Ikegami: [
                        {title: 'Report', prop: 'Result_Ikegami_ENIGMA_BRCA12_Functional_Assays', core: true},
                    ],
                    Mesman: [
                        {title: 'Complementation', prop: 'Complementation_Mesman_ENIGMA_BRCA12_Functional_Assays', core: true},
                        {title: 'HDR', prop: 'HDR_Mesman_ENIGMA_BRCA12_Functional_Assays', core: true},
                        {title: 'Cisplatin', prop: 'Cisplatin_Mesman_ENIGMA_BRCA12_Functional_Assays', core: true},
                        {title: 'Report', prop: 'Result_Mesman_ENIGMA_BRCA12_Functional_Assays', core: true},
                    ],
                    Petitalot: [
                        {title: 'Control Group', prop: 'Control_Group_Petitalot_ENIGMA_BRCA12_Functional_Assays', core: true},
                        {title: 'Report', prop: 'Result_Petitalot_ENIGMA_BRCA12_Functional_Assays', core: true},
                    ],
                    Richardson: [
                        {title: 'HDR', prop: 'HDR_Richardson_ENIGMA_BRCA12_Functional_Assays', core: true},
                        {title: 'Report', prop: 'Result_Richardson_ENIGMA_BRCA12_Functional_Assays', core: true},
                    ],
                    Starita: [
                        {title: 'Report', prop: 'Result_Starita_ENIGMA_BRCA12_Functional_Assays', core: true},
                    ]
                }]
            }
        ]
    },
];

// subColumns populate the column selection checkboxes.
// They should match the variant detail groupings unless hideFromColumnSelection is true.
const subColumns = _.map(_.filter(researchModeGroups, function(group) { return !group.hideFromColumnSelection; }), function (group) {
    if (Object.prototype.hasOwnProperty.call(group, 'innerCols')) {
        return {
            subColTitle: group.groupTitle,
            // hide dummy columns from column selection
            subColList: _.map(_.filter(group.innerCols, ({dummy}) => !dummy), function (innerCol) {
                return {
                    title: innerCol.title,
                    prop: innerCol.prop,
                    render: renderCell
                };
            })
        };
    } else if (Object.prototype.hasOwnProperty.call(group, 'reportBinding') && Object.prototype.hasOwnProperty.call(group.reportBinding, 'cols')) {
        return {
            subColTitle: group.groupTitle,
            // hide dummy columns from column selection
            subColList: _.map(_.filter(group.reportBinding.cols, ({dummy}) => !dummy), function (col) {
                return {
                    title: col.title,
                    prop: col.prop,
                    render: renderCell
                };
            })
        };
    } else if (Object.prototype.hasOwnProperty.call(group, 'innerGroups')) {
        let subCols = group.innerGroups.map(function(innerGroup) {
            // hide dummy columns from column selection
            return _.map(_.filter(innerGroup.data, ({dummy}) => !dummy), function (col) {
                            return {
                                title: col.title,
                                prop: col.prop,
                                render: renderCell
                            };
                        });
        });
        // flatten array of arrays
        let flattenedSubColList = [].concat.apply([], subCols);
        return {
            subColTitle: group.groupTitle,
            subColList: flattenedSubColList
        };
    }
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
    {title: 'ClinGen Allele Registry', prop: 'CA_ID'},
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
    {title: 'Genome (GRCh37)', prop: 'Genomic_Coordinate_hg37'},
    {title: 'Genome (GRCh38)', prop: 'Genomic_Coordinate_hg38'},
    {title: 'Mutation category (BIC)', prop: 'Mutation_type_BIC'},
    {title: 'BIC Variant Identifier', prop: 'BIC_Nomenclature'},
    {title: 'ClinGen Allele Registry', prop: 'CA_ID'},
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
    {title: 'Variant Data Type (LOVD)', prop: 'Genetic_origin_LOVD'},
    {title: 'Clinical Classification (LOVD)', prop: 'Classification_LOVD'},
    {title: 'Individuals (LOVD)', prop: 'Individuals_LOVD'},
    {title: 'Created Date (LOVD)', prop: 'Created_date_LOVD'},
    {title: 'Edited Date (LOVD)', prop: 'Edited_date_LOVD'},
    {title: 'Variant Remarks (LOVD)', prop: 'Remarks_LOVD'},
    {title: 'Variant ID (LOVD)', prop: 'DBID_LOVD'},
    {title: 'Allele Origin (ClinVar)', prop: 'Allele_Origin_ClinVar'},
    {title: 'Allele Origin (ENIGMA)', prop: 'Allele_origin_ENIGMA'},
    {title: 'Ethnicity (BIC)', prop: 'Ethnicity_BIC'},
    {title: 'Allele Origin (BIC)', prop: 'Germline_or_Somatic_BIC'},
    {title: 'Patient Nationality (BIC)', prop: 'Patient_nationality_BIC'},
    {title: 'Variant Haplotype (LOVD)', prop: 'Variant_haplotype_LOVD'},
    {title: 'Family members carrying this variant (BIC)', prop: 'Number_of_family_member_carrying_mutation_BIC'},
    {title: 'Co-occurrence likelihood (ExUV)', prop: 'Co_occurrence_LR_exLOVD'},
    {title: 'Pathology likelihood (ExUV)', prop: 'Pathology_LR_exLOVD'},
    {title: 'Case-Control likelihood (ExUV)', prop: 'Case_control_LR_exLOVD'},
    {title: 'Prior probability of pathogenicity (ExUV)', prop: 'Combined_prior_probablility_exLOVD'},
    {
        title: 'Missense analysis probability of pathogenicity (ExUV)',
        prop: 'Missense_analysis_prior_probability_exLOVD'
    },
    {title: 'Probability of pathogenicity (ExUV)', prop: 'Posterior_probability_exLOVD'},
    {title: 'Segregation Likelihood Ratio (ExUV)', prop: 'Segregation_LR_exLOVD'},
    {title: 'Summary Family History Likelihood Ratio (ExUV)', prop: 'Sum_family_LR_exLOVD'},
    {title: 'Assertion Method (ENIGMA)', prop: 'Assertion_method_citation_ENIGMA'},
    {title: 'Clinical Significance Citation (ENIGMA)', prop: 'Clinical_significance_citations_ENIGMA'},
    {title: 'Literature Reference (BIC)', prop: 'Literature_citation_BIC'},
    {title: 'Literature Reference (ExUV)', prop: 'Literature_source_exLOVD'},
    {title: 'Pathogenicity', prop: 'Pathogenicity_all'},
    {title: 'Assertion Method (ENIGMA)', prop: 'Assertion_method_ENIGMA'},
    {title: 'Clinical Significance (ClinVar)', prop: 'Clinical_Significance_ClinVar'},
    {title: 'Clinical Significance (ENIGMA)', prop: 'Clinical_significance_ENIGMA'},
    {title: 'Collection Method (ENIGMA)', prop: 'Collection_method_ENIGMA'},
    {title: 'Comment on Clinical Significance (ENIGMA)', prop: 'Comment_on_clinical_significance_ENIGMA'},
    {title: 'Date last evaluated (ENIGMA)', prop: 'Date_last_evaluated_ENIGMA'},
    {title: 'Date last updated (ClinVar)', prop: 'Date_Last_Updated_ClinVar'},
    {title: 'Date Significance Last Evaluated (ClinVar)', prop: 'DateSignificanceLastEvaluated_ClinVar'},
    {title: 'Functional Analysis Result (LOVD)', prop: 'Functional_analysis_result_LOVD'},
    {title: 'Functional Analysis Method (LOVD)', prop: 'Functional_analysis_technique_LOVD'},
    {title: 'Analysis Method (ClinVar)', prop: 'Method_ClinVar'},
    {title: 'Summary Evidence (ClinVar)', prop: 'Summary_Evidence_ClinVar'},
    {title: 'Supporting Observations (ClinVar)', prop: 'Description_ClinVar'},
    {title: 'Review Status (ClinVar)', prop: 'Review_Status_ClinVar'},
    {title: 'ClinVar Accession', prop: 'ClinVarAccession_ENIGMA'},
    {title: 'Condition Category (ENIGMA)', prop: 'Condition_category_ENIGMA'},
    {title: 'Condition ID Type (ENIGMA)', prop: 'Condition_ID_type_ENIGMA'},
    {title: 'Condition ID Value (ENIGMA)', prop: 'Condition_ID_value_ENIGMA'},
    {title: 'Submitter (ClinVar)', prop: 'Submitter_ClinVar'},
    {title: 'URL (ENIGMA)', prop: 'URL_ENIGMA'},
    {title: 'Allele Frequency', prop: 'Allele_Frequency'},
    {title: 'Maximum Allele Frequency', prop: 'Max_Allele_Frequency'},
    {title: 'Allele frequency (gnomAD V3.1 Genomes)', prop: 'Allele_frequency_genome_GnomADv3'},
    {title: 'AFR Allele frequency (gnomAD V3.1 Genomes)', prop: 'Allele_frequency_genome_AFR_GnomADv3'},
    {title: 'AMR Allele frequency (gnomAD V3.1 Genomes)', prop: 'Allele_frequency_genome_AMR_GnomADv3'},
    {title: 'ASJ Allele frequency (gnomAD V3.1 Genomes)', prop: 'Allele_frequency_genome_ASJ_GnomADv3'},
    {title: 'EAS Allele frequency (gnomAD V3.1 Genomes)', prop: 'Allele_frequency_genome_EAS_GnomADv3'},
    {title: 'FIN Allele frequency (gnomAD V3.1 Genomes)', prop: 'Allele_frequency_genome_FIN_GnomADv3'},
    {title: 'NFE Allele frequency (gnomAD V3.1 Genomes)', prop: 'Allele_frequency_genome_NFE_GnomADv3'},
    {title: 'OTH Allele frequency (gnomAD V3.1 Genomes)', prop: 'Allele_frequency_genome_OTH_GnomADv3'},
    {title: 'SAS Allele frequency (gnomAD V3.1 Genomes)', prop: 'Allele_frequency_genome_SAS_GnomADv3'},
    {title: 'MID Allele frequency (gnomAD V3.1 Genomes)', prop: 'Allele_frequency_genome_MID_GnomADv3'},
    {title: 'AMI Allele frequency (gnomAD V3.1 Genomes)', prop: 'Allele_frequency_genome_AMI_GnomADv3'},
    {title: 'Allele frequency (gnomAD V2.1 Exomes)', prop: 'Allele_frequency_exome_GnomAD'},
    {title: 'AFR Allele frequency (gnomAD V2.1 Exomes)', prop: 'Allele_frequency_exome_AFR_GnomAD'},
    {title: 'AMR Allele frequency (gnomAD V2.1 Exomes)', prop: 'Allele_frequency_exome_AMR_GnomAD'},
    {title: 'ASJ Allele frequency (gnomAD V2.1 Exomes)', prop: 'Allele_frequency_exome_ASJ_GnomAD'},
    {title: 'EAS Allele frequency (gnomAD V2.1 Exomes)', prop: 'Allele_frequency_exome_EAS_GnomAD'},
    {title: 'FIN Allele frequency (gnomAD V2.1 Exomes)', prop: 'Allele_frequency_exome_FIN_GnomAD'},
    {title: 'NFE Allele frequency (gnomAD V2.1 Exomes)', prop: 'Allele_frequency_exome_NFE_GnomAD'},
    {title: 'OTH Allele frequency (gnomAD V2.1 Exomes)', prop: 'Allele_frequency_exome_OTH_GnomAD'},
    {title: 'SAS Allele frequency (gnomAD V2.1 Exomes)', prop: 'Allele_frequency_exome_SAS_GnomAD'},
    {title: 'Allele frequency (gnomAD Genomes)', prop: 'Allele_frequency_genome_GnomAD'},
    {title: 'AFR Allele frequency (gnomAD Genomes)', prop: 'Allele_frequency_genome_AFR_GnomAD'},
    {title: 'AMR Allele frequency (gnomAD Genomes)', prop: 'Allele_frequency_genome_AMR_GnomAD'},
    {title: 'ASJ Allele frequency (gnomAD Genomes)', prop: 'Allele_frequency_genome_ASJ_GnomAD'},
    {title: 'EAS Allele frequency (gnomAD Genomes)', prop: 'Allele_frequency_genome_EAS_GnomAD'},
    {title: 'FIN Allele frequency (gnomAD Genomes)', prop: 'Allele_frequency_genome_FIN_GnomAD'},
    {title: 'NFE Allele frequency (gnomAD Genomes)', prop: 'Allele_frequency_genome_NFE_GnomAD'},
    {title: 'OTH Allele frequency (gnomAD Genomes)', prop: 'Allele_frequency_genome_OTH_GnomAD'},
    {title: 'SAS Allele frequency (gnomAD Genomes)', prop: 'Allele_frequency_genome_SAS_GnomAD'},
    {title: 'Allele frequency (gnomAD Exomes)', prop: 'Allele_frequency_exome_GnomAD'},
    {title: 'AFR Allele frequency (gnomAD Exomes)', prop: 'Allele_frequency_exome_AFR_GnomAD'},
    {title: 'AMR Allele frequency (gnomAD Exomes)', prop: 'Allele_frequency_exome_AMR_GnomAD'},
    {title: 'ASJ Allele frequency (gnomAD Exomes)', prop: 'Allele_frequency_exome_ASJ_GnomAD'},
    {title: 'EAS Allele frequency (gnomAD Exomes)', prop: 'Allele_frequency_exome_EAS_GnomAD'},
    {title: 'FIN Allele frequency (gnomAD Exomes)', prop: 'Allele_frequency_exome_FIN_GnomAD'},
    {title: 'NFE Allele frequency (gnomAD Exomes)', prop: 'Allele_frequency_exome_NFE_GnomAD'},
    {title: 'OTH Allele frequency (gnomAD Exomes)', prop: 'Allele_frequency_exome_OTH_GnomAD'},
    {title: 'SAS Allele frequency (gnomAD Exomes)', prop: 'Allele_frequency_exome_SAS_GnomAD'},
    {title: 'Variant Frequency (LOVD)', prop: 'Variant_frequency_LOVD'},
    {title: 'In Silico Prior Probability', prop: 'applicablePrior'},
    {title: 'Protein-level Estimation', prop: 'proteinPrior'},
    {title: 'Donor Impact', prop: 'refDonorPrior'},
    {title: 'De Novo Donor', prop: 'deNovoDonorPrior'},
    {title: 'Acceptor Impact', prop: 'refAccPrior'},
    {title: 'HGVS Nucleotide', prop: 'HGVS_Nucleotide_ENIGMA_BRCA12_Functional_Assays'},
    {title: 'Chromosomal Variant', prop: 'Chromosomal_Variant_ENIGMA_BRCA12_Functional_Assays'},
    {title: 'Report Findlay', prop: 'Result_Findlay_ENIGMA_BRCA12_Functional_Assays'},
    {title: 'Functional Enrichment Findlay', prop: 'Functional_Enrichment_Findlay_ENIGMA_BRCA12_Functional_Assays'},
    {title: 'RNA Score Findlay', prop: 'RNA_Score_Findlay_ENIGMA_BRCA12_Functional_Assays'},
    {title: 'RNA Class Findlay', prop: 'RNA_Class_Findlay_ENIGMA_BRCA12_Functional_Assays'},
    {title: 'Report Starita', prop: 'Result_Starita_ENIGMA_BRCA12_Functional_Assays'},
    {title: 'Control Group Petitalot', prop: 'Control_Group_Petitalot_ENIGMA_BRCA12_Functional_Assays'},
    {title: 'Report Petitalot', prop: 'Result_Petitalot_ENIGMA_BRCA12_Functional_Assays'},
    {title: 'Selection Bouwman1', prop: 'Selection_Bouwman1_ENIGMA_BRCA12_Functional_Assays'},
    {title: 'Report Bouwman1', prop: 'Result_Bouwman1_ENIGMA_BRCA12_Functional_Assays'},
    {title: 'Cisplatin Bouwman2', prop: 'Cisplatin_Bouwman2_ENIGMA_BRCA12_Functional_Assays'},
    {title: 'Olaparib Bouwman2', prop: 'Olaparib_Bouwman2_ENIGMA_BRCA12_Functional_Assays'},
    {title: 'DRGFP Bouwman2', prop: 'DRGFP_Bouwman2_ENIGMA_BRCA12_Functional_Assays'},
    {title: 'Report Bouwman2', prop: 'Result_Bouwman2_ENIGMA_BRCA12_Functional_Assays'},
    {title: 'Report Fernandes', prop: 'Result_Fernandes_ENIGMA_BRCA12_Functional_Assays'},
    {title: 'Complementation Mesman', prop: 'Complementation_Mesman_ENIGMA_BRCA12_Functional_Assays'},
    {title: 'HDR Mesman', prop: 'HDR_Mesman_ENIGMA_BRCA12_Functional_Assays'},
    {title: 'Cisplatin Mesman', prop: 'Cisplatin_Mesman_ENIGMA_BRCA12_Functional_Assays'},
    {title: 'Report Mesman', prop: 'Result_Mesman_ENIGMA_BRCA12_Functional_Assays'},
    {title: 'HDR Richardson', prop: 'HDR_Richardson_ENIGMA_BRCA12_Functional_Assays'},
    {title: 'Report Richardson', prop: 'Result_Richardson_ENIGMA_BRCA12_Functional_Assays'},
    {title: 'Report Ikegami', prop: 'Result_Ikegami_ENIGMA_BRCA12_Functional_Assays'},
    {title: 'Olaparib fClass Ikegami', prop: 'Olaparib_fClass_Ikegami_ENIGMA_BRCA12_Functional_Assays'},
    {title: 'Niraparif fClass Ikegami', prop: 'Niraparif_fClass_Ikegami_ENIGMA_BRCA12_Functional_Assays'},
    {title: 'Rucaparib fClass Ikegami', prop: 'Rucaparib_fClass_Ikegami_ENIGMA_BRCA12_Functional_Assays'},
    {title: 'CBDCA fClass Ikegami', prop: 'CBDCA_fClass_Ikegami_ENIGMA_BRCA12_Functional_Assays'},
    {title: 'Cell Survival Biwas', prop: 'Cell_Survival_Biwas_ENIGMA_BRCA12_Functional_Assays'},
    {title: 'Drug Sensitivity Biwas', prop: 'Drug_Sensitivity_Biwas_ENIGMA_BRCA12_Functional_Assays'},
    {title: 'HAT+DS Score Biwas', prop: 'HAT_DS_Score_Biwas_ENIGMA_BRCA12_Functional_Assays'},
    {title: 'Report Biwas', prop: 'Result_Biwas_ENIGMA_BRCA12_Functional_Assays'},
    {title: 'Result BayesDel', prop: 'BayesDel_nsfp33a_noAF'},
    {title: 'Result SpliceAI', prop: 'result_spliceai'},
    {title: 'Delta Score Acceptor Gain SpliceAI', prop: 'DS_AG_spliceAI'},
    {title: 'Delta Score Acceptor Loss SpliceAI', prop: 'DS_AL_spliceAI'},
    {title: 'Delta Score Donor Gain SpliceAI', prop: 'DS_DG_spliceAI'},
    {title: 'Delta Score Donor Loss SpliceAI', prop: 'DS_DL_spliceAI'},
    {title: 'Delta Position Acceptor Gain SpliceAI', prop: 'DP_AG_spliceAI'},
    {title: 'Delta Position Acceptor Loss SpliceAI', prop: 'DP_AL_spliceAI'},
    {title: 'Delta Position Donor Gain SpliceAI', prop: 'DP_DG_spliceAI'},
    {title: 'Delta Position Donor Loss SpliceAI', prop: 'DP_DL_spliceAI'}
];




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

class Table extends React.Component {
    constructor(props) {
	super(props);
	this.dataTableRef = React.createRef();
    }

    render() {
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
                ref={this.dataTableRef}
                className='row-clickable data-table table-grayheader'
                {...opts}
                columnSelection={columnSelection}
                sourceSelection={sourceSelection}
                buildRowOptions={r => ({
                    className: r['Change_Type_id'] === 2 ? 'warning data-table-row' : 'data-table-row',
                    title: 'click for details',
                    onMouseUp: (e) => hasSelection() ? null : onRowClick(r, e)
                })}
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
}

var ResearchVariantTableSupplier = function (Component) {
    class ResearchVariantTableComponent extends React.Component {
        state = (() => {
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
            var selectedColumns = _.object(_.map(this.getColumns(),
                                    c => [c.prop, true])
                                );
            var selectedSources = getAllSources();

            // Get query params.
            const urlParams = this.getQuery();
            const useQueryParams = Object.prototype.hasOwnProperty.call(urlParams, "hide") || Object.prototype.hasOwnProperty.call(urlParams, "hideSources");

            if (useQueryParams) {
                // If query params are present, use them for settings.
                if (Object.prototype.hasOwnProperty.call(urlParams, "hide")) {
                    const columnsToHide = urlParams.hide;
                    for (let i = 0; i < columnsToHide.length; i++) {
                        selectedColumns[columnsToHide[i]] = false;
                    }
                }
                if (Object.prototype.hasOwnProperty.call(urlParams, "hideSources")) {
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

            // initialize open state for each sub-column panel from localStorage (RB0.33 defaultExpanded logic)
            const openSubcols = {};
            subColumns.forEach(({ subColTitle }) => {
                const key = "collapse-subcol_" + subColTitle;
                openSubcols[subColTitle] = localStorage.getItem(key) !== "true";
            });

            return {
                sourceSelection: selectedSources,
                columnSelection: selectedColumns,
                changeInProgress: false,
                openSubcols
            };
        })();
        componentDidUpdate(prevProps) {
            // Change is now complete (has propagated all the way
            // down and back up through the parent component)
	    if (this.state.changeInProgress && prevProps !== this.props) {
                this.setState({ changeInProgress: false });
            }
        }

        getQuery() {
            const search =
                (this.props.location && this.props.location.search) ||
                (typeof window !== "undefined" ? window.location.search : "") ||
                "";

            const params = new URLSearchParams(search);
            const out = {};

            const readMulti = (key) => {
                const vals = params.getAll(key);
                if (!vals || vals.length === 0) return null;
                return vals
                    .flatMap(v => String(v).split(','))
                    .map(s => s.trim())
                    .filter(Boolean);
            };

            const hide = readMulti("hide");
            if (hide) out.hide = hide;

            const hideSources = readMulti("hideSources");
            if (hideSources) out.hideSources = hideSources;

            return out;
        }

        toggleColumns = (prop) => {
            let {columnSelection} = this.state,
                val = columnSelection[prop],
                cs = {...columnSelection, [prop]: !val};
            localStorage.setItem('columnSelection', JSON.stringify(cs));
            this.setState({columnSelection: cs, changeInProgress: true});
        };

        setSource = (prop, event) => {
            // this function uses 1, 0 and -1 to accommodate excluding sources as well as not-including them
            // currently only uses 1 and 0 because exclusion is not being used
            let {sourceSelection} = this.state;
            let value = event.target.checked ? 1 : 0;
            let ss = {...sourceSelection, [prop]: value};
            localStorage.setItem('sourceSelection', JSON.stringify(ss));
            this.setState({sourceSelection: ss, changeInProgress: true});
        };

        filterFormCols = (subColList, columnSelection) => {
            return _.map(subColList, ({title, prop}) =>
                (<ColumnCheckbox onChange={() => this.toggleColumns(prop)} key={prop || title} label={prop || title} title={title}
                                initialCheck={columnSelection}/>));
        };

        toggleSubcolOpen = (subColTitle) => {
	   const key = "collapse-subcol_" + subColTitle;
	   this.setState((prev) => {
	     const isOpen = !!prev.openSubcols[subColTitle];
	     const nextOpen = !isOpen;
	     // store "true" when collapsed (matches old semantics)
	     localStorage.setItem(key, nextOpen ? "false" : "true");
	     return {
	       openSubcols: { ...prev.openSubcols, [subColTitle]: nextOpen }
	     };
	   });
        };

        getColumnSelectors() {
            var filterFormSubCols = _.map(subColumns, ({subColTitle, subColList}) =>
                (<Col sm={6} md={4} key={subColTitle}>
                    <Card className="mb-3">
                        <Card.Header
                            onClick={() => this.toggleSubcolOpen(subColTitle)}
                            role="button"
                            aria-expanded={!!this.state.openSubcols[subColTitle]}
                            className="d-flex justify-content-between align-items-center"
                        >
                            <span className="fw-bold">{subColTitle}</span>
                            <i className={`fa fa-chevron-${this.state.openSubcols[subColTitle] ? 'down' : 'right'}`} aria-hidden="true" />
                        </Card.Header>
                        <Collapse in={!!this.state.openSubcols[subColTitle]}>
                            <div>
                                <Card.Body>
                                    {this.filterFormCols(subColList, this.state.columnSelection)}
                                </Card.Body>
                            </div>
                        </Collapse>
                    </Card>
                </Col>)
            );
            return (<div className='control-label'>
                <Card>
                    <Card.Header className="fw-bold">Column Selection</Card.Header>
                    <Card.Body>
		    	<Row>
                    	    {filterFormSubCols}
		    	</Row>
                    </Card.Body>
                </Card>
            </div>);
        }

        getSourceName(name) {
            let source = name.substring(11).replace(/_/g, " ");
            if (source.toLowerCase() === "exlovd") {
                source = "ExUV";
            } else if (source.toLowerCase() === "enigma_brca12_functional_assays") {
                source = "ENIGMA BRCA12 Functional Assays";
            } else if (source.toLowerCase() === "gnomad") {
                source = "gnomAD 2.1 Exomes";
            } else if (source.toLowerCase() === "gnomadv3") {
                source = "gnomAD 3.1 Genomes";
            } else if (source.toLowerCase() === "esp") {
                source = "ESP (deprecated)";
            } else if (source.toLowerCase() === "exac") {
                source = "ExAC (deprecated)";
            } else if (source.toLowerCase() === "1000 genomes") {
                source = "1000 Genomes (deprecated)";
            }
            return source;
        }

        getFilters() {
            var sourceCheckboxes = _.map(this.state.sourceSelection, (value, name) =>
                (<Col sm={6} md={3} key={name}>
                    <Form.Check
                        type="checkbox"
                        onChange={v => this.setSource(name, v)}
                        checked={value > 0}
                        label={this.getSourceName(name)}
                    />
                </Col>)
            );
            return (<div className='source-filters'>
                <Card className="top-buffer">
                    <Card.Header className="fw-bold">Source Selection</Card.Header>
                    <Card.Body>
		    	<Row>
                    	    {sourceCheckboxes}
		    	</Row>
                    </Card.Body>
                </Card>
            </div>);
        }

        getDownloadButton(callback) {
            return <Button variant="secondary" className="rgt-buffer" download="variants.tsv" href={callback()}>Download</Button>;
        }

        getColumns() {
            return researchModeColumns;
        }

        getDefaultColumnSelections() {
            return _.object(_.map(researchModeColumns,
                c => _.contains(getDefaultResearchColumns(), c.prop) ? [c.prop, true] : [c.prop, false])
            );
        }

        getDefaultSourceSelections() {
            return getAllSources();
        }

        researchVariantTableRestoreDefaults = (callback) => {
            const columnSelection = this.getDefaultColumnSelections();
            const sourceSelection = this.getDefaultSourceSelections();
            this.setState({columnSelection: columnSelection,
                           sourceSelection: sourceSelection},
                           function() {
                                this.props.restoreDefaults(callback);
                           });
        };

        render() {
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
                    changeInProgress={this.state.changeInProgress}/>
            );
        }
    }
    return ResearchVariantTableComponent;
};

var VariantTableSupplier = function (Component) {
    class VariantTableComponent extends React.PureComponent {
        getColumns() {
            return columns;
        }
        expertVariantTableRestoreDefaults = (callback) => {
            this.props.restoreDefaults(callback);
        };
        render() {
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
                    downloadButton={()=> null}/>
            );
        }
    }
    return VariantTableComponent;
};

const VariantTable = VariantTableSupplier(Table);
const ResearchVariantTable = ResearchVariantTableSupplier(Table);

export {
    VariantTable,
    ResearchVariantTable,
    researchModeColumns,
    columns,
    researchModeGroups,
    expertModeGroups
};
