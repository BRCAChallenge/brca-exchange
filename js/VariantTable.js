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

//function renderClinVarLink(val) {
//    return (
//        <a title="View on ClinVar"
//            onClick={ev => ev.stopPropagation()}
//            href={"http://www.ncbi.nlm.nih.gov/clinvar/?term=" + val}>{val}</a>
//    );
//}

function renderCell(val) {
    return <span>{val}</span>;
}

var filterColumns = [
    {name: 'Gene', prop: 'Gene_Symbol', values: ['BRCA1', 'BRCA2']},
//    {name: 'Exon', values: ['Any', 1, 2, 3, 4, 5]}, // XXX needs refgene to get exon count
    {name: 'Pathogenicity', prop: 'Pathogenicity_default', values: ['Pathogenic', 'Benign / Little Clinical Significance', 'Not Yet Classified']}
];

// XXX duplicate this functionality on the server, perhaps
// by having the client pass in order_by of Genomic_Coordinate
// for hgvs columns.
//var strPropCmpFn = prop => (a, b) => {
//    var ap = a[prop],
//        bp = b[prop];
//    if (ap == null && bp == null || ap === bp) {
//        return 0;
//    }
//    if (bp == null || bp < ap) {
//        return 1;
//    }
//    return -1;
//};
//
//var posCmpFn = strPropCmpFn('Genomic_Coordinate');
//
//function sortColumns(columns, {prop, order}, data) {
//    var sortFn = _.findWhere(columns, {prop: prop}).sortFn || strPropCmpFn(prop),
//        sorted = data.slice(0).sort(sortFn);
//    if (order === 'descending') {
//        sorted.reverse();
//    }
//    return sorted;
//}

var columns = [
    {title: 'Gene', prop: 'Gene_Symbol'},
    {title: 'Genomic (GRCh38)', prop: 'Genomic_Coordinate_hg38'},
    {title: 'Nucleotide', prop: 'HGVS_cDNA'},
    {title: 'Protein', prop: 'HGVS_Protein'},
    {title: 'Pathogenicity', prop: 'Pathogenicity_default'},
    {title: 'Source URL(s)', prop: 'Source_URL'},
];

var research_mode_columns = [

    {title: 'Gene Symbol', prop: 'Gene_Symbol'},
    {title: 'Genome (GRCh36)', prop: 'Genomic_Coordinate_hg36'},
    {title: 'Genome (GRCh37)', prop: 'Genomic_Coordinate_hg37'},
    {title: 'Genome (GRCh38)', prop: 'Genomic_Coordinate_hg38'},

    {title: 'Mutation category (BIC)', prop: 'Mutation_type_BIC'},
    {title: 'PolyPhen score', prop: 'PolyPhen_VEP'},
    {title: 'SIFT score', prop: 'SIFT_VEP'},


    {title: 'BIC Variant Identifier', prop: 'BIC_Identifier'},
    {title: 'Nucleotide', prop: 'HGVS_cDNA'},
    {title: 'Protein', prop: 'HGVS_Protein'},
    {title: 'Analysis Method (ClinVar)', prop: 'SCV_ClinVar'},
    {title: 'Source(s)', prop: 'Source'},
    {title: 'Source URL(s)', prop: 'Source_URL'},
    {title: 'Synonyms', prop: 'Synonyms'},
    {title: 'Protein Amino Acid Change', prop: 'Protein_Change'},
    {title: 'Reference cDNA Sequence', prop: 'Reference_Sequence'},

    {title: 'SCV Accession (ClinVar)', prop: 'Allele_Origin_ClinVar'},
    {title: 'Allele Origin (ENIGMA)', prop: 'Allele_origin_ENIGMA'},
    {title: 'Ethnicity (BIC)', prop: 'Ethnicity_BIC'},
    {title: 'Allele Origin (BIC)', prop: 'Germline_or_Somatic_BIC'},
    {title: 'Allele Origin (LOVD)', prop: 'Origin_of_variant_LOVD'},
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
    {title: 'Pathogenicity', prop: 'Pathogenicity_research'},

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
    {title: 'Allele Origin (ClinVar)', prop: 'Method_ClinVar'},

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
            {title: 'BIC Variant Identifier', prop: 'BIC_Identifier', render: renderCell},
            {title: 'Protein', prop: 'HGVS_Protein', render: renderCell},
            {title: 'Analysis Method (ClinVar)', prop: 'SCV_ClinVar', render: renderCell},
            {title: 'HGVS Nucleotide', prop: 'HGVS_cDNA', render: renderCell},
            {title: 'Protein Amino Acid Change', prop: 'Protein_Change', render: renderCell},
            {title: 'Reference cDNA Sequence', prop: 'Reference_Sequence', render: renderCell}
        ]
    },
    {
        subColTitle: "Origin",
        subColList: [
            {title: 'SCV Accession (ClinVar)', prop: 'Allele_Origin_ClinVar', render: renderCell},
            {title: 'Allele Origin (ENIGMA)', prop: 'Allele_origin_ENIGMA', render: renderCell},
            {title: 'Ethnicity (BIC)', prop: 'Ethnicity_BIC', render: renderCell},
            {title: 'Allele Origin (BIC)', prop: 'Germline_or_Somatic_BIC', render: renderCell},
            {title: 'Allele Origin (LOVD)', prop: 'Origin_of_variant_LOVD', render: renderCell},
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
            {title: 'Maximum Allele Frequency', prop: 'Max_Allele_Frequency',render: renderCell},
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
            {title: 'PolyPhen score', prop: 'PolyPhen_VEP', render: renderCell},
            {title: 'SIFT score', prop: 'SIFT_VEP', render: renderCell}
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
            {title: 'Pathogenicity', prop: 'Pathogenicity_research', render: renderCell},
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
            {title: 'Allele Origin (ClinVar)', prop: 'Method_ClinVar', render: renderCell}
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

var defaultColumns = ['Gene_Symbol', 'Genomic_Coordinate_hg38', 'HGVS_cDNA', 'HGVS_Protein', 'Pathogenicity_default'];
var defaultResearchColumns = ['Gene_Symbol', 'Genomic_Coordinate_hg38', 'HGVS_cDNA', 'HGVS_Protein', 'Pathogenicity_research', 'Allele_Frequency'];

var allSources = {
    Variant_in_ENIGMA: 1,
    Variant_in_ClinVar: 1,
    Variant_in_1000_Genomes: 1,
    Variant_in_ExAC: 1,
    Variant_in_LOVD: 1,
    Variant_in_BIC: 1,
    Variant_in_ESP: 1,
    Variant_in_exLOVD: 1
};
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
        var {data, onHeaderClick, onRowClick, hiddenSources,...opts} = this.props;
        return (
            <DataTable
                ref='table'
                className='row-clickable'
                {...opts}
                buildRowOptions={r => ({title: 'click for details', onClick: () => hasSelection() ? null : onRowClick(r)})}
                buildHeader={title => buildHeader(onHeaderClick, title)}
                onRowClick={onRowClick}
                onHeaderClick={onHeaderClick}
                filterColumns={filterColumns}
                initialData={data}
                initialPageLength={20}
                initialSortBy={{prop: 'Gene_Symbol', order: 'descending'}}
                pageLengthOptions={[ 20, 50, 100 ]}/>
        );
    }
});

var ResearchVariantTableSupplier = function (Component) {
    var ResearchVariantTableComponent = React.createClass({
        mixins: [PureRenderMixin],

        getInitialState: function () {
            var defaultColumnSelection = _.object(
                _.map(this.getColumns(),
                    c => _.contains(this.getDefaultColumns(), c.prop) ? [c.prop, true] : [c.prop, false]));
            var columnSelectionQueryParams = this.props.initialState.columnSelection;

            return {
                sourceSelection: {...allSources, ...this.props.sourceSelection},
                columnSelection: {...defaultColumnSelection, ...columnSelectionQueryParams}
            };
        },
        toggleColumns: function (prop) {
            var {columnSelection} = this.state,
                val = columnSelection[prop],
                cs = {...columnSelection, [prop]: !val};
            this.setState({columnSelection: cs});
        },
        setSource: function (prop, event) {
            // this function uses 1, 0 and -1 to accommodate excluding sources as well as not-including them
            // currently only uses 1 and 0 because exclusion is not being used
            var {sourceSelection} = this.state
            var value = event.target.checked ? 1 : 0;
            var ss = {...sourceSelection, [prop]: value};
            this.setState({sourceSelection: ss});
        },
        filterFormCols: function (subColList, columnSelection) {
            return _.map(subColList, ({title, prop}) =>
                <ColumnCheckbox onChange={v => this.toggleColumns(prop)} key={prop} label={prop} title={title}
                                initialCheck={columnSelection}/>);
        },
        getAdvancedFilters() {
            var sourceCheckboxes = _.map(this.state.sourceSelection, (value, name) =>
                <Col sm={6} md={3} key={name}>
                    <Input type="checkbox"
                        onChange={v => this.setSource(name,v)}
                        label={name.substring(11).replace(/_/g," ")} // eg "Variant_in_1000_Genomes" => "1000 Genomes"
                        checked={value>0}>
                    </Input>
                </Col>
            );
            var filterFormSubCols = _.map(subColumns, ({subColTitle, subColList}) =>
                <Col sm={6} md={4} key={subColTitle}>
                    <Panel header={subColTitle}>
                        {this.filterFormCols(subColList, this.state.columnSelection)}
                    </Panel>
                </Col>
            );
            return (<label className='control-label'>
                <Panel className="top-buffer" header="Source Selection">
                    {sourceCheckboxes}
                </Panel>
                <Panel header="Column Selection">
                    {filterFormSubCols}
                </Panel>
            </label>);
        },
        getDownloadButton: function (callback) {
            return <Button className="btn-sm rgt-buffer" download="variants.csv" href={callback()}>Download</Button>;
        },
        getLollipopButton: function (callback, isOpen) {
            return <Button className="btn-sm rgt-buffer"
                           onClick={callback}>{(isOpen ? 'Hide' : 'Show' ) + ' Lollipop Chart'}</Button>
        },

        getColumns: function () {
            return research_mode_columns;
        },
        getDefaultColumns: function () {
            return defaultResearchColumns;
        },
        render: function () {
            var sourceSelection = this.state.sourceSelection;
            var columnSelection = this.state.columnSelection;
            return (
                <Component
                    {...this.props}
                    columns={this.getColumns()}
                    advancedFilters={this.getAdvancedFilters()}
                    sourceSelection={sourceSelection}
                    columnSelection={columnSelection}
                    downloadButton={this.getDownloadButton}
                    lollipopButton={this.getLollipopButton}/>
            );
        }
    });
    return ResearchVariantTableComponent;
};

var VariantTableSupplier = function (Component) {
    var ResearchVariantTableComponent = React.createClass({
        mixins: [PureRenderMixin],
        getColumns: function () {
            return columns;
        },
        getDefaultColumns: function () {
            return defaultColumns;
        },
        render: function () {
            var columnSelection = _.object(
                _.map(this.getColumns(),
                    c => _.contains(this.getDefaultColumns(), c.prop) ? [c.prop, true] : [c.prop, false]));
            var sourceSelection = allSources
            return (
                <Component
                    {...this.props}
                    columns={this.getColumns()}
                    columnSelection={columnSelection}
                    sourceSelection={sourceSelection}
                    downloadButton={()=> null}
                    lollipopButton={()=> null}/>
            );
        }
    });
    return ResearchVariantTableComponent;
};


module.exports = ({
    VariantTable: VariantTableSupplier(Table),
    ResearchVariantTable: ResearchVariantTableSupplier(Table),
    research_mode_columns: research_mode_columns,
    columns: columns
});
