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
    {name: 'Gene', prop: 'Gene_symbol', values: ['BRCA1', 'BRCA2']},
//    {name: 'Exon', values: ['Any', 1, 2, 3, 4, 5]}, // XXX needs refgene to get exon count
    {name: 'Pathogenicity', prop: 'Clinical_significance', values: ['Pathogenic', 'Benign']}
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
    {title: 'Gene', prop: 'Gene_symbol', render: renderCell},
    {title: 'Genomic Coordinate', prop: 'Genomic_Coordinate', render: renderCell},
    {title: 'HGVS cDNA', prop: 'HGVS_cDNA', /*sortFn: posCmpFn, */render: renderCell},
    {title: 'HGVS protein', prop: 'HGVS_protein', /*sortFn: posCmpFn, */render: renderCell},
    {title: 'HGVS protein (Abbrev.)', prop: "Abbrev_AA_change", render: renderCell},
    {title: 'BIC nucleotide', prop: "BIC_Nomenclature", render: renderCell},
    {title: 'Pathogenicity', prop: 'Clinical_significance', render: renderCell},
    {title: 'Allele frequency (1000 Genomes)', prop: 'Allele_frequency_1000_Genomes', render: renderCell},
    {title: 'SAS Allele frequency (1000 Genomes)', prop: 'SAS_Allele_frequency_1000_Genomes', render: renderCell},
    {title: 'EAS Allele frequency (1000 Genomes)', prop: 'EAS_Allele_frequency_1000_Genomes', render: renderCell},
    {title: 'AMR Allele frequency (1000 Genomes)', prop: 'AMR_Allele_frequency_1000_Genomes', render: renderCell},
    {title: 'EUR Allele frequency (1000 Genomes)', prop: 'EUR_Allele_frequency_1000_Genomes', render: renderCell},
    {title: 'AFR Allele frequency (1000 Genomes)', prop: 'AFR_Allele_frequency_1000_Genomes', render: renderCell},
    {title: 'Allele origin (ClinVar)', prop: 'Allele_origin_ClinVar', render: renderCell},
    {title: 'Variant clinical significance (ClinVar)', prop: 'Variant_clinical_significance_ClinVar', render: renderCell},
    {title: 'HGVS genomic (LOVD)', prop: 'HGVS_genomic_LOVD', render: renderCell},
    {title: 'Origin of variant (LOVD)', prop: 'Origin_of_variant_LOVD', render: renderCell},
    {title: 'HGVS protein (LOVD)', prop: 'HGVS_protein_LOVD', render: renderCell},
    {title: 'Variant frequency (LOVD)', prop: 'Variant_frequency_LOVD', render: renderCell},
    {title: 'HGVS cDNA (LOVD)', prop: 'HGVS_cDNA_LOVD', render: renderCell},
    {title: 'Variant affecting protein (LOVD)', prop: 'Variant_affecting_protein_LOVD', render: renderCell},
    {title: 'Variant haplotype (LOVD)', prop: 'Variant_haplotype_LOVD', render: renderCell},
    {title: 'VEP Gene (ExAC)', prop: 'VEP_Gene_ExAC', render: renderCell},
    {title: 'Allele frequency (ExAC)', prop: 'Allele_frequency_ExAC', render: renderCell},
    {title: 'VEP HGVSc (ExAC)', prop: 'VEP_HGVSc_ExAC', render: renderCell},
    {title: 'VEP Consequence (ExAC)', prop: 'VEP_Consequence_ExAC', render: renderCell},
    {title: 'VEP HGVSp (ExAC)', prop: 'VEP_HGVSp_ExAC', render: renderCell},
    {title: 'Exon number (exLOVD)', prop: 'Exon_number_exLOVD', render: renderCell},
    {title: 'IARC class (exLOVD)', prop: 'IARC_class_exLOVD', render: renderCell},
    {title: 'BIC (exLOVD)', prop: 'BIC_exLOVD', render: renderCell},
    {title: 'HGVS cDNA (exLOVD)', prop: 'HGVS_cDNA_exLOVD', render: renderCell},
    {title: 'Literature source (exLOVD)', prop: 'Literature_source_exLOVD', render: renderCell},
    {title: 'HGVS protein (exLOVD)', prop: 'HGVS_protein_exLOVD', render: renderCell}
];

var subColumns = [
    {subColTitle: "ENIGMA",
     subColList: [
        {title: 'Gene', prop: 'Gene_symbol', render: renderCell},
        {title: 'Genomic Coordinate', prop: 'Genomic_Coordinate', render: renderCell},
        {title: 'HGVS cDNA', prop: 'HGVS_cDNA', /*sortFn: posCmpFn, */render: renderCell},
        {title: 'HGVS protein', prop: 'HGVS_protein', /*sortFn: posCmpFn, */render: renderCell},
        {title: 'HGVS protein (Abbrev.)', prop: "Abbrev_AA_change", render: renderCell},
        {title: 'BIC nucleotide', prop: "BIC_Nomenclature", render: renderCell},
        {title: 'Pathogenicity', prop: 'Clinical_significance', render: renderCell}
     ]
    },
    {subColTitle: "1000 Genomes",
     subColList: [
        {title: 'Allele frequency', prop: 'Allele_frequency_1000_Genomes', render: renderCell},
        {title: 'SAS Allele frequency', prop: 'SAS_Allele_frequency_1000_Genomes', render: renderCell},
        {title: 'EAS Allele frequency', prop: 'EAS_Allele_frequency_1000_Genomes', render: renderCell},
        {title: 'AMR Allele frequency', prop: 'AMR_Allele_frequency_1000_Genomes', render: renderCell},
        {title: 'EUR Allele frequency', prop: 'EUR_Allele_frequency_1000_Genomes', render: renderCell},
        {title: 'AFR Allele frequency', prop: 'AFR_Allele_frequency_1000_Genomes', render: renderCell}
     ]
    },
    {subColTitle: 'ClinVar',
     subColList: [
        {title: 'Allele origin', prop: 'Allele_origin_ClinVar', render: renderCell},
        {title: 'Variant clinical significance', prop: 'Variant_clinical_significance_ClinVar', render: renderCell}
     ]
    },
    {subColTitle: 'LOVD',
     subColList: [
        {title: 'HGVS genomic', prop: 'HGVS_genomic_LOVD', render: renderCell},
        {title: 'Origin of variant', prop: 'Origin_of_variant_LOVD', render: renderCell},
        {title: 'HGVS protein', prop: 'HGVS_protein_LOVD', render: renderCell},
        {title: 'Variant frequency', prop: 'Variant_frequency_LOVD', render: renderCell},
        {title: 'HGVS cDNA', prop: 'HGVS_cDNA_LOVD', render: renderCell},
        {title: 'Variant affecting protein', prop: 'Variant_affecting_protein_LOVD', render: renderCell},
        {title: 'Variant haplotype', prop: 'Variant_haplotype_LOVD', render: renderCell}
     ]
    },
    {subColTitle: 'ExAC',
     subColList: [
        {title: 'VEP Gene', prop: 'VEP_Gene_ExAC', render: renderCell},
        {title: 'Allele frequency', prop: 'Allele_frequency_ExAC', render: renderCell},
        {title: 'VEP HGVSc', prop: 'VEP_HGVSc_ExAC', render: renderCell},
        {title: 'VEP Consequence', prop: 'VEP_Consequence_ExAC', render: renderCell},
        {title: 'VEP HGVSp', prop: 'VEP_HGVSp_ExAC', render: renderCell}
     ]
    },
    {subColTitle: 'exLOVD',
     subColList: [
        {title: 'Exon number', prop: 'Exon_number_exLOVD', render: renderCell},
        {title: 'IARC class', prop: 'IARC_class_exLOVD', render: renderCell},
        {title: 'BIC', prop: 'BIC_exLOVD', render: renderCell},
        {title: 'HGVS cDNA', prop: 'HGVS_cDNA_exLOVD', render: renderCell},
        {title: 'Literature source', prop: 'Literature_source_exLOVD', render: renderCell},
        {title: 'HGVS protein', prop: 'HGVS_protein_exLOVD', render: renderCell}
     ]
    }
];

var columnSelection = {
    Gene_symbol: {selectVal: true},
    Genomic_Coordinate: {selectVal: true},
    HGVS_cDNA: {selectVal: true},
    HGVS_protein: {selectVal: true},
    Abbrev_AA_change: {selectVal: true},
    BIC_Nomenclature: {selectVal: true},
    Clinical_significance: {selectVal: true},
    Allele_frequency_1000_Genomes: {selectVal: false},
    SAS_Allele_frequency_1000_Genomes: {selectVal: false},
    EAS_Allele_frequency_1000_Genomes: {selectVal: false},
    AMR_Allele_frequency_1000_Genomes: {selectVal: false},
    EUR_Allele_frequency_1000_Genomes: {selectVal: false},
    AFR_Allele_frequency_1000_Genomes: {selectVal: false},
    Allele_origin_ClinVar: {selectVal: false},
    Variant_clinical_significance_ClinVar: {selectVal: false},
    HGVS_genomic_LOVD: {selectVal: false},
    Origin_of_variant_LOVD: {selectVal: false},
    HGVS_protein_LOVD: {selectVal: false},
    Variant_frequency_LOVD: {selectVal: false},
    HGVS_cDNA_LOVD: {selectVal: false},
    Variant_affecting_protein_LOVD: {selectVal: false},
    Variant_haplotype_LOVD: {selectVal: false},
    VEP_Gene_ExAC: {selectVal: false},
    Allele_frequency_ExAC: {selectVal: false},
    VEP_HGVSc_ExAC: {selectVal: false},
    VEP_Consequence_ExAC: {selectVal: false},
    VEP_HGVSp_ExAC: {selectVal: false},
    Exon_number_exLOVD: {selectVal: false},
    IARC_class_exLOVD: {selectVal: false},
    BIC_exLOVD: {selectVal: false},
    HGVS_cDNA_exLOVD: {selectVal: false},
    Literature_source_exLOVD: {selectVal: false},
    HGVS_protein_exLOVD: {selectVal: false}
};

var sources = {
    Variant_in_ENIGMA: true,
    Variant_in_ClinVar: true,
    Variant_in_1000_Genomes: true,
    Variant_in_ExAC: true,
    Variant_in_LOVD: true,
    Variant_in_BIC: true
}
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

var VariantTable = React.createClass({
    mixins: [PureRenderMixin],
    getData: function () {
        return this.refs.table.state.data;
    },
    render: function () {
        var {data, onHeaderClick, onRowClick, ...opts} = this.props;
        return (
            <DataTable
                ref='table'
                className='row-clickable'
                {...opts}
                buildRowOptions={r => ({title: 'click for details', onClick: () => hasSelection() ? null : onRowClick(r)})}
                buildHeader={title => buildHeader(onHeaderClick, title)}
                filterColumns={filterColumns}
                columns={columns}
                subColumns={subColumns}
                columnSelection={columnSelection}
                source={sources}
                initialData={data}
                initialPageLength={20}
                initialSortBy={{prop: 'Abbrev_AA_change', order: 'descending'}}
                pageLengthOptions={[ 20, 50, 100 ]}
            />
        );
    }
});


module.exports = VariantTable;
