/*global module: false, require: false, window: false */
'use strict';

var React = require('react');
var PureRenderMixin = require('./PureRenderMixin');
var DataTable = require('./DataTable');
require('react-data-components-bd2k/css/table-twbs.css');
var _ = require('underscore');
var {utils} = require('react-data-components-bd2k');

function buildHeader(onClick, title) {
	return (
		<span>
			{title}
			<span onClick={ev => {ev.stopPropagation(); onClick(title); }}
				className='help glyphicon glyphicon-question-sign superscript'/>
		</span>
	);
}

function renderClinVarLink(val) {
	return (
		<a title="View on ClinVar"
			onClick={ev => ev.stopPropagation()}
			href={"http://www.ncbi.nlm.nih.gov/clinvar/?term=" + val}>{val}</a>
	);
}

function renderCell(val) {
	return <span>{val}</span>;
}

var filterColumns = [
	{name: 'Gene', prop: 'Gene_symbol', values: ['BRCA1', 'BRCA2']},
//	{name: 'Exon', values: ['Any', 1, 2, 3, 4, 5]}, // XXX needs refgene to get exon count
	{name: 'Pathogenicity', prop: 'Clinical_significance', values: ['Pathogenic', 'Benign']}
];



// This callback is used to apply all active filters. We override the
// one in react-data-components.utils, which performs a union of all
// matches, with this one which does an intersection.
var applyFilters = (filters, filterValues, data) => {
	return _.filter(data, row => _.every(filterValues, utils.filterPass(filters, row)));
};

// react-data-components filters are an object with keys for each filter, and
// values being a object with props 'filter', and optional 'prop' if the filter
// applies to just one property. Here we have filters 'visibleSearch', which does
// a 'string contains' filter on the visible columns, plus filters for each column
// that has a filter UI control.
// {visibleSearch: {filter: ...}, 'Gene symbol': {prop: 'Gene symbol', filter: ...}}
function filters(columns) {
	var visible = _.object(_.map(columns, c => [c.prop, true]));
	var colFilters = _.object(_.map(filterColumns, c =>
		[c.prop, {
			filter: (fv, val) => _.isNull(fv) || fv === val,
			prop: c.prop
		}]
	));
	return _.extend({
		visibleSearch: {
			filter: (filterValue, value, key) => visible[key] &&
				value.toLowerCase().indexOf(filterValue.toLowerCase()) !== -1
		}
	}, colFilters);
}

var strPropCmpFn = prop => (a, b) => {
	var ap = a[prop],
		bp = b[prop];
	if (ap == null && bp == null || ap === bp) {
		return 0;
	}
	if (bp == null || bp < ap) {
		return 1;
	}
	return -1;
};

var posCmpFn = strPropCmpFn('Genomic_Coordinate');

function sortColumns(columns, {prop, order}, data) {
	var sortFn = _.findWhere(columns, {prop: prop}).sortFn || strPropCmpFn(prop),
		sorted = data.slice(0).sort(sortFn);
	if (order === 'descending') {
		sorted.reverse();
	}
	return sorted;
}

var columns = [
	{title: 'Gene', prop: 'Gene_symbol', render: renderCell},
	{title: 'HGVS cDNA', prop: 'HGVS_cDNA', sortFn: posCmpFn, render: renderCell},
	{title: 'HGVS protein', prop: 'HGVS_protein', sortFn: posCmpFn, render: renderCell},
	{title: 'BIC', prop: "BIC_Nomenclature", render: renderCell},
	{title: 'Abbrev AA Change', prop: "Abbrev_AA_change", render: renderCell},
	{title: 'Genomic Coordinate', prop: 'Genomic_Coordinate', render: renderCell},
	{title: 'Pathogenicity', prop: 'Clinical_significance', render: renderCell},
	{title: 'ClinVar Link', prop: 'ClinVarAccession', render: renderClinVarLink}
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
				sort={(sb, d) => sortColumns(columns, sb, d)}
				filter={applyFilters}
				filters={filters(columns)}
				filterColumns={filterColumns}
				columns={columns}
				initialData={data}
				initialPageLength={20}
                initialSortBy={{prop: 'Abbrev_AA_change', order: 'descending'}}
                pageLengthOptions={[ 20, 50, 100 ]}
            />
		);
	}
});


module.exports = VariantTable;
