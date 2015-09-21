/*global module: false, require: false */
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
		<a href={"http://www.ncbi.nlm.nih.gov/clinvar/?term=" + val}>{val}</a>
	);
}

var filterColumns = [
	{name: 'Gene', prop: 'Gene symbol', values: ['BRCA1', 'BRCA2']},
//	{name: 'Exon', values: ['Any', 1, 2, 3, 4, 5]}, // XXX needs refgene to get exon count
	{name: 'Pathogenicity', prop: 'Clinical significance', values: ['Pathogenic', 'Benign']}
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

var posCmpFn = strPropCmpFn('Genomic Coordinate');

function sortColumns(columns, {prop, order}, data) {
	var sortFn = _.findWhere(columns, {prop: prop}).sortFn || strPropCmpFn(prop),
		sorted = data.slice(0).sort(sortFn);
	if (order === 'descending') {
		sorted.reverse();
	}
	return sorted;
}

var columns = [
	{title: 'Gene', prop: 'Gene symbol'},
	{title: 'HGVS cDNA', prop: 'HGVS_cDNA', sortFn: posCmpFn},
	{title: 'HGVS protein', prop: 'HGVS_protein', sortFn: posCmpFn},
	{title: 'Genomic Coordinate', prop: 'Genomic Coordinate'},
	{title: 'Pathogenicity', prop: 'Clinical significance'},
	{title: 'Classification method', prop: 'Classification method'},
	{title: 'ClinVar Link', prop: 'ClinVarAccession', render: renderClinVarLink}
];

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
				{...opts}
				buildRowOptions={r => ({title: 'click for details', onClick: () => onRowClick(r)})}
				buildHeader={title => buildHeader(onHeaderClick, title)}
				sort={(sb, d) => sortColumns(columns, sb, d)}
				filter={applyFilters}
				filters={filters(columns)}
				filterColumns={filterColumns}
				columns={columns}
				initialData={data}
				initialPageLength={10}
                initialSortBy={{prop: 'Gene symbol', order: 'descending'}}
                pageLengthOptions={[ 10, 50, 100 ]}
            />
		);
	}
});


module.exports = VariantTable;
