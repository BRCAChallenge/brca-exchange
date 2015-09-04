/*global module: false, require: false, console: false */
'use strict';

var React = require('react');
var _ = require('underscore');
var DataTable = require('react-data-components').DataTable;
var DataGrid = require('react-datagrid');
require('react-datagrid/index.css');
require('react-data-components/css/table-twbs.css');

var VariantTable = React.createClass({
	onSortChange: function (info) {
		var [{prop}] = info, {data} = this.state;
		this.setState({data: _.sortBy(data, r => r[prop])});
	},
	render: function () {
		var {data, ...opts} = this.props,
			columns = [
				{title: 'Gene symbol', prop: 'Gene symbol'},
				{title: 'HGVS', prop: 'HGVS'},
				{title: 'Alternate', prop: 'Alternate'},
				{title: 'Clinical significance', prop: 'Clinical significance'},
				{title: 'Allele origin', prop: 'Allele origin'},
				{title: 'ClinVarAcession', prop: 'ClinVarAcession'}
			];
		return (
			<DataTable
				{...opts}
				columns={columns}
				initialData={data}
				initialPageLength={5}
                initialSortBy={{ title: 'Gene', prop: 'Gene', order: 'descending' }}
                pageLengthOptions={[ 5, 20, 50 ]}
                keys={['id']}
            />
		);
	}
});


module.exports = VariantTable;
