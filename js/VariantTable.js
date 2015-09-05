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
				{title: 'Gene', prop: 'Gene symbol'},
				{title: '  HGVS  ', prop: 'HGVS'},
				{title: 'Pathogenicity', prop: 'Clinical significance'},
				{title: 'Allele origin', prop: 'Allele origin'},
				{title: 'CVA', prop: 'ClinVarAccession'}
			];
		return (
			<DataTable
				{...opts}
				columns={columns}
				initialData={data}
				initialPageLength={10}
                initialSortBy={{ title: 'Gene', prop: 'Gene', order: 'descending' }}
                pageLengthOptions={[ 10, 50, 100 ]}
                keys={['id']}
            />
		);
	}
});


module.exports = VariantTable;
