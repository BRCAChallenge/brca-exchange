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
				{prop: 'GENE', title: 'Gene', width: 80},
				{prop: 'INFO$HGVS_G', title: 'HGVS g', width: 250},
				{prop: 'INFO$HGVS_C', title: 'HGVS c', width: 200},
				{prop: 'INFO$HGVS_P', title: 'HGVS p', width: 200},
				{prop: 'INFO$BIC_N', title: 'BIC n', width: 70},
				{prop: 'INFO$BIC_P', title: 'BIC p', width: 70},
				{prop: 'INFO$DBSource', title: 'Source', width: 80},
				{prop: 'INFO$MUTTYPE', title: 'Type', width: 100},
				{prop: 'PROB', title: 'Posterior prob', width: 120},
				{prop: 'INFO$FREQ', title: 'Allele freq', width: 100},
				{prop: 'REFS', title: 'References', width: 150},
				{prop: 'INFO$IARC', title: 'IARC Classification', width: 120},
				{prop: 'PATH', title: 'Pathogenicity', width: 120}
			];
		return (
			<DataTable
				{...opts}
				style={{height: "30em"}}
				columns={columns}
				initialData={data}
				initialPageLength={5}
                initialSortBy={{ prop: 'Gene', order: 'descending' }}
                pageLengthOptions={[ 5, 20, 50 ]}
                keys={['id']}
            />
		);
	}
});


module.exports = VariantTable;
