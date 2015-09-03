/*global module: false, require: false, console: false */
'use strict';

var React = require('react');
var _ = require('underscore');

var DataGrid = require('react-datagrid');
require('react-datagrid/index.css');

var VariantTable = React.createClass({
	onSortChange: function (info) {
		var [{name}] = info, {data} = this.state;
		this.setState({data: _.sortBy(data, r => r[name])});
	},
	render: function () {
		var {data, ...opts} = this.props,
			columns = [
				{name: 'GENE', title: 'Gene', width: 80},
				{name: 'INFO$HGVS_G', title: 'HGVS g', width: 250},
				{name: 'INFO$HGVS_C', title: 'HGVS c', width: 200},
				{name: 'INFO$HGVS_P', title: 'HGVS p', width: 200},
				{name: 'INFO$BIC_N', title: 'BIC n', width: 70},
				{name: 'INFO$BIC_P', title: 'BIC p', width: 70},
				{name: 'INFO$DBSource', title: 'Source', width: 80},
				{name: 'INFO$MUTTYPE', title: 'Type', width: 100},
				{name: 'PROB', title: 'Posterior prob', width: 120},
				{name: 'INFO$FREQ', title: 'Allele freq', width: 100},
				{name: 'REFS', title: 'References', width: 150},
				{name: 'INFO$IARC', title: 'IARC Classification', width: 120},
				{name: 'PATH', title: 'Pathogenicity', width: 120}
			];
		return (
			<DataGrid
				{...opts}
				style={{height: "30em"}}
				columns={columns}
				dataSource={data}
				onSortChange={this.onSortChange}
				idProperty="id"/>
		);
	}
});

module.exports = VariantTable;
