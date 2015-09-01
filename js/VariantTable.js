/*global module: false, require: false, console: false */
'use strict';

var React = require('react');
var _ = require('underscore');

var DataGrid = require('react-datagrid');
require('react-datagrid/index.css');

var merge = (...objs) => _.extend({}, ...objs);

function mergeInfo(row) {
	var info = _.object(_.map(_.pairs(row.INFO), ([k, v]) => ['INFO$' + k, v]));
	return merge(info, _.omit(row, ['INFO']));
}

// XXX hard-coded GENE, PROB, REFS, PATH for now
function sanitize(data) {
	return _.map(data, (r, i) => merge({id: i, GENE: 'BRCA1', PROB: 0.23, REFS: 'sciencemag.org/content', PATH: 'pathogenic'}, mergeInfo(r)));
}

var VariantTable = React.createClass({
	getInitialState: function () {
		var {data} = this.props,
			cleaned = sanitize(data.records);
		// get initial sort of data, same as passed in.
		return {data: cleaned};
	},
	onSortChange: function (info) {
		var [{name}] = info, {data} = this.state;
		this.setState({data: _.sortBy(data, r => r[name])});
	},
	render: function () {
		var {data} = this.state,
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
				style={{height: "30em"}}
				columns={columns}
				dataSource={data}
				onSortChange={this.onSortChange}
				idProperty="id"/>
		);
	}
});

module.exports = VariantTable;
