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

// XXX hard-coded GENE for now
function sanitize(data) {
	return _.map(data, (r, i) => merge({id: i, GENE: 'BRCA1'}, mergeInfo(r)));
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
				{name: 'GENE', title: 'Gene'},
				{name: 'INFO$HGVS_G', title: 'HGVS g'},
				{name: 'INFO$HGVS_C', title: 'HGVS c'},
				{name: 'INFO$HGVS_P', title: 'HGVS p'},
				{name: 'INFO$BIC_N', title: 'BIC n'},
				{name: 'INFO$BIC_P', title: 'BIC p'},
				{name: 'INFO$DBSource', title: 'Source'},
				{name: 'INFO$MUTTYPE', title: 'Type'},
				{name: 'INFO$IARC', title: 'IARC Classification'}
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
