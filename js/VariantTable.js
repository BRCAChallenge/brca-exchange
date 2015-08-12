/*global module: false, require: false, console: false */
'use strict';

var React = require('react');
var _ = require('underscore');

var DataGrid = require('react-datagrid');
require('react-datagrid/index.css');

var merge = (...objs) => _.extend({}, ...objs);

function sanitize(data, fields) {
	var select = _.map(data, r => _.pick(r, fields));
	return _.map(select, (r, i) => merge({id: i}, r));
}

var VariantTable = React.createClass({
	getInitialState: function () {
		var {data} = this.props,
			cleaned = sanitize(data.records, ["CHROM", "POS", "ID", "REF", "ALT"]);
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
				{name: "ID", title: "ID"},
				{name: "CHROM", title: "Chrom"},
				{name: "POS", title: "Position"},
				{name: "REF", title: "Reference"},
				{name: "ALT", title: "Alternate"}
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
