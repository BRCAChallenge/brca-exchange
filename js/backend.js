/*global require: false, module: false */
'use strict';

var Rx = require('rx');
require('rx-dom');
var _ = require('underscore');

var databaseUrl = "http://localhost:8000";

function filters(filterValues) {
	return _.map(_.pick(filterValues, v => v != null),
				(v, k) => `filter=${k}&filterValue=${v}`).join('&');
}

function searchColumns(columns) {
	return _.map(columns, c => `search_column=${c}`).join('&');
}

// XXX these defaults might produce odd user experience, since they
// are not reflected in the UI.
function data(opts) {
	var {
		filterValues = {},
		source = '',
	    sortBy: {prop = 'Gene_symbol', order = 'ascending'},
		pageLength = 100,
		page = 1,
		search = '',
		searchColumn = ['Variant_Source', 'Gene_symbol']} = opts,
		// XXX use a proper url escape API!!
		query = `${databaseUrl}/data?${filters(filterValues)}&source=${source}&order_by=${prop}&direction=${order}&page_size=${pageLength}&page_num=${page}&search_term=${search}&${searchColumns(searchColumn)}`;

	return Rx.DOM.get(query).map(xhr => JSON.parse(xhr.responseText));
}

module.exports = {
	data
};
