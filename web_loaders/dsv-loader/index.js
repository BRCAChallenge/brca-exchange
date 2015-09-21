/*global require: false, module: false */
'use strict';

var loaderUtils = require('loader-utils');
var tsv = require('d3-dsv').tsv;
var _ = require('underscore');

function indexOfOrThrow(a, v) {
	var r = _.indexOf(a, v);
	if (r === -1) {
		throw new Error("No key named " + v);
	}
	return r;
}

function keyValue(indices, row) {
	return JSON.stringify(_.pick(row, indices));
}

module.exports = function(text) {
	this.cacheable();
	var query = loaderUtils.parseQuery(this.query),
		key = query.key,
		rows = tsv.parseRows(text),
		data = {
			header: rows[0],
			rows: rows.slice(1)
		};

	 if (key) {
		 var indices = _.map(key, _.partial(indexOfOrThrow, data.header)),
			 allKeys = _.map(data.rows, _.partial(keyValue, indices));
		 if (allKeys.length !== _.uniq(allKeys).length) {
			 throw new Error("Duplicate primary keys in tsv file. See webpack.config.js.");
		 }
	 }
	 return JSON.stringify(data);
};
