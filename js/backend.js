/*global require: false, module: false */
'use strict';

var Rx = require('rx');
require('rx-dom');
var _ = require('underscore');
var qs = require('qs');

// XXX Need a config mechanism for this. For now, uncomment the localhost
// setting when working locally.
var databaseUrl = "http://brcaexchange.cloudapp.net/backend";
//var databaseUrl = "http://localhost:8000";

var transpose = a => _.zip.apply(_, a);

// XXX these defaults might produce odd user experience, since they
// are not reflected in the UI.
function url(opts) {
    var {
        format = 'json',
        filterValues = {},
        sortBy: {prop = 'Gene_symbol', order = 'ascending'} = {},
        pageLength = 100,
        page = 0,
        search = '',
        source
        } = opts,

        [filter, filterValue] = transpose(_.pairs(_.pick(filterValues, v => v)));

    return `${databaseUrl}/data/?${qs.stringify(_.pick({
        format,
        filter,
        filterValue,
        'order_by': prop,
        direction: order,
        'page_size': pageLength,
        'page_num': page,
        'search_term': search,
        'source': source
    }, v => v != null), {arrayFormat: 'repeat'})}`;
}

function data(opts) {
    return Rx.DOM.get(url(opts)).map(xhr => JSON.parse(xhr.responseText));
}

module.exports = {
    data,
    url
};
