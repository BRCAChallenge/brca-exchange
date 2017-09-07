/*global require: false, module: false */
'use strict';

var Rx = require('rx');
require('rx-dom');
var _ = require('underscore');
var config = require('./config');
var qs = require('qs');

var transpose = a => _.zip.apply(_, a);

// URIs have a 2083 character size limit and some search terms exceed that.
// Limit the length and cut at a semicolon if possible to ensure the search works
function trimSearchTerm(search) {
    var maxLength = 200;
    if (search.length > maxLength) {
        search = search.slice(0, maxLength);
        var lastColonPosition = search.lastIndexOf(":");
        if (lastColonPosition !== -1) {
            search = search.slice(0, lastColonPosition);
        }
    }
    return search;
}

// XXX these defaults might produce odd user experience, since they
// are not reflected in the UI.
function url(opts) {
    var {
        format = 'json',
        filterValues = {},
        sortBy: {prop = 'Gene_Symbol', order = 'ascending'} = {},
        pageLength = 100,
        page = 0,
        search = '',
        column,
        include,
        exclude,
        release,
        changeTypes,
        showDeleted
        } = opts,

        [filter, filterValue] = transpose(_.pairs(_.pick(filterValues, v => v)));
    search = trimSearchTerm(search);

    return `${config.backend_url}/data/?${qs.stringify(_.pick({
        format,
        filter,
        filterValue,
        'order_by': prop,
        direction: order,
        'page_size': pageLength,
        'page_num': page,
        'search_term': search,
        'column': column,
        'include': include,
        'exclude': exclude,
        'release': release,
        'change_types': changeTypes,
        'show_deleted': showDeleted
    }, v => v != null), {arrayFormat: 'repeat'})}`;
}

function data(opts) {
    return Rx.DOM.get(url(opts)).map(xhr => JSON.parse(xhr.responseText));
}

function variant(variant) {
    return Rx.DOM.get(`${config.backend_url}/data/variant/?variant_id=${variant}`).map(xhr => JSON.parse(xhr.responseText));
}

function variantCounts() {
    return Rx.DOM.get(`${config.backend_url}/data/variantcounts`).map(xhr => JSON.parse(xhr.responseText));
}

function releases() {
    return Rx.DOM.get(`${config.backend_url}/data/releases`).map(xhr => JSON.parse(xhr.responseText));
}

function release(id) {
    return Rx.DOM.get(`${config.backend_url}/data/releases?release_id=${id}`).map(xhr => JSON.parse(xhr.responseText));
}

function users(opts) {
    var {page, pageLength, search} = opts;
    var usersUrl = `${config.backend_url}/accounts/users/?${qs.stringify(_.pick({
        'page_num': page,
        'page_size': pageLength,
        'search': search
    }, v => v != null))}`;
    return Rx.DOM.get(usersUrl).map(xhr => JSON.parse(xhr.responseText));
}

function userLocations(search, roles) {
    var url = `${config.backend_url}/accounts/user_locations/?${qs.stringify({search: search, 'roles[]': roles}, {indices: false})}`;
    return Rx.DOM.get(url).map(xhr => JSON.parse(xhr.responseText));
}

function lollipopData(opts) {
    opts.pageLength = 0;
    opts.format = 'json';
    opts.column = ['id', 'Genomic_Coordinate_hg38', 'Pathogenicity_expert', 'Allele_Frequency'];
    return Rx.DOM.get(url(opts)).map(xhr => JSON.parse(xhr.responseText));
}

module.exports = {
    data,
    variant,
    variantCounts,
    releases,
    release,
    users,
    userLocations,
    lollipopData,
    url,
    trimSearchTerm
};
