/*global require: false, module: false */
'use strict';

// RxJS 7 imports - CHANGED
import { ajax } from 'rxjs/ajax';
import { map } from 'rxjs/operators';

import _ from 'underscore';
import config from './config';
import qs from 'qs';
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
    // CHANGED: Use ajax.getJSON() which automatically parses JSON
    return ajax.getJSON(url(opts));
}

function variant(variant) {
    // CHANGED: Fixed template literal syntax and use ajax.getJSON()
    return ajax.getJSON(`${config.backend_url}/data/variant/?variant_id=${variant}`);
}

function variantReports(variant) {
    // CHANGED: Fixed template literal syntax and use ajax.getJSON()
    return ajax.getJSON(`${config.backend_url}/data/variant/${variant}/reports`);
}

function variantPapers(variant) {
    // CHANGED: Fixed template literal syntax and use ajax.getJSON()
    return ajax.getJSON(`${config.backend_url}/data/variantpapers/?variant_id=${variant}`);
}

function variantCounts() {
    // CHANGED: Fixed template literal syntax and use ajax.getJSON()
    return ajax.getJSON(`${config.backend_url}/data/variantcounts`);
}

function releases() {
    // CHANGED: Fixed template literal syntax and use ajax.getJSON()
    return ajax.getJSON(`${config.backend_url}/data/releases`);
}

function release(id) {
    // CHANGED: Fixed template literal syntax and use ajax.getJSON()
    return ajax.getJSON(`${config.backend_url}/data/releases?release_id=${id}`);
}

function users(opts) {
    var {page, pageLength, search} = opts;
    var usersUrl = `${config.backend_url}/accounts/users/?${qs.stringify(_.pick({
        'page_num': page,
        'page_size': pageLength,
        'search': search
    }, v => v != null))}`;
    
    // CHANGED: Use ajax.getJSON()
    return ajax.getJSON(usersUrl);
}

function userLocations(search, roles) {
    var url = `${config.backend_url}/accounts/user_locations/?${qs.stringify({search: search, 'roles[]': roles}, {indices: false})}`;
    
    // CHANGED: Use ajax.getJSON()
    return ajax.getJSON(url);
}

export default {
    data,
    variant,
    variantReports,
    variantCounts,
    variantPapers,
    releases,
    release,
    users,
    userLocations,
    url,
    trimSearchTerm
};
