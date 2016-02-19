// HGVS support.
//
// For now, just a quick hack to recognize some reference ids.
// The filters method returns fetch options.  When we recognize a hgvs
// coordinate with reference id, we set the gene filter to the appropriate
// value, and cut the reference id from the free-text search.
//
// The effect of this is that <reference-id>:<any> will match
// any field on the variant, while restricting the results to
// the given gene. This may not be ideal, but is a compromise given
// the schedule.
//
// If the inferred gene conflicts with the user's gene filter, we
// pass the free text unmodified. The effect of this is usually
// zero matches, which is fine.

/*global module: false, require: false */
'use strict';

var _ = require('underscore');

var hgvsPatterns = [
    {pat: /^(NM_007294\.3(\(BRCA1\))?|BRCA1)(:|$)/i, gene: 'BRCA1'},
    {pat: /^(NM_000059\.3(\(BRCA2\))?|BRCA2)(:|$)/i, gene: 'BRCA2'}
];

var gene = 'Gene_symbol';

function filters(search, filterValues) {
    var hgvs = _.find(hgvsPatterns, p => search.match(p.pat));
    return hgvs && (filterValues[gene] == null || filterValues[gene] === hgvs.gene) ?
        {
            search: search.replace(hgvs.pat, ''),
            filterValues: {
                ...filterValues,
                [gene]: hgvs.gene
            }
        } : {
            search,
            filterValues
        };
}

module.exports = {
    filters: filters
};
