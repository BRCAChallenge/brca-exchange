// HGVS support.
//
// For now, just a quick hack to recognize some reference ids.
// The filters method returns filterValues settings as expected
// by react-data-components DataMixin. The hgvsGene filter will
// filter by gene. The visibleSearch setting will filter by
// free-text match on all fields. When we recognize a hgvs
// coordinate with reference id, we set the hgvsGene filter to
// the appropriate value, and cut the reference id from the
// free-text search.
//
// The effect of this is that <reference-id>:<any> will match
// any field on the variant, while restricting the results to
// the given gene. This may not be ideal, but is a compromise given
// the schedule.

/*global module: false, require: false */
'use strict';

var _ = require('underscore');

var hgvsPatterns = [
	{pat: /^(NM_007294\.3(\(BRCA1\))?|BRCA1):/i, gene: 'BRCA1'},
	{pat: /^(NM_000059\.3(\(BRCA2\))?|BRCA2):/i, gene: 'BRCA2'}
];

function filters(text) {
	var hgvs = _.find(hgvsPatterns, p => text.match(p.pat));
	return hgvs ? {
		hgvsGene: hgvs.gene,
		visibleSearch: text.replace(hgvs.pat, '')
	} : {
		hgvsGene: null,
		visibleSearch: text
	};
}

module.exports = {
	filters: filters
};
