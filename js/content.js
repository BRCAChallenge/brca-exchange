/*eslint-env browser */
/*global require: false, module: false */
'use strict';

var content = {
	home: require('../content/home.md'),
	history: require('../content/history.md'),
	brca1_2: require('../content/brca1_2.md'), //eslint-disable-line camelcase
	variation: require('../content/variationAndCancer.md'),
	help: require('../content/help.md')
};

var lollipopDomains = {
	brca1: JSON.stringify(require('../content/brca1LollipopDomain.json')),
	brca2: JSON.stringify(require('../content/brca2LollipopDomain.json'))
};


module.exports = {
	pages: content,
	lollipop: lollipopDomains
};
