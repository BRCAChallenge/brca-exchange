/*eslint-env browser */
/*global require: false, module: false */
'use strict';

var content = {
	home: require('../content/home.md'),
	history: require('../content/history.md'),
	brca1_2: require('../content/brca1_2.md'), //eslint-disable-line camelcase
	variation: require('../content/variationAndCancer.md'),
	help: require('../content/help.md'),
    disclaimer: require('../content/disclaimer.md')
};

module.exports = {
	pages: content
};
