/*eslint-env browser */
/*global require: false, module: false */
'use strict';

var content = {
    home: require('../content/home.md'),
    history: require('../content/history.md'),
    variation: require('../content/variationAndCancer.md'),
    help: require('../content/help.md'),
    help_research: require('../content/help_research.md'),
    disclaimer: require('../content/disclaimer.md'),
    thisSite: require('../content/thisSite.md')
};

module.exports = {
    pages: content
};
