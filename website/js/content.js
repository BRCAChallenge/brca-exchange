/*eslint-env browser */
/*global require: false, module: false */
'use strict';

var content = {
    home: require('../content/home.md'),
    history: require('../content/history.md'),
    variation: require('../content/variationAndCancer.md'),
    help: require('../content/help.md'),
    helpResearch: require('../content/help_research.md'),
    disclaimer: require('../content/disclaimer.md'),
    thisSite: require('../content/thisSite.md'),
    dataUse: require('../content/dataUse.md'),
    api: require('../content/api.md'),
    variantsDefault: require('../content/variantsDefault.md'),
    variantsResearch: require('../content/variantsResearch.md'),
    researchWarning: require('../content/researchWarning.md'),
    signupMessage: require('../content/signupMessage.md')
};

module.exports = {
    pages: content
};
