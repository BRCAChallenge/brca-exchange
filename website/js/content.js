/*eslint-env browser */
/*global require: false, module: false */
'use strict';

var jQuery = require('jquery');
var slugify = require('slugify');

const content = {
    home: require('../content/home.md'),
    history: require('../content/history.md'),
    variation: require('../content/variationAndCancer.md'),
    help: require('../content/help.md'),
    helpResearch: require('../content/help_research.md'),
    disclaimer: require('../content/disclaimer.md'),
    thisSite: require('../content/thisSite.md'),
    dataSubmissionPolicy: require('../content/dataSubmissionPolicy.md'),
    api: require('../content/api.md'),
    variantsDefault: require('../content/variantsDefault.md'),
    variantsResearch: require('../content/variantsResearch.md'),
    researchWarning: require('../content/researchWarning.md'),
    signupMessage: require('../content/signupMessage.md')
};

const FAQContent = [
    {
        "question": "What is the BRCA Exchange?",
        "content": require('../content/faqs/what-is-the-brca-exchange.md')
    },
    {
        "question": "What is ENIGMA and how does it determine variant classifications?",
        "content": require('../content/faqs/what-is-enigma-and-how-does-it-determine-variant-classifications.md')
    },
    {
        "question": "What does it mean if a variant is \"Not Yet Reviewed\"?",
        "content": require('../content/faqs/what-does-it-mean-if-a-variant-is-not-yet-reviewed.md')
    },
    {
        "question": "Why can't I find the variant I searched?",
        "content": require('../content/faqs/why-cant-i-find-the-variant-i-searched.md')
    },
    {
        "question": "What genome build is this site using?",
        "content": require('../content/faqs/what-genome-build-is-this-site-using.md')
    },
    {
        "question": "How is the data on this site updated?",
        "content": require('../content/faqs/how-is-the-data-on-this-site-updated.md')
    },
    {
        "question": "What makes brca exchange different from other public databases?",
        "content": require('../content/faqs/what-makes-brca-exchange-different-from-other-public-databases.md')
    }
];

const mupitStructures = [
    {
        "name": "1t15",
        "image": require('../content/mupit/1t15.png'),
        "url": "http://mupit.icm.jhu.edu/MuPIT_Interactive?structure_id=1t15&addtlinfo=brca",
        "humanReadableName": "BRCA1 BRCT Domain"
    },
    {
        "name": "1jm7",
        "image": require('../content/mupit/1jm7.png'),
        "url": "http://mupit.icm.jhu.edu/MuPIT_Interactive?structure_id=1jm7&addtlinfo=brca",
        "humanReadableName": "BRCA1 Ring Domain"
    },
    {
        "name": "4igk",
        "image": require('../content/mupit/4igk.png'),
        "url": "http://mupit.icm.jhu.edu/MuPIT_Interactive?structure_id=4igk&addtlinfo=brca",
        "humanReadableName": "BRCA1 BRCT Domain"
    },
    {
        "name": "fENSP00000380152_7",
        "image": require('../content/mupit/fENSP00000380152_7.png'),
        "url": "http://mupit.icm.jhu.edu/MuPIT_Interactive/?gene=BRCA2&addtlinfo=brca",
        "humanReadableName": "BRCA2 Homology Model"
    }
];

/**
 * Scrapes the help documentation to extract tooltips for fields displayed on the Variant Details page.
 * @param isResearchMode specifies whether to scrape the research mode or expert-reviewed help docs
 * @returns {*} a mapping from slugified field names to HTML help text
 */
function parseTooltips(isResearchMode) {
    // extract help text depending on the research mode
    const helpContent = isResearchMode ? content.helpResearch : content.help;

    const helpElem = document.createElement('html');
    helpElem.innerHTML = helpContent;

    let extracted = {};

    if (isResearchMode) {
        // the glossary's kind of implemented as h4 tags with ids followed by some text within the same parent
        extracted = jQuery("h4", helpElem).map((idx, x) => {
            const $x = jQuery(x);
            const helpText = $x.siblings().html();

            return helpText ? {name: $x.attr("id"), text: helpText} : null;
        }).toArray();
    }
    else {
        // for the non-research help text, "tooltippable" entries need to be manually annotated with that class, e.g.
        // <span class="tooltippable"><em>field name</em>: (any html)</span>, with the resulting key being "field-name"
        extracted = jQuery(".tooltippable", helpElem).map((idx, x) => {
            // extract the initial <em> indicating the term name, then everything after the : as the text
            const $term = jQuery("em:first-child", x);
            const fullText = jQuery(x).html();
            const helpText = fullText.substring(fullText.indexOf(":") + 1);

            return helpText ? {name: slugify($term.text()), text: helpText} : null;
        }).toArray();
    }

    // creates an object {name1: val1, ...} from our [{name1,val1}, ...] array for faster access
    return extracted.reduce((c, x) => { c[x.name] = x.text; return c; }, {});
}

module.exports = {
    pages: content,
    faqs: FAQContent,
    mupitStructures: mupitStructures,
    parseTooltips: parseTooltips
};
