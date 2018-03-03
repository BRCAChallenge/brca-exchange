/*eslint-env browser */
/*global require: false, module: false */
'use strict';

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

module.exports = {
    pages: content,
    faqs: FAQContent
};
