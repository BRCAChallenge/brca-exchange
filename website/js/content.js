/*eslint-env browser */
/*global require: false, module: false */
'use strict';

var jQuery = require('jquery');

import flattenDeep from 'lodash/flattenDeep';
import _ from 'lodash';

const content = {
    home: require('../content/home.md'),
    history: require('../content/history.md'),
    variation: require('../content/variationAndCancer.md'),
    disclaimer: require('../content/disclaimer.md'),
    thisSite: require('../content/thisSite.md'),
    dataSubmissionPolicy: require('../content/dataSubmissionPolicy.md'),
    api: require('../content/api.md'),
    variantsDefault: require('../content/variantsDefault.md'),
    variantsResearch: require('../content/variantsResearch.md'),
    researchWarning: require('../content/researchWarning.md'),
    signupMessage: require('../content/signupMessage.md'),
    insilicoScoring: require('../content/help/research/insilico/insilico-scoring.md'),
    moreOnNaming: require('../content/help/default/variant-naming/more-on-naming.md'),
    app: require('../content/about/app.md'),
    whyDonate: require('../content/whyDonate.md'),
    whyDonateSubtext: require('../content/whyDonateSubtext.md'),
    fundraisingDetails: require('../content/fundraisingDetails.md')
};

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

const helpContentDefault = [
    {
	section: "Getting Started: Variants and How to Find Them",
        tiles: [
           {
                name: "How do I search for a variant?",
                contents: require("../content/help/default/how-do-i-search-for-a-variant.md")
           },
           {
                name: "What does a variant's classification mean?",
                list: [
                    {
                        name: "Pathogenic",
                        contents: require("../content/help/default/variant-classifications-pathogenic.md")
                    },
                    {
                        name: "Benign/Little Clinical Significance",
                        contents: require("../content/help/default/variant-classifications-benign.md")
                    },
                    {
                        name: "Variant of Uncertain Significance",
                        contents: require("../content/help/default/variant-classifications-VUS.md")
                    },
                    {
                        name: "Not Yet Reviewed",
                        contents: require("../content/help/default/variant-classifications-not-yet-reviewed.md")
                    },
                  ]
           },
           {
		   name: "What does it mean if a variant is \"Not Yet Reviewed\"?",
		   contents: require("../content/help/default/what-does-it-mean-if-a-variant-is-not-yet-reviewed.md")
           },
           {
                   name: "How do variants get their names?",
                   contents: require("../content/help/default/variant-naming/demystifying-variant-naming.md")
           },
           {
                   name: "Why can't I find the variant I searched?",
                   contents: require("../content/help/default/why-cant-i-find-the-variant-i-searched.md")
           },
        ]
    },
    {
	section: "Variant Details: Understanding the Data",
        tiles: [
           {
                name: "What do the fields in the Variant Details Page mean?",
                contents: require("../content/help/default/variant-details-fields.md")
           },
           {
		   name: "How is the data on this site updated?",
                contents: require("../content/help/default/how-is-the-data-on-this-site-updated.md")
	   },
	]
    },
    {
	section: "Who We Are: Understanding the BRCA Exchange",
        tiles: [
            {
                name: "What is the BRCA Exchange?",
                contents: require("../content/help/default/what-is-the-brca-exchange.md")
            },
            {
                name: "What makes BRCA Exchange different from other public databases?",
                contents: require("../content/help/default/what-makes-brca-exchange-different-from-other-public-databases.md")
            },
            {
                name: "What is ENIGMA and how does it determine variant classifications?",
                contents: require("../content/help/default/what-is-enigma-and-how-does-it-determine-variant-classifications.md")
            },
        ]
    },
    {
            section: "More Resources",
        tiles: [
            {
                name: "Where can I find additional information on BRCA Genes, Genetic Testing, and Cancer Risk?",
                contents: require("../content/help/default/outside-resources-on-brca-genes.md")
            },
        ]
    },
];

const helpContentResearch = [
/* Examples:

    {
        section: "Section name goes here",
        tiles: [
            {
                name: "This is the heading of a tile",
                contents: require("../content/help/research/name-of-markdown-file.md")
            },
            {
                name: "This tile has a reference link in the heading",
                // you can specify a custom ID for the header here.
                // to get the ID you need, you can just click a hover link on the site, and then look at the fragment (after the #) in the URL
                id: "custom-id-to-match-hover-links",
                contents: require("../content/help/research/another-tile-content.md"),
                reference: "http://example.com/this-is-a-reference-link-for-the-tile"
            },
            {
                name: "This tile contains an expandable list",
                content: require("../content/help/research/this-content-appears-before-the-list.md"),
                list: [
                    {
                        name: "List item 1",
                        contents: require("../content/help/research/a-list-item-contents.md")
                    },
                    {
                        name: "List item two",
                        id: "custom-id-to-match-hover-links",
                        contents: require("../content/help/research/and-another-list-item.md")
                    },
                ]
            },
        ]
    },

*/
    {
	section: "Getting Started\: Variants and How to Find Them",
        tiles: [
            {
                name: "How do I search for a variant?",
                contents: require("../content/help/research/how-do-i-search-for-a-variant.md")
            },
            {
                name: "Why can't I find the variant I searched?",
                contents: require("../content/help/research/why-cant-i-find-the-variant-i-searched.md")
            },
            {
                name: "How do I use column selectors?",
                contents: require("../content/help/research/using-column-selectors.md")
            },
            /* No content yet
            {
                name: "How do I use Filters?",
                contents: require("../content/help/research/using-filters.md")
            },
            */
        ]
    },
    {
	section: "Variant Details: Understanding the Data",
        tiles: [
            {
                name: "What are the fields in the Variant Nomenclature Tile?",
                id: "variant-nomenclature",
                contents: require("../content/help/research/variant-nomenclature.md")
            },
            {
                name: "What data sources and fields are provided in Clinical Significance Tiles?",
                list: [
                    {
                        name: "ENIGMA",
                        id: "clinical-significance-enigma",
                        contents: require("../content/help/research/clinical-significance-enigma.md")
                    },
                    {
                        name: "ClinVar",
                        id: "clinical-significance-clinvar",
                        contents: require("../content/help/research/clinical-significance-clinvar.md")
                    },
                    {
                        name: "Leiden Open Variation Database (LOVD)",
                        id: "clinical-significance-lovd",
                        contents: require("../content/help/research/clinical-significance-lovd.md")
                    },
                ]
            },
            {
                name: "What information is displayed in the Transcript Visualization?",
                id: "transcript-visualization",
                contents: require("../content/help/research/transcript-visualization.md")
            },
            {
                name: "What data is used in Multifactorial Likelihood Analysis?",
                id: "multifactorial-likelihood-analysis",
                contents: require("../content/help/research/multifactorial-likelihood-analysis.md"),
                // reference: "https://www.ncbi.nlm.nih.gov/pubmed/21990134"
            },
            {
                name: "What information is provided in the Allele Frequency Reference Sets tile?",
                id: "allele-frequency-reference-sets",
                contents: require("../content/help/research/allele-frequency-reference-sets.md"),
                list: [
                    {
                        name: "gnomAD: Genome Aggregation Database (non-cancer cohort)",
                        contents: require("../content/help/research/allele-frequency-gnomad.md")
                    },
                ]
            },
            {
                name: "What information is provided in the ACMG Variant Evidence Codes, Provisional Assignment tile?",
                id: "acmg-variant-evidence-codes-provisional-assignment",
                contents: require("../content/help/research/acmg-variant-evidence-codes-provisional-assignment.md"),
                list: [
                    {
                        name: "Population Frequency",
                        id: "acmg-variant-evidence-codes-population-frequency",
                        contents: require("../content/help/research/acmg-variant-evidence-codes-population-frequency.md")
                    },
                    {
                        name: "Computational Prediction",
                        id: "acmg-variant-evidence-codes-computational-prediction",
                        contents: require("../content/help/research/acmg-variant-evidence-codes-computational-prediction.md")
                    },
                ]
            },
            {
                name: "What sources are available in the Functional Assay Results tile?",
                id: "functional-assay-results",
                contents: require("../content/help/research/functional-assay-disclaimer.md"),
                list: [
                    {
                        name: "Functional Assay Scores",
                        contents: require("../content/help/research/functional-assays.md")
                    },
                ]
            },
            {
                name: "What sources are available in the Computational Predictions tile?",
                id: "computational-predictions",
                contents: require("../content/help/research/computational-predictions.md"),
                list: [
                    {
                        name: "Computational Predictions",
                        contents: require("../content/help/research/computational-predictions.md")
                    },
                ]
            },
            {
                name: "What are In Silico Prior Probabilities of Pathogenicity?",
                id: "in-silico-prior-probabilities-of-pathogenicity",
                contents: require("../content/help/research/insilico/insilico-pred.md")
            },
            {
                name: "How do I interpret the CRAVAT/MuPIT Interactive Protein Structure Viewer?",
                id: "cravat-mupit-3d-protein-view",
                contents: require("../content/help/research/cravat-mupit.md")
            },
            {
                name: "How does the BRCA Exchange provide Literature Search Results for variants?",
                contents: require("../content/help/research/literature-search.md"),
                isBeta: true
            },
            {
                name: "What genome build is this site using?",
                contents: require("../content/help/research/what-genome-build-is-this-site-using.md")
            },
            {
                name: "How is the data on this site updated?",
                contents: require("../content/help/research/how-is-the-data-on-this-site-updated.md"),
            },
            {
                name: "How do I Download Variant Data?",
                contents: require("../content/help/research/downloading-variant-data.md")
            }
        ]
    },
    {
        section: "Who We Are: Understanding the BRCA Exchange",
        tiles: [
            {
                name: "What is the BRCA Exchange?",
                contents: require("../content/help/default/what-is-the-brca-exchange.md")
            },
            {
                name: "What makes BRCA Exchange different from other public databases?",
                contents: require("../content/help/default/what-makes-brca-exchange-different-from-other-public-databases.md")
            },
            {
                name: "What is ENIGMA and how does it determine variant classifications?",
                contents: require("../content/help/default/what-is-enigma-and-how-does-it-determine-variant-classifications.md")
            },
        ]
    },
];

/**
 * Recursively descends into the object 'head', looking for fields named 'content'. Returns an array of the values of these fields.
 * @param head the object in which to look for content nodes
 * @returns {*} an array of content blobs, unless called on a single object
 */
function findContentNodes(head) {
    if (Array.isArray(head)) {
        return head.map(x => findContentNodes(x));
    }
    else {
        if (Array.isArray(head.tiles)) {
            return head.tiles.map(x => findContentNodes(x));
        }
        else if (Array.isArray(head.list)) {
            return head.list.map(x => findContentNodes(x));
        }
        else if (head.contents !== undefined) {
            return head.contents;
        }
    }
}

// helper for parseContentForTips()
function extractNonHeaders(x) {
    const result = x.parent().clone();
    // remove the field name
    result.children('h4,h5').remove();
    // get only the first element (typically a paragraph)
    result.children().slice(1).remove();
    // unwrap paragraphs to remove weird bootstrap styling
    result.find('p').replaceWith(function() { return jQuery(this).html(); });

    return result.html().trim();
}

function parseContentForTips(helpContent) {
    // we enclose the payload in a div because otherwise jQuery can't find top-level elements in the blob
    const helpElem = jQuery.parseHTML("<div>" + helpContent + "</div>");

    // the glossary's formatted as lis with nested h4s or h5s with ids followed by some text within the same parent
    // in rare cases where the terms aren't in a list, we also look for an explicit <span class="term_entry"> tag
    const extracted = jQuery("li,span.term_entry", helpElem).find("h4,h5").map((idx, x) => {
        const $x = jQuery(x);
        const helpText = extractNonHeaders($x);

        return helpText ? {name: $x.attr("id"), text: helpText} : null;
    }).toArray();

    // creates an object {name1: val1, ...} from our [{name1,val1}, ...] array for faster access
    return extracted.reduce((c, x) => { c[x.name] = x.text; return c; }, {});
}


/**
 * Scrapes the help documentation to extract tooltips for fields displayed on the Variant Details page.
 * @param isResearchMode specifies whether to scrape the research mode or expert-reviewed help docs
 * @returns {*} a mapping from slugified field names to HTML help text
 */
function parseTooltips(isResearchMode) {
    // extract help text depending on the research mode
    // this recursively searches the nested help structure for 'contents' markdown nodes to scrape
    const nodes = flattenDeep(findContentNodes(isResearchMode ? helpContentResearch : helpContentDefault));

    // merge all the nodes' respective dictionaries into one master dictionary
    return _.reduce(nodes.map(node => parseContentForTips(node)), _.extend);
}

module.exports = {
    pages: content,
    mupitStructures,
    parseTooltips,
    helpContentDefault,
    helpContentResearch,
};
