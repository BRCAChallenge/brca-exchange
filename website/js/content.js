'use strict';

import jQuery from 'jquery';
import flattenDeep from 'lodash/flattenDeep';
import _ from 'lodash';

// Content pages
import homeMd from '../content/home.md';
import historyMd from '../content/history.md';
import variationMd from '../content/variationAndCancer.md';
import disclaimerMd from '../content/disclaimer.md';
import thisSiteMd from '../content/thisSite.md';
import dataSubmissionPolicyMd from '../content/dataSubmissionPolicy.md';
import apiMd from '../content/api.md';
import variantsDefaultMd from '../content/variantsDefault.md';
import variantsResearchMd from '../content/variantsResearch.md';
import researchWarningMd from '../content/researchWarning.md';
import signupMessageMd from '../content/signupMessage.md';
import insilcoScoringMd from '../content/help/research/insilico/insilico-scoring.md';
import moreOnNamingMd from '../content/help/default/variant-naming/more-on-naming.md';
import appMd from '../content/about/app.md';
import whyDonateMd from '../content/whyDonate.md';
import whyDonateSubtextMd from '../content/whyDonateSubtext.md';
import fundraisingDetailsMd from '../content/fundraisingDetails.md';

// MuPIT images
import mupit1t15 from '../content/mupit/1t15.png';
import mupit1jm7 from '../content/mupit/1jm7.png';
import mupit4igk from '../content/mupit/4igk.png';
import mupitENSP from '../content/mupit/fENSP00000380152_7.png';

// Help content - Default
import howDoISearchMd from '../content/help/default/how-do-i-search-for-a-variant.md';
import variantClassificationsPathogenicMd from '../content/help/default/variant-classifications-pathogenic.md';
import variantClassificationsBenignMd from '../content/help/default/variant-classifications-benign.md';
import variantClassificationsVUSMd from '../content/help/default/variant-classifications-VUS.md';
import variantClassificationsNotYetReviewedMd from '../content/help/default/variant-classifications-not-yet-reviewed.md';
import whatDoesItMeanNotYetReviewedMd from '../content/help/default/what-does-it-mean-if-a-variant-is-not-yet-reviewed.md';
import demystifyingVariantNamingMd from '../content/help/default/variant-naming/demystifying-variant-naming.md';
import whyCantIFindVariantMd from '../content/help/default/why-cant-i-find-the-variant-i-searched.md';
import variantDetailsFieldsMd from '../content/help/default/variant-details-fields.md';
import howIsDataUpdatedMd from '../content/help/default/how-is-the-data-on-this-site-updated.md';
import whatIsBRCAExchangeMd from '../content/help/default/what-is-the-brca-exchange.md';
import whatMakesBRCADifferentMd from '../content/help/default/what-makes-brca-exchange-different-from-other-public-databases.md';
import whatIsEnigmaMd from '../content/help/default/what-is-enigma-and-how-does-it-determine-variant-classifications.md';
import outsideResourcesMd from '../content/help/default/outside-resources-on-brca-genes.md';

// Help content - Research
import howDoISearchResearchMd from '../content/help/research/how-do-i-search-for-a-variant.md';
import whyCantIFindVariantResearchMd from '../content/help/research/why-cant-i-find-the-variant-i-searched.md';
import usingColumnSelectorsMd from '../content/help/research/using-column-selectors.md';
import variantNomenclatureMd from '../content/help/research/variant-nomenclature.md';
import clinicalSignificanceEnigmaMd from '../content/help/research/clinical-significance-enigma.md';
import clinicalSignificanceClinvarMd from '../content/help/research/clinical-significance-clinvar.md';
import clinicalSignificanceLovdMd from '../content/help/research/clinical-significance-lovd.md';
import transcriptVisualizationMd from '../content/help/research/transcript-visualization.md';
import multifactorialLikelihoodMd from '../content/help/research/multifactorial-likelihood-analysis.md';
import alleleFrequencyRefSetsMd from '../content/help/research/allele-frequency-reference-sets.md';
import alleleFrequencyGnomadMd from '../content/help/research/allele-frequency-gnomad.md';
import acmgVariantEvidenceCodesMd from '../content/help/research/acmg-variant-evidence-codes-provisional-assignment.md';
import acmgPopulationFrequencyMd from '../content/help/research/acmg-variant-evidence-codes-population-frequency.md';
import acmgComputationalPredictionMd from '../content/help/research/acmg-variant-evidence-codes-computational-prediction.md';
import functionalAssayDisclaimerMd from '../content/help/research/functional-assay-disclaimer.md';
import functionalAssaysMd from '../content/help/research/functional-assays.md';
import computationalPredictionsMd from '../content/help/research/computational-predictions.md';
import insilicoPredMd from '../content/help/research/insilico/insilico-pred.md';
import cravatMupitMd from '../content/help/research/cravat-mupit.md';
import literatureSearchMd from '../content/help/research/literature-search.md';
import whatGenomeBuildMd from '../content/help/research/what-genome-build-is-this-site-using.md';
import howIsDataUpdatedResearchMd from '../content/help/research/how-is-the-data-on-this-site-updated.md';
import downloadingVariantDataMd from '../content/help/research/downloading-variant-data.md';

const content = {
    home: homeMd,
    history: historyMd,
    variation: variationMd,
    disclaimer: disclaimerMd,
    thisSite: thisSiteMd,
    dataSubmissionPolicy: dataSubmissionPolicyMd,
    api: apiMd,
    variantsDefault: variantsDefaultMd,
    variantsResearch: variantsResearchMd,
    researchWarning: researchWarningMd,
    signupMessage: signupMessageMd,
    insilicoScoring: insilcoScoringMd,
    moreOnNaming: moreOnNamingMd,
    app: appMd,
    whyDonate: whyDonateMd,
    whyDonateSubtext: whyDonateSubtextMd,
    fundraisingDetails: fundraisingDetailsMd
};

export const mupitStructures = [
    {
        "name": "1t15",
        "image": mupit1t15,
        "url": "http://mupit.icm.jhu.edu/MuPIT_Interactive?structure_id=1t15&addtlinfo=brca",
        "humanReadableName": "BRCA1 BRCT Domain"
    },
    {
        "name": "1jm7",
        "image": mupit1jm7,
        "url": "http://mupit.icm.jhu.edu/MuPIT_Interactive?structure_id=1jm7&addtlinfo=brca",
        "humanReadableName": "BRCA1 Ring Domain"
    },
    {
        "name": "4igk",
        "image": mupit4igk,
        "url": "http://mupit.icm.jhu.edu/MuPIT_Interactive?structure_id=4igk&addtlinfo=brca",
        "humanReadableName": "BRCA1 BRCT Domain"
    },
    {
        "name": "fENSP00000380152_7",
        "image": mupitENSP,
        "url": "http://mupit.icm.jhu.edu/MuPIT_Interactive/?gene=BRCA2&addtlinfo=brca",
        "humanReadableName": "BRCA2 Homology Model"
    }
];

export const helpContentDefault = [
    {
	section: "Getting Started: Variants and How to Find Them",
        tiles: [
           {
                name: "How do I search for a variant?",
                contents: howDoISearchMd
           },
           {
                name: "What does a variant's classification mean?",
                list: [
                    {
                        name: "Pathogenic",
                        contents: variantClassificationsPathogenicMd
                    },
                    {
                        name: "Benign/Little Clinical Significance",
                        contents: variantClassificationsBenignMd
                    },
                    {
                        name: "Variant of Uncertain Significance",
                        contents: variantClassificationsVUSMd
                    },
                    {
                        name: "Not Yet Reviewed",
                        contents: variantClassificationsNotYetReviewedMd
                    },
                  ]
           },
           {
		   name: "What does it mean if a variant is \"Not Yet Reviewed\"?",
		   contents: whatDoesItMeanNotYetReviewedMd
           },
           {
                   name: "How do variants get their names?",
                   contents: demystifyingVariantNamingMd
           },
           {
                   name: "Why can't I find the variant I searched?",
                   contents: whyCantIFindVariantMd
           },
        ]
    },
    {
	section: "Variant Details: Understanding the Data",
        tiles: [
           {
                name: "What do the fields in the Variant Details Page mean?",
                contents: variantDetailsFieldsMd
           },
           {
		   name: "How is the data on this site updated?",
                contents: howIsDataUpdatedMd
	   },
	]
    },
    {
	section: "Who We Are: Understanding the BRCA Exchange",
        tiles: [
            {
                name: "What is the BRCA Exchange?",
                contents: whatIsBRCAExchangeMd
            },
            {
                name: "What makes BRCA Exchange different from other public databases?",
                contents: whatMakesBRCADifferentMd
            },
            {
                name: "What is ENIGMA and how does it determine variant classifications?",
                contents: whatIsEnigmaMd
            },
        ]
    },
    {
            section: "More Resources",
        tiles: [
            {
                name: "Where can I find additional information on BRCA Genes, Genetic Testing, and Cancer Risk?",
                contents: outsideResourcesMd
            },
        ]
    },
];

export const helpContentResearch = [
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
	section: "Getting Started: Variants and How to Find Them",
        tiles: [
            {
                name: "How do I search for a variant?",
                contents: howDoISearchResearchMd
            },
            {
                name: "Why can't I find the variant I searched?",
                contents: whyCantIFindVariantResearchMd
            },
            {
                name: "How do I use column selectors?",
                contents: usingColumnSelectorsMd
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
                contents: variantNomenclatureMd
            },
            {
                name: "What data sources and fields are provided in Clinical Significance Tiles?",
                list: [
                    {
                        name: "ENIGMA",
                        id: "clinical-significance-enigma",
                        contents: clinicalSignificanceEnigmaMd
                    },
                    {
                        name: "ClinVar",
                        id: "clinical-significance-clinvar",
                        contents: clinicalSignificanceClinvarMd
                    },
                    {
                        name: "Leiden Open Variation Database (LOVD)",
                        id: "clinical-significance-lovd",
                        contents: clinicalSignificanceLovdMd
                    },
                ]
            },
            {
                name: "What information is displayed in the Transcript Visualization?",
                id: "transcript-visualization",
                contents: transcriptVisualizationMd
            },
            {
                name: "What data is used in Multifactorial Likelihood Analysis?",
                id: "multifactorial-likelihood-analysis",
                contents: multifactorialLikelihoodMd,
                // reference: "https://www.ncbi.nlm.nih.gov/pubmed/21990134"
            },
            {
                name: "What information is provided in the Allele Frequency Reference Sets tile?",
                id: "allele-frequency-reference-sets",
                contents: alleleFrequencyRefSetsMd,
                list: [
                    {
                        name: "gnomAD: Genome Aggregation Database (non-cancer cohort)",
                        contents: alleleFrequencyGnomadMd
                    },
                ]
            },
            {
                name: "What information is provided in the ACMG Variant Evidence Codes, Provisional Assignment tile?",
                id: "acmg-variant-evidence-codes-provisional-assignment",
                contents: acmgVariantEvidenceCodesMd,
                list: [
                    {
                        name: "Population Frequency",
                        id: "acmg-variant-evidence-codes-population-frequency",
                        contents: acmgPopulationFrequencyMd
                    },
                    {
                        name: "Computational Prediction",
                        id: "acmg-variant-evidence-codes-computational-prediction",
                        contents: acmgComputationalPredictionMd
                    },
                ]
            },
            {
                name: "What sources are available in the Functional Assay Results tile?",
                id: "functional-assay-results",
                contents: functionalAssayDisclaimerMd,
                list: [
                    {
                        name: "Functional Assay Scores",
                        contents: functionalAssaysMd
                    },
                ]
            },
            {
                name: "What sources are available in the Computational Predictions tile?",
                id: "computational-predictions",
                contents: computationalPredictionsMd,
                list: [
                    {
                        name: "Computational Predictions",
                        contents: computationalPredictionsMd
                    },
                ]
            },
            {
                name: "What are In Silico Prior Probabilities of Pathogenicity?",
                id: "in-silico-prior-probabilities-of-pathogenicity",
                contents: insilicoPredMd
            },
            {
                name: "How do I interpret the CRAVAT/MuPIT Interactive Protein Structure Viewer?",
                id: "cravat-mupit-3d-protein-view",
                contents: cravatMupitMd
            },
            {
                name: "How does the BRCA Exchange provide Literature Search Results for variants?",
                contents: literatureSearchMd,
                isBeta: true
            },
            {
                name: "What genome build is this site using?",
                contents: whatGenomeBuildMd
            },
            {
                name: "How is the data on this site updated?",
                contents: howIsDataUpdatedResearchMd,
            },
            {
                name: "How do I Download Variant Data?",
                contents: downloadingVariantDataMd
            }
        ]
    },
    {
        section: "Who We Are: Understanding the BRCA Exchange",
        tiles: [
            {
                name: "What is the BRCA Exchange?",
                contents: whatIsBRCAExchangeMd
            },
            {
                name: "What makes BRCA Exchange different from other public databases?",
                contents: whatMakesBRCADifferentMd
            },
            {
                name: "What is ENIGMA and how does it determine variant classifications?",
                contents: whatIsEnigmaMd
            },
        ]
    },
];

export const resources = [
    {
        section: "General Resources",
        tiles: [
            {
                name: "Where can I find additional information on BRCA Genes, Genetic Testing, and Cancer Risk?",
                contents: outsideResourcesMd
            },
        ]
    }
]

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

export function parseContentForTips(helpContent) {
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
export function parseTooltips(isResearchMode) {
    // extract help text depending on the research mode
    // this recursively searches the nested help structure for 'contents' markdown nodes to scrape
    const nodes = flattenDeep(findContentNodes(isResearchMode ? helpContentResearch : helpContentDefault));

    // merge all the nodes' respective dictionaries into one master dictionary
    return _.reduce(nodes.map(node => parseContentForTips(node)), _.extend);
}

// Named export for pages
export const pages = content;

// Default export for backward compatibility
export default {
    pages: content,
    mupitStructures,
    parseTooltips,
    helpContentDefault,
    helpContentResearch,
    resources,
};
