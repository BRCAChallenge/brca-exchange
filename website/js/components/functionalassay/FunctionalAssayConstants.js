/*global module: false, require: false, URL: false, Blob: false */
'use strict';

const FunctionalAssayConstants = [{
        name: "Findlay",
        loc: "Both",
        author: "Findlay et al, 2018",
        publication: "https://pubmed.ncbi.nlm.nih.gov/30209399",
        previousPublications: "",
        resultDescription: "FUNC: Functional, INT: Intermediate, LOF: Loss of function"
    },
    {
        name: "Starita",
        loc: "Protein",
        author: "Starita et al, 2018",
        publication: "https://pubmed.ncbi.nlm.nih.gov/30219179",
        previousPublications: "",
        resultDescription: "Counts depleted: 0, 1, 4"
    },
    {
        name: "Petitalot",
        loc: "Protein",
        author: "Petitalot et al, 2019",
        publication: "https://pubmed.ncbi.nlm.nih.gov/30257991",
        previousPublications: "PMID:20516115,PMID:20378548,PMID:21473589,PMID:23867111,PMID:25748678,PMID:28781887",
        resultDescription: "1P: no impact, 2P: intermediate impact, 3P: severe impact"
    },
    {
        name: "Bouwman1",
        loc: "Protein",
        author: "Bouwman et al, 2013",
        publication: "https://pubmed.ncbi.nlm.nih.gov/23867111",
        previousPublications: "",
        resultDescription: "Neutral, Likely Neutral, Intermediate, Not clear, Likely Deleterious, Deleterious"
    },
    {
        name: "Bouwman2",
        loc: "Protein",
        author: "Bouwman et al, 2020",
        publication: "https://pubmed.ncbi.nlm.nih.gov/32546644",
        previousPublications: "",
        resultDescription: "Neutral, Likely Neutral, Intermediate, Not clear, Likely Deleterious, Deleterious"
    },
    {
        name: "Fernandes",
        loc: "Protein",
        author: "Fernandes et al, 2019",
        publication: "https://pubmed.ncbi.nlm.nih.gov/30765603",
        previousPublications: "PMID:11157798,PMID:12496476,PMID:15172985,PMID:15689452,PMID:17020472,PMID:17308087,PMID:17311832,PMID:18285836,PMID:18992264,PMID:20516115,PMID:21447777,PMID:23613828,PMID:24845084,PMID:28781887",
        resultDescription: "fClass: Functional Class 1, 2, 3, 4, 5"
    },
    {
        name: "Mesman",
        loc: "Both",
        author: "Mesman et al, 2019",
        publication: "https://pubmed.ncbi.nlm.nih.gov/29988080",
        previousPublications: "PMID:2514691",
        resultDescription: "Functional, Intermediate function, Loss of function"
    },
    {
        name: "Richardson",
        loc: "Protein",
        author: "Richardson et al, 2021",
        publication: "https://pubmed.ncbi.nlm.nih.gov/33609447",
        previousPublications: "PMID:29884841,PMID:2310813,PMID:21990134,PMID:29394989",
        resultDescription: "Neutral, Damaging"
    },
    {
        name: "Ikegami",
        loc: "Protein",
        author: "Ikegami et al, 2020",
        publication: "https://pubmed.ncbi.nlm.nih.gov/32444794",
        previousPublications: "",
        resultDescription: "fClass: Functional Class 1, 2, 3, 4, 5"
    },
    {
        name: "Biwas",
        loc: "Protein",
        author: "Biwas et al, 2020",
        publication: "https://pubmed.ncbi.nlm.nih.gov/33293522",
        previousPublications: "",
        resultDescription: "(F): Functional, (I): Intermediate, (NF): Nonfunctional"
    }
];

module.exports = FunctionalAssayConstants;
