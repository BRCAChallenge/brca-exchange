'use strict';

var React = require('react');
var moment = require('moment');
var _ = require('underscore');

function isEmptyField(value) {
    if (Array.isArray(value)) {
        value = value[0];
    }

    if (value === null || (typeof value === 'undefined')) {
        return true;
    }

    var v = value.trim();
    return v === '' || v === '-' || v === 'None';
}


// attempts to parse the given date string using a variety of formats,
// returning the formatted result as something like '08 September 2016'.
// just returns the input if every pattern fails to match
function normalizeDateFieldDisplay(value) {
    // extend this if there are more formats in the future
    const formats = ["MM/DD/YYYY", "YYYY-MM-DD"];

    for (let i = 0; i < formats.length; i++) {
        const q = moment(value, formats[i]);

        if (q.isValid()) {
            return q.format("DD MMMM YYYY");
        }
    }

    return value;
}


// replaces commas with comma-spaces to wrap long lines better, removes blank entries from comma-delimited lists,
// and normalizes blank/null values to a single hyphen
function normalizedFieldDisplay(value) {
    if (value) {
        // replace any number of underscores with spaces
        // make sure commas, if present, wrap
        value = value
            .split(/_+/).join(" ")
            .split(",")
            .map(x => x.trim())
            .filter(x => x && x !== '-')
            .join(", ");

        // ensure that blank entries are always normalized to hyphens
        if (value.trim() === "") {
            value = "-";
        }
    }
    else {
        // similar to above, normalize blank entries to a hyphen
        value = "-";
    }

    return value;
}


function generateLinkToGenomeBrowser(prop, value) {
    let hgVal = (prop === "Genomic_Coordinate_hg38") ? '38' : '19';
    let genomicCoordinate = value;
    let genomicCoordinateElements = genomicCoordinate.split(':');
    let ref = genomicCoordinateElements[2].split('>')[0];
    let position = parseInt(genomicCoordinateElements[1].split('.')[1]);
    let positionRangeStart = position - 1;
    let positionRangeEnd = position + ref.length + 1;
    let positionParameter = (genomicCoordinate.length > 1500) ? positionRangeStart + '-' + positionRangeEnd : genomicCoordinate;
    let genomeBrowserUrl = 'http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg' + hgVal + '&position=' + positionParameter + '&hubUrl=http://brcaexchange.org/trackhubs/hub.txt';
    return <a target="_blank" href={genomeBrowserUrl}>{value}</a>;
}


function getFormattedFieldByProp(prop, variant) {
    let rowItem;

    if (prop === "Gene_Symbol") {
        rowItem = <i>{variant[prop]}</i>;
    } else if (prop === "URL_ENIGMA") {
        if (variant[prop].length) {
            rowItem = <a target="_blank" href={variant[prop]}>link to multifactorial analysis</a>;
        }
    } else if (prop === "SCV_ClinVar" && variant[prop].toLowerCase().indexOf("scv") !== -1) {
        // Link all clinvar submissions back to clinvar
        let accessions = variant[prop].split(',');
        rowItem = [];
        for (let i = 0; i < accessions.length; i++) {
            if (i < (accessions.length - 1)) {
                rowItem.push(<span><a target="_blank" href={"http://www.ncbi.nlm.nih.gov/clinvar/?term=" + accessions[i].trim()}>{accessions[i]}</a>, </span>);
            } else {
                // exclude trailing comma
                rowItem.push(<a target="_blank" href={"http://www.ncbi.nlm.nih.gov/clinvar/?term=" + accessions[i].trim()}>{accessions[i]}</a>);
            }
        }
    } else if (prop === "DBID_LOVD" && variant[prop].toLowerCase().indexOf("brca") !== -1) { // Link all dbid's back to LOVD
        let ids = variant[prop].split(',');
        rowItem = [];
        for (let i = 0; i < ids.length; i++) {
            if (i < (ids.length - 1)) {
                rowItem.push(<span><a target="_blank" href={"http://lovd.nl/" + ids[i].trim()}>{ids[i]}</a>, </span>);
            } else {
                // exclude trailing comma
                rowItem.push(<a target="_blank" href={"http://lovd.nl/" + ids[i].trim()}>{ids[i]}</a>);
            }
        }
    } else if (prop === "Assertion_method_citation_ENIGMA") {
        rowItem = <a target="_blank" href="https://enigmaconsortium.org/library/general-documents/">Enigma Rules version Mar 26, 2015</a>;
    } else if (prop === "Source_URL") {
        if (variant[prop].startsWith("http://hci-exlovd.hci.utah.edu")) {
            rowItem = <a target="_blank" href={variant[prop].split(',')[0]}>link to multifactorial analysis</a>;
        }
    } else if (prop === "Comment_on_clinical_significance_ENIGMA" || prop === "Clinical_significance_citations_ENIGMA") {
        const pubmed = "http://ncbi.nlm.nih.gov/pubmed/";
        rowItem = _.map(variant[prop].split(/PMID:? ?([0-9]+)/), piece =>
            (/^[0-9]+$/.test(piece)) ? <a target="_blank" href={pubmed + piece}>PMID: {piece}</a> : piece);
    } else if (prop === "HGVS_cDNA") {
        rowItem = variant[prop].split(":")[1];
    } else if (prop === "HGVS_Protein") {
        rowItem = variant[prop].split(":")[1];
    } else if (prop === "Date_last_evaluated_ENIGMA" && !isEmptyField(variant[prop])) {
        // try a variety of formats until one works, or just display the value if not?
        rowItem = normalizeDateFieldDisplay(variant[prop]);
    } else if (/Allele_frequency_.*_ExAC/.test(prop)) {
        let count = variant[prop.replace("frequency", "count")],
            number = variant[prop.replace("frequency", "number")];
        rowItem = [variant[prop], <small style={{float: 'right'}}>({count} of {number})</small>];
    } else if (prop === "Genomic_Coordinate_hg38" || prop === "Genomic_Coordinate_hg37") {
        rowItem = generateLinkToGenomeBrowser(prop, variant[prop]);
    } else {
        rowItem = normalizedFieldDisplay(variant[prop]);
    }

    return rowItem;
}


function abbreviatedSubmitter(originalSubmitter) {
    return originalSubmitter
        .replace('Evidence-based Network for the Interpretation of Germline Mutant Alleles (ENIGMA)', 'ENIGMA')
        .replace('Breast Cancer Information Core (BIC)', 'BIC');
}


function sentenceCase(str) {
    return str.replace(/\b\S/g, (t) => t.toUpperCase() );
}


module.exports = {
    isEmptyField,
    normalizeDateFieldDisplay,
    normalizedFieldDisplay,
    generateLinkToGenomeBrowser,
    getFormattedFieldByProp,
    abbreviatedSubmitter,
    sentenceCase
};
