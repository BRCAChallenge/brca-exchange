'use strict';

var React = require('react');
var moment = require('moment');
var _ = require('underscore');

// keys that contain date values that need reformatting for the ui
const dateKeys = [
    "Date_Last_Updated_ClinVar",
    "DateSignificanceLastEvaluated_ClinVar",
    "Date_last_evaluated_ENIGMA",
    "Edited_date_LOVD",
    "Created_date_LOVD"
];

const AminoAcids = {
    'ala': 'a',
    'arg': 'r',
    'asn': 'n',
    'asp': 'd',
    'asx': 'b',
    'cys': 'c',
    'glu': 'e',
    'gln': 'q',
    'glx': 'z',
    'gly': 'g',
    'his': 'h',
    'ile': 'i',
    'leu': 'l',
    'lys': 'k',
    'met': 'm',
    'phe': 'f',
    'pro': 'p',
    'ser': 's',
    'thr': 't',
    'trp': 'w',
    'tyr': 'y',
    'val': 'v'
};


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

function isNumeric(n) {
    return !isNaN(parseFloat(n)) && isFinite(n);
}

function sentenceCase(str) {
    return str.replace(/\b\S/g, (t) => t.toUpperCase() );
}

function capitalize(w) {
    return w.charAt(0).toUpperCase() + w.substr(1);
}

function extractValInsideParens(str) {
    const regExp = /\(([^)]+)\)/;
    return regExp.exec(str)[1];
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
    let genomeBrowserUrl = 'http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg' + hgVal + '&position=' + positionParameter + '&hubUrl=https://brcaexchange.org/trackhubs/hub.txt';
    return <a target="_blank" href={genomeBrowserUrl}>{value}</a>;
}


function reformatDate(date) { //handles single dates or an array of dates
    if (isEmptyField(date)) {
        return date;
    }
    if (!Array.isArray(date)) {
        date = date.split(',');
    }
    return date.map(function(d) {
        return normalizeDateFieldDisplay(d);
    }).join();
}

function formatConditionLink(db, id) {
    let formattedDbId;
    if (db === "MedGen") {
        formattedDbId = "https://www.ncbi.nlm.nih.gov/medgen/" + id;
    } else if (db === "OMIM") {
        formattedDbId = "http://www.omim.org/entry/" + id;
    } else if (db === "Orphanet") {
        formattedDbId = "http://www.orpha.net/consor/cgi-bin/OC_Exp.php?lng=EN&Expert=" + id;
    } else {
        // No url for other sources
        return db;
    }
    return <a target="_blank" href={formattedDbId}>{db}</a>;
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
        let versions = variant["SCV_Version_ClinVar"] ? variant["SCV_Version_ClinVar"].split(',') : null;
        rowItem = [];
        for (let i = 0; i < accessions.length; i++) {
            let displayText = accessions[i];

            if (versions && i < versions.length && versions[i] !== '-') {
                // appending accession version if available
                displayText = accessions[i].concat('.').concat(versions[i]);
            }

            if (i < (accessions.length - 1)) {
                rowItem.push(<span><a target="_blank" href={"http://www.ncbi.nlm.nih.gov/clinvar/?term=" + accessions[i].trim()}>{displayText}</a>,</span>);
            } else {
                // exclude trailing comma
                rowItem.push(<a target="_blank" href={"http://www.ncbi.nlm.nih.gov/clinvar/?term=" + accessions[i].trim()}>{displayText}</a>);
            }
        }
    } else if (prop === "Condition_Value_ClinVar" && !isEmptyField(variant['Condition_DB_ID_ClinVar'])) {
        let dbIds = variant['Condition_DB_ID_ClinVar'].split(',');
        rowItem = [normalizedFieldDisplay(variant['Condition_Value_ClinVar'])];
        rowItem.push(' [');
        for (let i = 0; i < dbIds.length; i++) {
            let dbId = dbIds[i];
            let splitDbId = dbId.split('_');
            let db = splitDbId[0];
            let id = splitDbId[1];
            if (i === (dbIds.length - 1)) {
                let formattedDbId = formatConditionLink(db, id);
                rowItem.push(formattedDbId);
            } else {
                let formattedDbId = formatConditionLink(db, id);
                rowItem.push(formattedDbId);
                rowItem.push(" | ");
            }
        }
        rowItem.push(']');
    }  else if (prop === "Condition_ID_value_ENIGMA" && !isEmptyField(variant['Condition_ID_type_ENIGMA'])) {
        let db = variant['Condition_ID_type_ENIGMA'];
        let id = extractValInsideParens(variant['Condition_ID_value_ENIGMA']);
        let conditionValue = sentenceCase(normalizedFieldDisplay(variant['Condition_ID_value_ENIGMA'].split(';')[0]).toLowerCase());
        rowItem = [conditionValue];
        rowItem.push(' [');
        rowItem.push(formatConditionLink(db, id));
        rowItem.push(']');
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
    } else if (/Allele_frequency_.*_ExAC/.test(prop)) {
        let count = variant[prop.replace("frequency", "count")],
            number = variant[prop.replace("frequency", "number")];
        rowItem = [variant[prop], <small style={{float: 'right'}}>({count} of {number})</small>];
    } else if (prop === "Genomic_Coordinate_hg38" || prop === "Genomic_Coordinate_hg37") {
        rowItem = generateLinkToGenomeBrowser(prop, variant[prop]);
    } else if (prop === "Synonyms") {
        let syns = variant[prop].split(',');
        let synsNoWhitespace = _.map(syns, s => s.replace(' ', '_'));
        rowItem = synsNoWhitespace.join(", ");
    } else {
        rowItem = normalizedFieldDisplay(variant[prop]);
    }

    if (_.contains(dateKeys, prop)) {
        rowItem = reformatDate(variant[prop]);
    }

    return rowItem;
}


function abbreviatedSubmitter(originalSubmitter) {
    return originalSubmitter
        .replace('Evidence-based Network for the Interpretation of Germline Mutant Alleles (ENIGMA)', 'ENIGMA')
        .replace('Breast Cancer Information Core (BIC)', 'BIC');
}


function getAminoAcidCode(hgvsProtein) {
    let trimmedHgvs = hgvsProtein.replace(/[0-9()]/g, '');
    if (trimmedHgvs.length < 3) {
        return false;
    } else {
        let lastThreeChars = trimmedHgvs.substr(trimmedHgvs.length - 3).toLowerCase();
        if (!(lastThreeChars in AminoAcids)) {
            return false;
        } else {
            return AminoAcids[lastThreeChars];
        }
    }
}


module.exports = {
    getAminoAcidCode,
    isEmptyField,
    isNumeric,
    normalizeDateFieldDisplay,
    normalizedFieldDisplay,
    generateLinkToGenomeBrowser,
    getFormattedFieldByProp,
    abbreviatedSubmitter,
    sentenceCase,
    reformatDate,
    dateKeys,
    capitalize
};
