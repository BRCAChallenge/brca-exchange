'use strict';

import moment from 'moment';
import _ from 'underscore';

// keys that contain date values that need reformatting for the ui
export const dateKeys = [
    "Date_Last_Updated_ClinVar",
    "DateSignificanceLastEvaluated_ClinVar",
    "Date_last_evaluated_ENIGMA",
    "Edited_date_LOVD",
    "Created_date_LOVD"
];

export function isEmptyField(value) {
    if (Array.isArray(value)) {
        value = value[0];
    }

    if (value === null || (typeof value === 'undefined')) {
        return true;
    }

    var v = value.trim();
    return v === '' || v === '-' || v === 'None';
}

export function isNumeric(n) {
    return !isNaN(parseFloat(n)) && isFinite(n);
}

export function sentenceCase(str) {
    return str.replace(/\b\S/g, (t) => t.toUpperCase() );
}

export function capitalize(w) {
    return w.charAt(0).toUpperCase() + w.substr(1);
}

function extractValInsideParens(str) {
    const regExp = /\(([^)]+)\)/;
    return regExp.exec(str)[1];
}


// attempts to parse the given date string using a variety of formats,
// returning the formatted result as something like '08 September 2016'.
// just returns the input if every pattern fails to match
export function normalizeDateFieldDisplay(value) {
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
export function normalizedFieldDisplay(value, prop) {
    if (value) {
        // leave underscores in Refence Sequence field
        if (prop !== "Reference_Sequence" && prop !== "VR_ID") {
            value = value.split(/_+/).join(" ");
        }
        // replace any number of underscores with spaces
        // make sure commas, if present, wrap
        value = value
            .split(",")
            .map(x => x.trim())
            .filter(x => x && x !== '-')
            .join(", ");

        // ensure that blank entries are always normalized to hyphens
        if (value.trim() === "") {
            value = "-";
        }
    } else {
        // similar to above, normalize blank entries to a hyphen
        value = "-";
    }

    return value;
}


export function generateLinkToGenomeBrowser(prop, value, hgvs) {
    let hgVal = (prop === "Genomic_Coordinate_hg38") ? '38' : '19';
    let genomicCoordinate = value;
    let genomicCoordinateElements = genomicCoordinate.split(':');
    let ref = genomicCoordinateElements[2].split('>')[0];
    let position = parseInt(genomicCoordinateElements[1].split('.')[1]);
    let positionRangeStart = position - 1;
    let positionRangeEnd = position + ref.length + 1;
    let positionParameter = (genomicCoordinate.length > 1500) ? positionRangeStart + '-' + positionRangeEnd : genomicCoordinate;
    let genomeBrowserUrl = 'http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg' + hgVal + '&position=' + positionParameter + '&hubUrl=https://brcaexchange.org/trackhubs/hub.txt';
    if (!isEmptyField(hgvs)) {
        value = hgvs;
    }
    return <a target="_blank" href={genomeBrowserUrl} rel="noreferrer">{value}</a>;
}


export function reformatDate(date) { //handles single dates or an array of dates
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
    return <a target="_blank" href={formattedDbId} rel="noreferrer">{db}</a>;
}


export function getFormattedFieldByProp(prop, variant) {
    let rowItem;

    if (prop === "Gene_Symbol") {
        rowItem = <i>{variant[prop]}</i>;
    } else if (prop === "URL_ENIGMA") {
        if (variant[prop].length) {
            rowItem = <a target="_blank" href={variant[prop]} rel="noreferrer">link to multifactorial analysis</a>;
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
                rowItem.push(<span key={`scv-${accessions[i].trim()}-${i}`}><a target="_blank" href={"http://www.ncbi.nlm.nih.gov/clinvar/?term=" + accessions[i].trim()} rel="noreferrer">{displayText}</a>,</span>);
            } else {
                // exclude trailing comma
                rowItem.push(<a key={`scv-${accessions[i].trim()}-${i}`} target="_blank" href={"http://www.ncbi.nlm.nih.gov/clinvar/?term=" + accessions[i].trim()} rel="noreferrer">{displayText}</a>);
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
                rowItem.push(<span key={`lovd-${ids[i].trim()}-${i}`}><a target="_blank" href={"http://lovd.nl/" + ids[i].trim()} rel="noreferrer">{ids[i]}</a>, </span>);
            } else {
                // exclude trailing comma
                rowItem.push(<a key={`lovd-${ids[i].trim()}-${i}`} target="_blank" href={"http://lovd.nl/" + ids[i].trim()} rel="noreferrer">{ids[i]}</a>);
            }
        }
    } else if (prop === "Assertion_method_citation_ENIGMA") {
        rowItem = <a target="_blank" href="https://enigmaconsortium.org/library/general-documents/" rel="noreferrer">Enigma Rules version Mar 26, 2015</a>;
    } else if (prop === "Source_URL") {
        if (variant[prop].startsWith("http://hci-exlovd.hci.utah.edu")) {
            rowItem = <a target="_blank" href={variant[prop].split(',')[0]} rel="noreferrer">link to multifactorial analysis</a>;
        }
    } else if (prop === "Comment_on_clinical_significance_ENIGMA" || prop === "Clinical_significance_citations_ENIGMA") {
        const pubmed = "http://ncbi.nlm.nih.gov/pubmed/";
        rowItem = _.map(variant[prop].split(/PMID:? ?([0-9]+)/), (piece, idx) =>
            (/^[0-9]+$/.test(piece)) ? <a key={`pmid-${piece}-${idx}`} target="_blank" href={pubmed + piece} rel="noreferrer">PMID: {piece}</a> : piece);
    } else if (prop === "HGVS_cDNA") {
        rowItem = variant[prop].split(":")[1];
    } else if (prop === "HGVS_Protein") {
        rowItem = variant[prop].split(":")[1];
    } else if (/Allele_frequency_.*_ExAC/.test(prop)) {
        let count = variant[prop.replace("frequency", "count")],
            number = variant[prop.replace("frequency", "number")];
        rowItem = [variant[prop], <small key={`${prop}-meta`} style={{float: 'right'}}>({count} of {number})</small>];
    } else if (prop === "Allele_frequency_genome_GnomADv3" || prop === "Allele_frequency_exome_GnomAD") {
        let flag;
        if (prop === "Allele_frequency_genome_GnomADv3") {
            flag = variant.Flags_GnomADv3;
        } else {
            flag = variant.Flags_GnomAD;
        }
        if (!isEmptyField(flag)) {
            rowItem = [variant[prop], <small key={`${prop}-flag`} style={{float: 'right'}}><span className="fa fa-flag gnomad-flag"><span>{flag}</span></span></small>];
        } else {
            rowItem = normalizedFieldDisplay(variant[prop]);
        }
    } else if (/Allele_frequency_.*_GnomAD/.test(prop)) {
        let count = variant[prop.replace("frequency", "count")],
            number = variant[prop.replace("frequency", "number")],
            hom = variant[prop.replace("frequency", "count_hom")];
        rowItem = [variant[prop], <small key={`${prop}-meta`} style={{float: 'right'}}>({count} of {number}, Hom={hom})</small>];
    } else if (/count.*_GnomAD/.test(prop) || /number.*_GnomAD/.test(prop)) {
        rowItem = variant[prop];
    } else if (prop === "faf95_popmax_genome_GnomADv3") {
        rowItem = [variant[prop], <small key={`${prop}-meta`} style={{float: 'right'}}>({variant.faf95_popmax_population_genome_GnomADv3})</small>];
    } else if (prop === "Genomic_Coordinate_hg38" || prop === "Genomic_Coordinate_hg37") {
        let hgvs;
        if (prop === "Genomic_Coordinate_hg38") {
            hgvs = variant.Genomic_HGVS_38;
        } else if (prop === "Genomic_Coordinate_hg37") {
            hgvs = variant.Genomic_HGVS_37;
        }
        rowItem = generateLinkToGenomeBrowser(prop, variant[prop], hgvs);
    } else if (prop === "Synonyms") {
        let syns = variant[prop].split(',');
        let synsNoWhitespace = _.map(syns, s => s.replace(' ', '_'));
        rowItem = synsNoWhitespace.join(", ");
    } else {
        rowItem = normalizedFieldDisplay(variant[prop], prop);
    }

    if (_.contains(dateKeys, prop)) {
        rowItem = reformatDate(variant[prop]);
    }

    return rowItem;
}

// Finds sentence-ending boundaries in `str` -- a run of '.', '!', or '?'
// followed by whitespace or end-of-string -- and returns the string split
// into chunks that each end at one of those boundaries, with any trailing
// unterminated text returned as a final chunk.
//
// This scans for boundaries rather than matching whole sentences in one
// pattern, so it never drops characters that appear before the first
// boundary. A pattern like /[^.!?]*[.!?]+(?:\s+|$)/g looks reasonable but is
// unsafe here: scientific text is full of periods that aren't followed by
// whitespace (e.g. "p.Glu23ValfsTer17", "gnomAD v2.1", ">0.00002"). Because
// that pattern's `[^.!?]*` can't span across those internal periods, a match
// attempt starting before one of them fails outright, and String.match()
// silently skips ahead one character at a time looking for a position where
// a full match succeeds -- discarding everything it skipped, including the
// start of the sentence. Scanning for boundaries with regex.exec() and
// always slicing from the previous boundary avoids that: every character
// ends up in some chunk no matter how many non-terminating periods it
// contains.
function splitStringAtSentenceBoundaries(str) {
    const chunks = [];
    const boundaryRe = /[.!?]+(?=\s|$)/g;
    let lastIndex = 0;
    let match;

    while ((match = boundaryRe.exec(str)) !== null) {
        const end = match.index + match[0].length;
        chunks.push(str.slice(lastIndex, end));
        lastIndex = end;
    }

    if (lastIndex < str.length) {
        chunks.push(str.slice(lastIndex));
    }

    return chunks;
}

// Splits `content` (a plain string, a single React node, or an array mixing
// strings and React nodes -- e.g. the output of getFormattedFieldByProp for a
// field with embedded PMID links) into an array of "sentence groups".
//
// Each group is itself an array of parts (strings/nodes) that together make up
// one sentence, in original order. Non-string parts (like links) are attached
// to whichever sentence they appear in; they don't count as sentence
// boundaries themselves. This is used to power ExpandableText's truncation.
export function splitPartsIntoSentences(content) {
    const parts = Array.isArray(content) ? content : [content];
    const groups = [];
    let current = [];

    parts.forEach((part) => {
        if (typeof part !== "string") {
            // non-text parts (e.g. PMID links) belong to the sentence in progress
            current.push(part);
            return;
        }

        // split the string into chunks that each end at a sentence boundary
        // (a run of '.', '!', or '?' followed by whitespace or end-of-string),
        // with any trailing unterminated text captured as a final chunk
        const chunks = splitStringAtSentenceBoundaries(part);

        chunks.forEach((chunk) => {
            current.push(chunk);
            if (/[.!?]\s*$/.test(chunk)) {
                groups.push(current);
                current = [];
            }
        });
    });

    if (current.length > 0) {
        groups.push(current);
    }

    return groups;
}

// Like splitPartsIntoSentences, but splits into individual words instead.
// Non-string parts (e.g. links) count as a single "word" each.
export function splitPartsIntoWords(content) {
    const parts = Array.isArray(content) ? content : [content];
    const groups = [];

    parts.forEach((part) => {
        if (typeof part !== "string") {
            groups.push([part]);
            return;
        }

        const words = part.match(/\s*\S+/g) || [];
        words.forEach((word) => groups.push([word]));
    });

    return groups;
}

export function abbreviatedSubmitter(originalSubmitter) {
    return originalSubmitter
        .replace('Evidence-based Network for the Interpretation of Germline Mutant Alleles (ENIGMA)', 'ENIGMA')
        .replace('Breast Cancer Information Core (BIC)', 'BIC');
}

// For backward compatibility with code that uses require()
export default {
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
    capitalize,
    splitPartsIntoSentences,
    splitPartsIntoWords
};
