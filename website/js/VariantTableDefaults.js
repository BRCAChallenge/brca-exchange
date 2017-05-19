/*global module: false, require: false, URL: false, Blob: false */
'use strict';
var _ = require('underscore');


const defaultExpertColumns = ['Gene_Symbol', 'HGVS_cDNA', 'HGVS_Protein', 'Protein_Change', 'BIC_Nomenclature', 'Pathogenicity_expert'];
const defaultResearchColumns = ['Gene_Symbol', 'Genomic_Coordinate_hg38', 'HGVS_cDNA', 'HGVS_Protein', 'Pathogenicity_all', 'Allele_Frequency'];

const allSources = {
    "Variant_in_ENIGMA": 1,
    "Variant_in_ClinVar": 1,
    "Variant_in_1000_Genomes": 1,
    "Variant_in_ExAC": 1,
    "Variant_in_LOVD": 1,
    "Variant_in_BIC": 1,
    "Variant_in_ESP": 1,
    "Variant_in_exLOVD": 1
};

const getDefaultExpertColumns = function() {
    // Returns clone of default expert columns.
    return defaultExpertColumns.slice(0);
};

const getDefaultResearchColumns = function() {
    // Returns clone of default research columns.
    return defaultResearchColumns.slice(0);
};

const getAllSources = function() {
    // Returns a clone of all sources object.
    return _.clone(allSources);
};

module.exports = {
    getDefaultExpertColumns: getDefaultExpertColumns,
    getDefaultResearchColumns: getDefaultResearchColumns,
    getAllSources: getAllSources
};
