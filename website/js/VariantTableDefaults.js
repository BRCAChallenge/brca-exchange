/*global module: false, require: false, URL: false, Blob: false */
'use strict';

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

module.exports = {
    defaultExpertColumns: defaultExpertColumns,
    defaultResearchColumns: defaultResearchColumns,
    allSources: allSources
};
