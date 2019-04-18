# GENOMIC VERSION:
VERSION = "hg38" # equivalent to GRCh38

# Specific columns in the output matrix
COLUMN_SOURCE = 0
COLUMN_GENE = 1
COLUMN_GENOMIC_HGVS = 2
COLUMN_VCF_CHR = 3
COLUMN_VCF_POS = 4
COLUMN_VCF_REF = 5
COLUMN_VCF_ALT = 6

# This is the string to be stored when a field is empty
DEFAULT_CONTENTS = "-"

# key value pair dictionaries of all extra fields in various databases to add
GENOME1K_FIELDS = {"Allele_frequency": "AF",
                   "EAS_Allele_frequency": "EAS_AF",
                   "EUR_Allele_frequency": "EUR_AF",
                   "AFR_Allele_frequency": "AFR_AF",
                   "AMR_Allele_frequency": "AMR_AF",
                   "SAS_Allele_frequency": "SAS_AF",
                   "BX_ID": "BX_ID"}

CLINVAR_FIELDS = {"HGVS": "HGVS",
                  "Submitter": "Submitter",
                  "Clinical_Significance": "ClinicalSignificance",
                  "Date_Last_Updated": "DateLastUpdated",
                  "DateSignificanceLastEvaluated":
                      "DateSignificanceLastEvaluated",
                  "SCV": "SCV",
                  "SCV_Version": "SCV_Version",
                  "Allele_Origin": "Origin",
                  "Protein": "Protein",
                  "Method": "Method",
                  "Description": "Description",
                  "Review_Status": "ReviewStatus",
                  "Synonyms": "Synonyms",
                  "Summary_Evidence": "SummaryEvidence",
                  "BX_ID": "BX_ID"}

LOVD_FIELDS = {"Variant_frequency": "frequency",
               "Functional_analysis_technique": "functional_analysis_technique",
               "Functional_analysis_result": "functional_analysis_result",
               "HGVS_cDNA": "cDNA",
               "HGVS_protein": "Protein",
               "Genetic_origin": "genetic_origin",
               "RNA": "RNA",
               "Variant_effect": "variant_effect",
               "Individuals": "individuals",
               "Submitters": "submitters",
               "DBID": "DBID",
               "BX_ID": "BX_ID",
               "Created_date": "created_date",
               "Edited_date": "edited_date",
               "Submission_ID": "submission_id"
               }

EX_LOVD_FIELDS = {"Combined_prior_probablility": "combined_prior_p",
                  "Segregation_LR": "segregation_lr",
                  "Sum_family_LR": "sum_family_lr",
                  "Co_occurrence_LR": "co_occurrence_lr",
                  "Missense_analysis_prior_probability": "missense_analysis_prior_p",
                  "Posterior_probability": "posterior_p",
                  "IARC_class": "iarc_class",
                  "BIC_Nomenclature": "bic_dna_change",
                  "Literature_source": "observational_reference",
                  "HGVS_cDNA": "dna_change",
                  "HGVS_protein": "protein_change",
                  "BX_ID": "BX_ID"}

BIC_FIELDS = {"Clinical_classification": "Category",
              "Number_of_family_member_carrying_mutation": "Number_Reported",
              "Patient_nationality": "Nationality",
              "Germline_or_Somatic": "G_or_S",
              "Mutation_type": "Mutation_Type",
              "BIC_Designation": "Designation",
              "Clinical_importance": "Clinically_Importance",
              "Ethnicity": "Ethnicity",
              "Literature_citation": "Reference",
              "BX_ID": "BX_ID"}

ESP_FIELDS = {"polyPhen2_result": "PH",
              "Minor_allele_frequency_percent": "MAF",
              "EA_Allele_Frequency": "BX_EAAF",
              "AA_Allele_Frequency": "BX_AAAF",
              "Allele_Frequency": "BX_AF",
              "BX_ID": "BX_ID"}


EXAC_FIELDS = {"Allele_frequency": "AF",
               "Allele_count_AFR": "AC_AFR",
               "Allele_number_AFR": "AN_AFR",
               "Allele_frequency_AFR": "AF_AFR",
               "Homozygous_count_AFR": "Hom_AFR",
               "Allele_count_AMR": "AC_AMR",
               "Allele_number_AMR": "AN_AMR",
               "Allele_frequency_AMR": "AF_AMR",
               "Homozygous_count_AMR": "Hom_AMR",
               "Allele_count_EAS": "AC_EAS",
               "Allele_number_EAS": "AN_EAS",
               "Allele_frequency_EAS": "AF_EAS",
               "Homozygous_count_EAS": "Hom_EAS",
               "Allele_count_FIN": "AC_FIN",
               "Allele_number_FIN": "AN_FIN",
               "Allele_frequency_FIN": "AF_FIN",
               "Homozygous_count_FIN": "Hom_FIN",
               "Allele_count_NFE": "AC_NFE",
               "Allele_number_NFE": "AN_NFE",
               "Allele_frequency_NFE": "AF_NFE",
               "Homozygous_count_NFE": "Hom_NFE",
               "Allele_count_OTH": "AC_OTH",
               "Allele_number_OTH": "AN_OTH",
               "Allele_frequency_OTH": "AF_OTH",
               "Homozygous_count_OTH": "Hom_OTH",
               "Allele_count_SAS": "AC_SAS",
               "Allele_number_SAS": "AN_SAS",
               "Allele_frequency_SAS": "AF_SAS",
               "Homozygous_count_SAS": "Hom_SAS",
               "BX_ID": "BX_ID"}

EXAC_SUBPOPULATIONS = ["AFR", "AMR", "EAS", "FIN", "NFE", "OTH", "SAS"]

GNOMAD_FIELDS = {"HGVS_cDNA": "hgvsc",
                 "HGVS_protein": "hgvsp",
                 "Allele_count": "ac",
                 "Allele_number": "an",
                 "Allele_frequency": "af",
                 "Allele_frequency_AFR": "AFR_ac",
                 "Allele_count_hom_AFR": "AFR_ac_hom",
                 "Allele_count_hemi_AFR": "AFR_ac_hemi",
                 "Allele_frequency_AFR": "AFR_af",
                 "Allele_number_AFR": "AFR_an",
                 "Allele_frequency_AMR": "AMR_ac",
                 "Allele_count_hom_AMR": "AMR_ac_hom",
                 "Allele_count_hemi_AMR": "AMR_ac_hemi",
                 "Allele_number_AMR": "AMR_an",
                 "Allele_frequency_AMR": "AMR_af",
                 "Allele_frequency_ASJ": "ASJ_ac",
                 "Allele_count_hom_ASJ": "ASJ_ac_hom",
                 "Allele_count_hemi_ASJ": "ASJ_ac_hemi",
                 "Allele_number_ASJ": "ASJ_an",
                 "Allele_frequency_ASJ": "ASJ_af",
                 "Allele_frequency_EAS": "EAS_ac",
                 "Allele_count_hom_EAS": "EAS_ac_hom",
                 "Allele_count_hemi_EAS": "EAS_ac_hemi",
                 "Allele_number_EAS": "EAS_an",
                 "Allele_frequency_EAS": "EAS_af",
                 "Allele_frequency_FIN": "FIN_ac",
                 "Allele_count_hom_FIN": "FIN_ac_hom",
                 "Allele_count_hemi_FIN": "FIN_ac_hemi",
                 "Allele_number_FIN": "FIN_an",
                 "Allele_frequency_FIN": "FIN_af",
                 "Allele_frequency_NFE": "NFE_ac",
                 "Allele_count_hom_NFE": "NFE_ac_hom",
                 "Allele_count_hemi_NFE": "NFE_ac_hemi",
                 "Allele_number_NFE": "NFE_an",
                 "Allele_frequency_NFE": "NFE_af",
                 "Allele_frequency_OTH": "OTH_ac",
                 "Allele_count_hom_OTH": "OTH_ac_hom",
                 "Allele_count_hemi_OTH": "OTH_ac_hemi",
                 "Allele_number_OTH": "OTH_an",
                 "Allele_frequency_OTH": "OTH_af",
                 "Allele_frequency_SAS": "SAS_ac",
                 "Allele_count_hom_SAS": "SAS_ac_hom",
                 "Allele_count_hemi_SAS": "SAS_ac_hemi",
                 "Allele_number_SAS": "SAS_an",
                 "Allele_frequency_SAS": "SAS_af",
                 "BX_ID": "BX_ID"}

GNOMAD_SUBPOPULATIONS = ["AFR", "AMR", "ASJ", "EAS", "FIN", "NFE", "OTH", "SAS"]

FINDLAY_BRCA1_RING_FUNCTION_SCORES_FIELDS = {
               "HGVS_Nucleotide": "hgvs_nucleotide",
               "Log_RNA_Depletion": "log_rna_depletion",
               "Functional_Enrichment_Score": "functional_enrichment_score",
               "BX_ID": "BX_ID"
               }

FIELD_DICT = {"1000_Genomes": GENOME1K_FIELDS,
              "ClinVar": CLINVAR_FIELDS,
              "LOVD": LOVD_FIELDS,
              "exLOVD": EX_LOVD_FIELDS,
              "ExAC": EXAC_FIELDS,
              "ESP": ESP_FIELDS,
              "BIC": BIC_FIELDS,
              "GnomAD": GNOMAD_FIELDS,
              "Findlay_BRCA1_Ring_Function_Scores": FINDLAY_BRCA1_RING_FUNCTION_SCORES_FIELDS}

LIST_TYPE_FIELDS = {
    "SCV", # Clinvar, treating it as list, to have the same order as with SCV_Version
    "SCV_Version"
}

ENIGMA_FILE = "enigma_from_clinvar.tsv"

GENOME1K_FILE = "1000G_brca.sorted.hg38.vcf"
CLINVAR_FILE = "ClinVarBrca.vcf"
LOVD_FILE = "sharedLOVD_brca12.sorted.hg38.vcf"
EX_LOVD_FILE = "exLOVD_brca12.sorted.hg38.vcf"
BIC_FILE = "bic_brca12.sorted.hg38.vcf"
EXAC_FILE = "exac.brca12.sorted.hg38.vcf"
ESP_FILE = "esp.sorted.hg38.vcf"
GNOMAD_FILE = "gnomAD.clean.sorted.hg38.vcf"

# Functional Assays
FINDLAY_BRCA1_RING_FUNCTION_SCORES_FIELDS_FILE = "findlay_BRCA1_ring_function_scores.clean.sorted.hg38.vcf"

