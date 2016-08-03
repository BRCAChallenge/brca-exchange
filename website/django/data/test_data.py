def existing_variant():
    return {"Variant_frequency_LOVD": "-", "Genomic_Coordinate_hg37": "NM_007294.3:chr17:41255102:A>G",
            "Literature_source_exLOVD": "-", "ClinVarAccession_ENIGMA": "SCV000244790", "Discordant": "Concordant",
            "Condition_category_ENIGMA": "Disease", "Date_Last_Updated_ClinVar": "2015-01-12",
            "Method_ClinVar": "germline", "AFR_Allele_frequency_1000_Genomes": "0.1566",
            "EUR_Allele_frequency_1000_Genomes": "0.3559", "Variant_in_LOVD": False, "Segregation_LR_exLOVD": "-",
            "Source_URL": "http://www.ncbi.nlm.nih.gov/clinvar/?term=SCV000244790",
            "Source": "ENIGMA,ClinVar,1000_Genomes", "Max_Allele_Frequency": "0.498000 (SAS from 1000 Genomes)",
            "Condition_ID_value_ENIGMA": "BREAST-OVARIAN CANCER, FAMILIAL, SUSCEPTIBILITY TO, 1; BROVCA1 (604370)",
            "Protein_Change": "", "HGVS_Protein": "p.?", "Co_occurrence_LR_exLOVD": "-",
            "SAS_Allele_frequency_1000_Genomes": "0.498", "Variant_in_exLOVD": False,
            "Allele_frequency_1000_Genomes": "0.335264", "Reference_Sequence": "NM_007294.3",
            "Clinical_Significance_ClinVar": "Benign", "Functional_analysis_result_LOVD": "-",
            "AMR_Allele_frequency_1000_Genomes": "0.3646", "Variant_in_ESP": False, "Variant_in_BIC": False,
            "Clinical_significance_ENIGMA": "Benign", "Synonyms": "", "Clinical_classification_BIC": "-",
            "Gene_Symbol": "BRCA1",
            "Comment_on_clinical_significance_ENIGMA": "Class 1 not pathogenic based on frequency >1% in an outbred sampleset. Frequency 0.3304 (Asian), 0.128 (African), 0.3575 (European), derived from 1000 genomes (2012-04-30).",
            "Sum_family_LR_exLOVD": "-", "Assertion_method_ENIGMA": "ENIGMA BRCA1/2 Classification Criteria (2015)",
            "Clinical_significance_citations_ENIGMA": "", "Variant_in_ENIGMA": True,
            "Posterior_probability_exLOVD": "-", "SIFT_VEP": "-", "Germline_or_Somatic_BIC": "-",
            "Missense_analysis_prior_probability_exLOVD": "-",
            "Genomic_Coordinate_hg36": "NM_007294.3:chr17:38508628:A>G", "BIC_Identifier": "IVS7+1037T>C",
            "Allele_Frequency": "0.335264 (1000 Genomes)", "Number_of_family_member_carrying_mutation_BIC": "-",
            "Pathogenicity_research": "Benign(ENIGMA); Benign (ClinVar)", "Allele_Origin_ClinVar": "SCV000244790",
            "Submitter_ClinVar": "curation", "Genomic_Coordinate_hg38": "NM_007294.3:chr17:43103085:A>G",
            "Patient_nationality_BIC": "-", "Literature_citation_BIC": "-",
            "HGVS_cDNA": "NM_007294.3(BRCA1):c.441+1037T>C", "Functional_analysis_technique_LOVD": "-",
            "Collection_method_ENIGMA": "Curation", "Allele_frequency_ExAC": "-", "Mutation_type_BIC": "-",
            "Combined_prior_probablility_exLOVD": "-", "Minor_allele_frequency_ESP": "-", "IARC_class_exLOVD": "-",
            "Clinical_importance_BIC": "-",
            "Assertion_method_citation_ENIGMA": "http://enigmaconsortium.org/documents/ENIGMA_Rules_2015-03-26.pdf",
            "EAS_Allele_frequency_1000_Genomes": "0.371", "Variant_haplotype_LOVD": "-",
            "SCV_ClinVar": "Evidence-based_Network_for_the_Interpretation_of_Germline_Mutant_Alleles_(ENIGMA)",
            "Condition_ID_type_ENIGMA": "OMIM", "Ethnicity_BIC": "-", "PolyPhen_VEP": "-",
            "Date_last_evaluated_ENIGMA": "1/12/15", "Pathogenicity_default": "Benign / Little Clinical Significance",
            "Variant_in_ExAC": False, "URL_ENIGMA": "", "Variant_in_ClinVar": True, "Allele_origin_ENIGMA": "Germline",
            "Origin_of_variant_LOVD": "-", "Variant_in_1000_Genomes": True}


def new_variant():
    v = existing_variant()
    v["Genomic_Coordinate_hg38"] = "chr17:999999:A>G"
    return v
