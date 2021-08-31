Functional assays are laboratory experiments that interrogate the impact of the variant on some aspect of BRCA protein function or abundance.

One of the most important functions of the BRCA1 and BRCA2 proteins is *Homologous DNA Repair (HDR)*, and HDR is the focus of many BRCA functional assays.  HDR is a molecular DNA repair process that enables the cell  to repair double-strand DNA breaks under certain circumstances.   Loss of HDR activity leads to an accumulation of DNA damage, which in turn can lead to cancer in humans and early cell death in the laboratory.  

*Cell survival* assays measure HDR activity by subjecting tissue with a given genetic variant to mutational processes to incur DNA damage, and measuring the rate of survival of the cells in that tissue, with a greater rate of cell death suggesting a harmful variant with greater impact.

*Drug sensitivity* assays measure HDR activity by exposing the cell to PARP inhibitors are a family of pharmaceutical drugs that inhibit the cell’s ability to repair single-strand DNA breaks.  These single-stranded DNA breaks lead to double-stranded DNA breaks, which lead to early cell death if not repaired.  For these reasons, the presence of harmful BRCA variants makes tissue sensitive to PARP inhibition, leading to early cell death.  PARP inhibitors include the drugs Olaparib, Rucaparib, Niraparib and Talazoparib.

The cellular processes that impact BRCA protein abundance often impact *RNA production*, or the abundance of BRCA mRNA.  Harmful variants can reduce the abundance of BRCA mRNA by impacting cellular processes including transcription, mRNA splicing and others.  When a variant impacts the abundance of BRCA RNA, that will impact the abundance of BRCA protein, and reduce the overall BRCA functional activity.  

For each functional assay publication, we present a set of general fields that are common to all variants in the publication (**Author**, **Publication**, **Previous Publications**, and **Report Descriptions**) as well as a set of fields that detail the reported experimental results on each each variant (this includes the **Reports** field with the authors’ report on overall functional impact, and additional fields reporting on specific experiments depending on the publication). We detail these fields below. For conciseness, we list the general fields just once.

----
#### General Fields

These fields are included for each functional assay publication, and are constant across all variants in the publication.

* #### Author ((Author))
    The lead author of the publication and the publication year.  This identifies the publication in a human-readable format.
* #### Publication ((Publication))
    PubMed ID and hyperlink to the publication in PubMed
* #### Previous Publications ((Previous Publications))
    This field lists any previous publications by the same authors. In general, these previous publications may provide foundational results that the current publication built on, with the current publication superseding earlier publications.
* #### Report Descriptions ((Report Descriptions))
    This field describes the content of the Report field by presenting the set of possible values.

----
#### Findlay et al, 2018

This MAVE assay (Multiplexed Assay of Variant Effects) interrogated roughly 3900 *BRCA1* variants, comprising most possible missense variants in the 13 exons that comprise the RING and BRCT protein domains. The authors measured the impact of the variants on cell survival and RNA production.

* #### Functional Enrichment Score ((Functional_Enrichment_Findlay_ENIGMA_BRCA12_Functional_Assays))
    Numerical score reflecting the impact of the variant on cell survivial
* #### RNA Score ((RNA_Score_Findlay_ENIGMA_BRCA12_Functional_Assays))
    Numerical score reflecting the impact of the variant on RNA production
* #### RNA Class ((RNA_Class_Findlay_ENIGMA_BRCA12_Functional_Assays))
    Authors' qualitative assessment of the impact of the variant on RNA production
* #### Report ((Result_Findlay_ENIGMA_BRCA12_Functional_Assays))
    Authors’ reported interpretation of the functional impact of the variant.
----
#### Starita et al, 2018

In this publication, the authors performed a multiplex reporter assay to assess the impact of 1056 *BRCA1* variants on homologous DNA repair.

* #### Report ((Result_Starita_ENIGMA_BRCA12_Functional_Assays))
    Authors' report on the functional impact of the variant given the results of the reporter assay
----
#### Petitalot et al, 2019

In this paper, the authors assessed the impact of 78 *BRCA1* variants on homology-detected DNA repair, cellular localization, protein stability and peptide binding.

* #### Control Group ((Control_Group_Petitalot_ENIGMA_BRCA12_Functional_Assays))
    Indicates whether the variant was part of a control group of either causal or neutral variants
* #### Report ((Result_Petitalot_ENIGMA_BRCA12_Functional_Assays))
    Authors' overall report on the functional impact of the variant, taking into consideration all of the individual assay results.
----
#### Bouwman et al, 2013

In this publication, the authors tested the impact of 86 *BRCA1* variants on cell proliferation and drug sensitivity.


* #### Selection ((Selection_Bouwman1_ENIGMA_BRCA12_Functional_Assays))
    Indicates whether the variant was designated as a neutral control, deleterious control, artificial, or VUS
* #### Report ((Result_Bouwman1_ENIGMA_BRCA12_Functional_Assays))
    Authors’ overall interpretation of the functional impact of the variant, taking into consideration both the cell proliferation and drug sensitivity experiments.
----
#### Bouwman et al, 2020

In this work, the authors tested the impact of 238 *BRCA1* variants on homologous DNA repair and drug sensitivity.

* #### Cisplatin ((Cisplatin_Bouwman2_ENIGMA_BRCA12_Functional_Assays))
    Impact of the variant on Cisplatin drug sensitivity
* #### Olaparib ((Olaparib_Bouwman2_ENIGMA_BRCA12_Functional_Assays))
    Impact of the variant on Olaparib drug sensitivity
* #### DRGFP ((DRGFP_Bouwman2_ENIGMA_BRCA12_Functional_Assays))
    Impact of the variant in a direct-repeat GFP homologous recombination DNA repair
* #### Report ((Result_Bouwman2_ENIGMA_BRCA12_Functional_Assays))
    In this publication, the authors did not report an overall interpretation of the functional impact of the variants, but presented many sets of assay results which each offer one measurement of functional impact.
----
#### Fernandes et al, 2019

The authors of this paper assessed 99 *BRCA1* variants on their impact on transcriptional activation in a yeast reporter system.

* #### Class ((Class_Fernandes_ENIGMA_BRCA12_Functional_Assays))
    Authors' interpretation of the functional class of the variant (1-5)
* #### Report ((Result_Fernandes_ENIGMA_BRCA12_Functional_Assays))
    Authors' report on the overall functional impact of the variant given the results of the yeast reporter system.
----

#### Mesman et al, 2019

The authors evaluated 79 *BRCA2* variants on their impact on homologous DNA repair and drug sensitivity.

* #### Complementation ((Complementation_Mesman_ENIGMA_BRCA12_Functional_Assays))
    Indicates whether the variant was able to complement the phenotype induced by loss of endogenous *BRCA2*
* #### HDR ((HDR_Mesman_ENIGMA_BRCA12_Functional_Assays))
    Impact of the variant on homolgous DNA repair
* #### Cisplatin ((Cisplatin_Mesman_ENIGMA_BRCA12_Functional_Assays))
    Impact of the variant to Cisplatin drug sensitivity
* #### Report ((Result_Mesman_ENIGMA_BRCA12_Functional_Assays))
    Authors' report on the functional impact of the variant, considering both homologous DNA repair and drug sensitivity
----

#### Richardson et al, 2021

In this publication, the authors evaluated the impact of 252 variants in the *BRCA2* DNA-binding domain on homology-directed DNA repair

* #### HDR ((HDR_Richardson_ENIGMA_BRCA12_Functional_Assays))
    Numeric measurement of the impact of the variant on homolous DNA repair
* #### Report ((Result_Richardson_ENIGMA_BRCA12_Functional_Assays))
    Authors' qualitative interpretation of the numeric functional impact of the variant
----

#### Ikegami et al, 2020

In this publication, the authors evaluated the impact of 239 *BRCA2* variants on homology-directed DNA repair and drug sensitivity.

* #### Discordant ((Discordant_Ikegami_ENIGMA_BRCA12_Functional_Assays))
    Indicates if the functional evidence on the variant was discordant with its clinical classification
* #### Olaparib ((Olaparib_fClass_Ikegami_ENIGMA_BRCA12_Functional_Assays))
    Impact of the variant on Olaparib drug sensitivity
* #### Niraparif ((Niraparif_fClass_Ikegami_ENIGMA_BRCA12_Functional_Assays))
    Impact of the variant on Nirapaif drug sensitivty
* #### Rucaparib ((Rucaparib_fClass_Ikegami_ENIGMA_BRCA12_Functional_Assays))
    Impact of the variant on Rucaparib drug sensitivity
* #### CBDCA ((CBDCA_fClass_Ikegami_ENIGMA_BRCA12_Functional_Assays))
    Impact of the variant on carboplatin (CBDCA) drug sensitivity
* #### Report ((Result_Ikegami_ENIGMA_BRCA12_Functional_Assays))
    Authors' report on the functional impact of the variant, taking into consideration all of the assays performed
----
#### Biwas et al, 2020

In this paper, the authors characterized the impact of 88 *BRCA2* variants on cell survival and drug sensitivity.

* #### Cell Survival ((Cell_Survival_Biwas_ENIGMA_BRCA12_Functional_Assays))
    Impact of the variant on cell survival in the HAT medium
* #### Drug Sensitivity ((Drug_Sensitivity_Biwas_ENIGMA_BRCA12_Functional_Assays))
    Impact of the variant on drug sensitivity
* #### HAT_DS Score ((HAT_DS_Score_Biwas_ENIGMA_BRCA12_Functional_Assays))
    Numeric score reflecting both cell survival and drug sensitivity analyzed together
* #### Report ((Result_Biwas_ENIGMA_BRCA12_Functional_Assays))
    Authors' qualitative interpretation of the functional impact of the variant, given the numeric score
