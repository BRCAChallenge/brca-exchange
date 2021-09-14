#### Functional Assay Scores

Functional assays are laboratory experiments that interrogate the impact of the variant on aspects of BRCA protein function or abundance.  To assess the functional effects of variants in *BRCA1* and *BRCA2*, assays have been developed to study the impact of the variant on experimental measurements that have been shown to have a high correlation with genetic and clinical data, including cell proliferation, cell survival, drug sensitivity and *Homology-Directed DNA Repair (HDR)*.

HDR, also called *homologous recombination (HR) mediated repair (HRR)*, is the major DNA repair pathway for error-free repair of DNA double-stranded breaks. The *BRCA1* and *BRCA2* proteins play a crucial role in in this pathway. A deficiency in HR has a severe impact on cell survival and leads to enhanced sensitivity to specific DNA damaging agents, including cisplatin or Poly (ADP-Ribose) polymerase (PARP).

For each functional assay publication, we present a set of general fields that are common to all variants in the publication (**Author**, **Publication**, **Previous Publications**, and **Report Descriptions**) as well as a set of fields that detail the reported experimental results on each each variant (this includes the **Reports** field with the authors’ report on overall functional impact, and additional fields reporting on specific experiments depending on the publication). We detail these fields below. For conciseness, we list the general fields just once.

----
#### General Fields

These fields are included for each publication, and are constant across all variants in the publication.

* #### Author ((Author))
    The lead author of the publication and the publication year.  This identifies the publication in a human-readable format.
* #### Publication ((Publication))
    PubMed ID and hyperlink to the publication in PubMed
* #### Previous Publications ((Previous Publications))
    This field lists any previous publications by the same lead or senior author. In general, these previous publications may provide foundational results that the current publication built on, with the current publication superseding earlier publications.
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
    Indicates whether the variant was designated as a neutral control, deleterious control, artificial variant, or VUS
* #### Report ((Result_Bouwman1_ENIGMA_BRCA12_Functional_Assays))
    Authors’ overall interpretation of the functional impact of the variant, taking into consideration both the cell proliferation and drug sensitivity experiments.
----
#### Bouwman et al, 2020

In this work, the authors tested the impact of 238 *BRCA1* variants on homologous recombination repair (HRR) and drug sensitivity.

* #### Cisplatin ((Cisplatin_Bouwman2_ENIGMA_BRCA12_Functional_Assays))
    Impact of the variant on Cisplatin drug sensitivity
* #### Olaparib ((Olaparib_Bouwman2_ENIGMA_BRCA12_Functional_Assays))
    Impact of the variant on Olaparib (PARP inhibitor) drug sensitivity
* #### DRGFP ((DRGFP_Bouwman2_ENIGMA_BRCA12_Functional_Assays))
     	Impact of the variant on homologous recombination repair measured in a direct-repeat GFP (DR-GFP) assay 
* #### Report ((Result_Bouwman2_ENIGMA_BRCA12_Functional_Assays))
    In this publication, the authors did not present an overall interpretation of the functional impact per variant, but presented the functional impact per variant per assay.
    
----
#### Fernandes et al, 2019

The authors of this paper assessed 99 *BRCA1* variants on their impact on transcriptional activation.

* #### Report ((Result_Fernandes_ENIGMA_BRCA12_Functional_Assays))
    Authors’ interpretation of the functional class of the variant using VarCall ([Iversen et al. 2011](https://pubmed.ncbi.nlm.nih.gov/21447777/)) (1-5)
----

#### Mesman et al, 2019

The authors evaluated 79 *BRCA2* variants on their impact on homology-directed DNA repair and drug sensitivity.

* #### Complementation ((Complementation_Mesman_ENIGMA_BRCA12_Functional_Assays))
    Indicates whether the variant was able to complement the cell lethal phenotype induced by loss of endogenous mouse *Brca2*
* #### HDR ((HDR_Mesman_ENIGMA_BRCA12_Functional_Assays))
    Impact of the variant on homology-directed DNA repair
* #### Cisplatin ((Cisplatin_Mesman_ENIGMA_BRCA12_Functional_Assays))
    Impact of the variant to Cisplatin drug sensitivity
* #### Report ((Result_Mesman_ENIGMA_BRCA12_Functional_Assays))
    Authors’ overall interpretation on the functional impact of the variant, considering complementation and homology directed repair
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
