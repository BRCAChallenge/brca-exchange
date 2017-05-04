# Searching

To search just start typing a variant in the search box and it will auto complete. Variants are identified using [HGVS nomenclature](http://www.hgvs.org/mutnomen/), a standard created by the Human Genome Variation Society.  Most genetic test results are reported in HGVS nomenclature. This nomenclature describes a variant by indicating (1) a reference sequence, (2) what kind of sequence it is (genomic, cDNA or protein), (3) the position of the variant in relation this reference sequence, and (4) the nucloetide or protein differences between the variant and the reference sequence.  For example, _NM_007294.3:c.5053A>G_ indicates a variant relative to reference sequence _NM_007294.3_, which is a cDNA sequence (as indicated by _c_), that the variant is in position 5053, and that it involves changing an _A_ to a _G_.

* * * * * * * * * * *

# Column Groups

#### Clinical Significance (ENIGMA)
Variant clinical classification and supporting information as provided by ENIGMA expert panel review.

#### Clinical Significance (ClinVar)
Fields extracted from ClinVar ([website](https://www.ncbi.nlm.nih.gov/clinvar/)) relating to variant clinical classification, for individual submitters to ClinVar (excludes ENIGMA and BIC - see other panels on the Variant Details page).

#### Clinical Significance (BIC)
Fields extracted from BIC ([website](https://research.nhgri.nih.gov/bic/)) relevant to variant clinical classification.

#### Clinical Significance (LOVD)
Fields extracted from LOVD ([website](http://www.lovd.nl/3.0/home)) relevant to variant clinical classification.

#### [Multifactorial Likelihood Analysis Classification](https://www.ncbi.nlm.nih.gov/pubmed/21990134)
Component likelihoods and final classification from multifactorial likelihood analysis, a quantitative integrated evaluation of variant pathogenicity. 

* * * * * * * * * * *

# Table Columns

### Variant Nomenclature

The Variant Nomenclature fields are forms of identification by which the variant can be labeled.

#### BIC Variant Identifier
Variant identifier in BIC nomenclature

#### Nucleotide
[HGVS coordinates](#HGVS) on the cDNA.

#### Protein
[HGVS coordinates](#HGVS) on the protein.

#### Reference cDNA Sequence
The reference sequence on which the variant was observed, corresponding to the HGVS coordinates indicated in the Nuceotide 
column, e.g. NM_000492.3, NG_016465.3, or NC_000007.13. 

#### SCV Accession (ClinVar)
SCV accession in the ClinVar repository

## Origin

The Origin fields describe how the variants are acquired or the samples in which they were observed

#### Allele Origin (BIC)
Origin of the allele, Germline (G) or Somatic (S).  Variants are classified as either _germline_ or _somatic_, depending on how they are acquired.  _Germline_ variants are genetic changes that we inherit from our parents.  _Somatic_ variants are DNA changes that we acquire over our lifetime, often through exposure to pollutants, toxins, radiation and other carcinogens.  
From BIC

#### Allele Origin (ClinVar)
Variants are classified as either _germline_ or _somatic_, depending on how they are acquired.  _Germline_ variants are genetic changes that we inherit from our parents.  _Somatic_ variants are DNA changes that we acquire over our lifetime, often through exposure to pollutants, toxins, radiation and other carcinogens.  

####  Allele Origin (ENIGMA)
Variants are classified as either _germline_ or _somatic_, depending on how they are acquired.  _Germline_ variants are genetic changes that we inherit from our parents.  _Somatic_ variants are DNA changes that we acquire over our lifetime, often through exposure to pollutants, toxins, radiation and other carcinogens.  

#### Allele Origin (LOVD)
Allele origin, from LOVD. Variants are classified as either _germline_ or _somatic_, depending on how they are acquired.  _Germline_ variants are genetic changes that we inherit from our parents.  _Somatic_ variants are DNA changes that we acquire over our lifetime, often through exposure to pollutants, toxins, radiation and other carcinogens.  

####  Ethnicity (BIC)
Ethnicities of the observed patients, from BIC

#### Patient Nationality (BIC)
Nationality of the patient associated with the variant.  From BIC

#### Variant Haplotype (BIC)
Haplotype (example: FA FANCD1_00037)

## Frequency

The Frequency columns describe how often the variant was observed overall, or in specific populations.

#### African Allele Frequency (1000 Genomes)

African-American minor allele frequency, from 1000 Genomes

#### Allele Frequency
Allele frequency, estimated from the combination of ExAC, ESP and 1000 Genomes allele frequences

#### Allele Frequency (1000 Genomes)
Overall allele frequency, from 1000 Genomes

#### Allele Frequencies: EA|AA|All (ESP)
Allele frequencies from ESP, expressed as EA (European)|AA (African American)|All

#### Allele Frequency (ExAC)
Minor allele frequency

#### AMR Allele Frequency (1000 Genomes)
Allele frequency in the Admixed American population, from 1000 Genomes

#### EAS Allele Frequency (1000 Genomes)
Allele frequency in the East Asian population, from 1000 Genomes

#### EUR Allele Frequency (1000 Genomes)
European allele frequency in the European population, from 1000 Genomes

#### Maximum Allele Frequency
Highest allele frequency reported by any data source, for any population.

#### SAS Allele Frequency (1000 Genomes)
Allele frequency in the Southeast Asian population, from 1000 Genomes

#### Variant Frequency (LOVD)
Allele frequency, from LOVD, expressed as a fraction relative to the reference population

#### Allele Frequency in reference sets
Variant allele frequency in non-cancer sample sets of different ethnicities. The term "minus TCGA" indicates that information from The Cancer Genome Atlas dataset was excluded when determining frequency Maximum allele frequency is the maximum allele frequency observed in any single ethnic group.

## Genomic

The Genomic columns describe genomic traits of the variant, such as the location of the variant in a reference genome


#### Gene Symbol
The Gene Symbol column displays the name of the gene on which the variant was found,
as named by [HGNC](http://www.genenames.org/).  This will be either _BRCA1_ or
_BRCA2_.

#### Genome (GRCh38)
Coordinate of the variant on the GRCh38 reference genome

#### Genome (GRCh36)
Coordinate of the variant on the GRCh36 reference genome

#### Genome (GRCh37)
Coordinate of the variant on the GRCh37 reference genome


## Bioinformatic Annotation

This set of columns offers bioinformatic analyses of the variant.

#### Mutation category (BIC)

Transcript context of the mutation.

#### PolyPhen score

PolyPhen score.  A score closer to 1 suggests a damaging variant

#### SIFT score

SIFT score, indicating if the variant is tolerated, deleterious or other

## Probability

The Probability columns describe the likelihood that the variant is pathogenic, given the variant data, and possibly given
a specific category of data on the variant.

#### Co-occurrence likelihood (exLOVD)
The likelihood ratio based on the frequency of co-occurrence between the variant of interest and clearly
 pathogenic variants in the same gene. From exLOVD.

#### Missense analysis probability of pathogenicity (exLOVD)
This prior probability estimate combines position in the protein with an evaluation of missense substitutions
 that fall in the proteins key functional domains.  From exLOVD.

#### Prior probability of pathogenicity (exLOVD)
The combined prior probability in favor of pathogenicity. is a combination of the missense analysis prior
 probability and the splicing analysis prior probability. Generally, it is the higher of these two prior probabilities.  From exLOVD.

#### Probability of pathogenicity (exLOVD)
Posterior probability of pathogenicity. From exLOVD.

#### Segregation Likelihood Ratio (exLOVD)
The likelihood ratio based on segregation analysis

#### Summary Family History Likelihood Ratio (exLOVD)
The likelihood ratio based on an analysis of the severity of summary family histories of breast and/ or ovarian
 cancer.

## Significance

The Significance columns describe the clinical significance of the variant.

#### Analysis Method (ClinVar)
Method by which the significance of the variant was determined.  From ClinVar

#### Assertion Method (ENIGMA)
Citation or URL describing the method and criteria used to make assertions of clinical significance.  

#### Clinical Importance (BIC)
Clinical significance (yes|no|pending|unknown)

#### Clinical Significance (BIC)
Clinical significance, from BIC.  This ranges from Class 1 to Class 5, where Class 1 is Benign and Class 5 is
 Pathogenic.

#### Clinical Significance (ClinVar)
The Clinical Significance columns indicate whether expert curators have determined if the variant is pathogenic or benign.

##### _What do these classifications mean?_
- *Pathogenic variants* confer an increased risk of disease.
- *Likely pathogenic variants* have good evidence to support an association with disease risk.
- *Likely benign variants* have good evidence to support no association with disease risk.
- *Benign variants* are not associated with any markedly increased risk of disease.
- *Variants of uncertain significance (VUS)* are those for which the evidence of disease risk is not clear yet, sometimes because there is not yet enough evidence to classify them as either pathogenic or benign.


#### Clinical Significance (ENIGMA)
The Clinical Significance columns indicates whether expert curators have determined if the variant is pathogenic or benign.

##### _What do these classifications mean?_
- *Pathogenic variants* confer an increased risk of disease.
- *Likely pathogenic variants* have good evidence to support an association with disease risk.
- *Likely benign variants* have good evidence to support no association with disease risk.
- *Benign variants* are not associated with any markedly increased risk of disease.
- *Variants of uncertain significance (VUS)* are those for which the evidence of disease risk is not clear yet, sometimes because there is not yet enough evidence to classify them as either pathogenic or benign.


#### Collection Method (ENIGMA)
Method used to collect the data that supports the assertion of clinical significance. Allowed values: case-control, clinical testing, literature only, reference population, research.

#### Comment on Clinical Significance (ENIGMA)
Comment on how the derivation of the IARC class, from ENIGMA

#### Date last evaluated (ENIGMA)
The date on which the clinical significance of the variant was last evaluated by ENIGMA

#### Date last updated (ClinVar)
Date the variant was last updated in ClinVar

#### Has Discordant Evidence
Indicates if there is evidence of a discordant classification for this variant

#### Functional Analysis Result (LOVD)
Functional analysis or prediction of the impact of this variant.  From LOVD

#### Functional Analysis Method (LOVD)
The method used in analyzing variant function.  From LOVD.

#### Pathogenicity
Pathogenicity indicates whether expert curators have determined if the variant is pathogenic or benign.  This column can contain
multiple values, if multiple curators have assessed this variant.

##### _What do these classifications mean?_
- *Pathogenic variants* confer an increased risk of disease.
- *Likely pathogenic variants* have good evidence to support an association with disease risk.
- *Likely benign variants* have good evidence to support no association with disease risk.
- *Benign variants* are not associated with any markedly increased risk of disease.
- *Variants of uncertain significance (VUS)* are those for which the evidence of disease risk is not clear yet, sometimes because there is not yet enough evidence to classify them as either pathogenic or benign.

## Pedigree

Pedigree columns describe the pedigree or ancestry of patients with this variant.

#### Family members carrying this variant (BIC)
Number of family members carrying this variant.  From BIC.

## Publications

Further information on this variant, from the scientific literature


#### Assertion Method (ENIGMA)
Describes the method and criteria that were used to make the assertion of clinical significance. 

#### Clinical Significance Citation (ENIGMA)
Citations documenting the clinical significance. Can be from PubMed, PubMedCentral, DOI, or NCBI Bookshelf.

#### Citations or URLs for clinical significance without database identifiers
Citations that require a URL, or that do not have an identifier in one of the resources indicated in the Clinical significance citations column.

#### Literature Reference (BIC)
Literature Reference(s).  From BIC

#### Literature Reference (exLOVD)
Literature Reference(s).  From exLOVD.


## Source

The Source columns indicate specific sources of data on this variant

#### ClinVar Accession
ClinVar SCV accession

#### Condition Category (ENIGMA)
Human-readable condition, describing the biological impact of the variation.

#### Condition ID Type (ENIGMA)
The Condition ID type, together with the Condition ID value, indicate the condition associated with this variant according to a standard ontology. Condition ID types can be: OMIM, MeSH, MedGen, UMLS, Orphanet, HPO.

#### Condition ID Value (ENIGMA)
The Condition ID value, together with the Condition ID type, indicate the condition associated with this variant according to a standard ontology.  The condition ID type specifies an ontology, while the condition ID value indicates the term within that ontology that identifes this condition.

#### Source(s)
List of repositories containing the variant

#### Source URL(s)
URL(s) pointing back to the original source data, with further information on this variant.

#### Submitter (ClinVar)
Submitting organization.  From ClinVar.

#### URL (ENIGMA)
URL listing the variant.  From ENIGMA.

* * * * * * * * * * *

# Downloading Variant Data

To download the variant data click on the 'Download' button located above the search bar within the Variants Table page. When you click on it the data that is downloaded is a comma delimited .csv file containing the set of variant details for the variants that matched the search and/or filtering criteria. The first row in the file represent the column labels.

* * * * * * * * * * *

# Lollipop Plots

Lollipop plots are a tool to visualize the chromosomal position and pathogenicity classification for each variant in a gene.  Here, each circle-capped 'lollipop' indicates whether a _BRCA1_ and _BRCA2_ variant is pathogenic (labeled in red), benign (labeled in light-blue), or the variants clinical significance is uncertain (labeled in dark-blue). The y-axis represents the pathogenicity classification of a variant. The x-axis represents the genomic coordinates of those variants in GRCh38 human reference space. The bottom box that runs along the x-axis of the diagram displays the position of each exon in the selected gene.

To open the lollipop chart click on the 'Show Lollipop Chart' button located at the top of the Variants page. To alternate between the BRCA1 and BRCA2 lollipop charts, click on the tab that contains the relevant gene name located in the upper-left hand corner next to the lollipop chart. To select a specific area of the gene to zoom in on a region, simply click and drag on the lower box to create a shaded box which displays which region you have selected. You can also select for specific exon ranges by clicking on the relevant colored exon box on the lower box-chart. Mouse scrolling on the main chart itself will also zoom in and out of a region. To de-select or zoom all the way back out, you can either click on an unshaded space in the lower box or click the 'Hide Lollipop Chart' button at the top twice to refresh and reset the chart. On the legend located within the chart you can see the relative percentages of pathogenicity classifications of variants that match the search and filter criteria for a particular gene. To hide/dim a class of variants set of lollipops, click on the colored box on the legend. You can reset this by clicking on the colored box again.

<div style="width: 100%; height: 0px; position: relative; padding-bottom: 59.552%;"><iframe src="https://streamable.com/e/zx9c?muted=1&amp;autoplay=1&amp;hd=1" frameborder="0" allowfullscreen webkitallowfullscreen mozallowfullscreen scrolling="no" style="width: 100%; height: 100%; position: absolute;"></iframe></div>
