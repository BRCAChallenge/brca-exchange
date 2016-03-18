# Searching

To search just start typing a variant in the search box and it will auto complete. Variants are identified using [HGVS nomenclature](http://www.hgvs.org/mutnomen/), a standard created by the Human Genome Variation Society.  Most genetic test results are reported in HGVS nomenclature. This nomenclature describes a variant by indicating (1) a reference sequence, (2) what kind of sequence it is (genomic, cDNA or protein), (3) the position of the variant in relation this reference sequence, and (4) the nucloetide or protein differences between the variant and the reference sequence.  For example, _NM_007294.3:c.5053A>G_ indicates a variant relative to reference sequence _NM_007294.3_, which is a cDNA sequence (as indicated by _c_), that the variant is in position 5053, and that it involves changing an _A_ to a _G_.

# Table Columns

#### Gene

The Gene column displays the name of the gene on which the variant was found.  This will be either BRCA1 or BRCA2.


#### Genomic (GRCh38)
Coordinate of the variant on the GRCh38 reference genome

#### Nucleotide
HGVS string that represents the variant at cDNA nucleotide level

#### Protein
The protein-level change (if any) that would be introduced by this variant.  HGVS notation indicates the position of the variant within the reference protein sequence, and the change that would be introduced by this variant.  For example, _p.(Tyr15His)_ indicates a change of _Tyr_ (Tyrosine) to _His_ (Histidine) at amino acid 15.  The notation _p.?_ indicates a variant is not within the protein-coding portions of the gene.

#### Pathogenicity
The Pathogenicity column indicates whether expert curators have determined if the variant is pathogenic or benign.

##### _What do these classifications mean?_
- *Pathogenic variants* confer an increased risk of disease.
- *Likely pathogenic variants* have good evidence to support an association with disease risk.
- *Likely benign variants* have good evidence to support no association with disease risk.
- *Benign variants* are not associated with any markedly increased risk of disease.
- *Variants of uncertain significance (VUS)* are those for which the evidence of disease risk is not clear yet, sometimes because there is not yet enough evidence to classify them as either pathogenic or benign.

##### _My variant is classified pathogenic. What do I do now?_
This website is not intended to provide a clinical diagnosis. Please contact your primary care provider to determine what steps may be necessary.


#### Source URL(s)
URL(s) pointing back to the original source data


# Variant Detail Page Glossary

#### Gene symbol
Gene name, as named by [HGNC](http://www.genenames.org/).

#### Genomic Coordinate
Coordinate of the variant on the genome.

#### Reference sequence
The reference sequence corresponding to the HGVS coordinates indicated in the HGVS_cDNA column, e.g. NM_000492.3, NG_016465.3, or NC_000007.13. 

#### HGVS_cDNA
[HGVS coordinates](#HGVS) on the cDNA.

#### BIC Nomenclature
BIC name for this variant.

#### Abbrev AA Change
The amino acid change introduced by the variant, using abbreviated (one-letter) amino acid names.

#### URL
A link with further information on this variant.

#### Condition ID type
The Condition ID type, together with the Condition ID value, indicate the condition associated with this variant according to a standard ontology. Condition ID types can be: OMIM, MeSH, MedGen, UMLS, Orphanet, HPO.

#### Condition ID value
The Condition ID value, together with the Condition ID type, indicate the condition associated with this variant according to a standard ontology.  The condition ID type specifies an ontology, while the condition ID value indicates the term within that ontology that identifes this condition.

#### Condition category
Human-readable condition, describing the biological impact of the variation.

#### Clinical significance
Possible values are Pathogenic, Likely pathogenic, Uncertain significance, Likely benign, Benign, association, drug response, confers sensitivity, protective, risk factor, other, not provided.

#### Date last evaluated
The date on which the clinical significance of the variant was last evaluated by the submitter.

#### Assertion method
Describes the method and criteria that were used to make the assertion of clinical significance. 

#### Assertion method citation
Citation or URL describing the method and criteria used to make assertions of clinical significance.  

#### Clinical significance citations
Citations documenting the clinical significance. Can be from PubMed, PubMedCentral, DOI, or NCBI Bookshelf.

#### Citations or URLs for clinical significance without database identifiers
Citations that require a URL, or that do not have an identifier in one of the resources indicated in the Clinical significance citations column.

#### Comment on clinical significance
Describes the rationale for the clinical significance.

#### Collection method
Method used to collect the data that supports the assertion of clinical significance. Allowed values: case-control, clinical testing, literature only, reference population, research.

#### Allele origin
Variants are classified as either _germline_ or _somatic_, depending on how they are acquired.  _Germline_ variants are genetic changes that we inherit from our parents.  _Somatic_ variants are DNA changes that we acquire over our lifetime, often through exposure to pollutants, toxins, radiation and other carcinogens.  

#### Affected status
Indicates whether or not the individual(s) had the condition.

#### ClinVarAccession
Accession number in the ClinVar database.

#### Novel or Update
Indicates whether the submission to ClinVar is novel (and accessions were reserved prior to submission), an update to an existing record.

#### HGVS_protein
[HGVS coordinates](#HGVS) on the protein.


# Lollipop Plots

Lollipop plots are a tool to visualize the chromosomal position and pathogenicity classification for each variant in a gene.  Here, each circle-capped 'lollipop' indicates whether a BRCA1 and BRCA2 variant is benign (labeled in blue) or pathogenic (labeled in red).  The diameter of the circle on top of each lollipop indicates the number of unique variants that have been observed at an individual locus.  The bottom bar that runs along the x-axis of the diagram displays the position of each exon in the selected gene.  To alternate between the BRCA1 and BRCA2 lollipop charts, click the *Select Gene* button located in the upper-left hand corner of the page.
