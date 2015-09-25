# Searching

Variants are identified using [HGVS nomenclature](http://www.hgvs.org/mutnomen/), a standard created by the Human Genome Variation Society.  Most genetic test results are reported in HGVS nomenclature. This nomenclature describes a variant by indicating (1) a reference sequence, (2) what kind of sequence it is (genomic, cDNA or protein), (3) the position of the variant in relation this reference sequence, and (4) the nucloetide or protein differences between the variant and the reference sequence.  For example, _NM_007294.3:c.5053A>G_ indicates a variant relative to reference sequence _NM_007294.3_, which is a cDNA sequence (as indicated by _c_), that the variant is in position 5053, and that it involves changing an _A_ to a _G_.

You can search for a variant by starting to enter its HGVS string.  As you type, the list of variants displayed will update automatically to show the variants that match your search string.

# Table Columns

#### Gene

The Gene column displays the name of the gene on which the variant was found.  This will be either BRCA1 or BRCA2.

#### HGVS

[HGVS](http://www.hgvs.org/mutnomen/) is a nomenclature proposed by the Human Genome Variation Society.  It describes a variant by specifying a reference sequence, the sequence type (genomic, cDNA or protein), the position of the variant within the reference sequence, and the changes that the variant introduces compared to the reference sequence.  For example, the string _NM_000059.3:c.8165C>G_ describes a variant relative to the reference sequence __NM_000059.3_, indicates that this sequence is a cDNA sequence (by the _c_), locates the variant at position 8165 in the reference sequence, and indicates that this variant involves changing a _C_ to a _G_.  On this site, we always use the same reference sequences: _NM_007294.3_ for BRCA1 and _NM_000059.3_ for BRCA2. For conciseness, we omit the name of the reference sequence from the HGVS string displayed.  So for the example of _NM_000059.3:c.8165C>G_, we would label the variant as _c.8165C>G_.

#### Pathogenicity

The Pathogenicity column indicates whether expert curators have determined if the variant is _pathogenic_ or _benign_.  _Pathogenic_ variants are potentially harmful, and have been shown to increase an individual's risk of disease.  _Benign_ variants are those that show no evidence of leading to disease risk, even after detailed examination.

#### Allele origin

Variants are classified as either _germline_ or _somatic_, depending on how they are acquired.  _Germline_ variants are genetic changes that we inherit from our parents.  _Somatic_ variants are DNA changes that we acquire over our lifetime, often through exposure to pollutants, toxins, radiation and other carcinogens.   

#### CVA

The variants shown on this site are stored in larger genomic variant databases.  One of the major variant datases is [ClinVar](http://www.ncbi.nlm.nih.gov/clinvar/).  This column gives the variant's ID in the ClinVar database, and contains hyperlinks showing the variant in ClinVar.
 
# Variant Detail Page Glossary

#### Gene symbol
This must be the preferred symbol from HGNC, if available. Separate multiple symbols with a semi-colon.

#### Reference sequence
Required for sequence variants if chromosome and some conbination of start and stop are not provided. The reference sequence for the HGVS expression provided in the next column, e.g. NM_000492.3, NG_016465.3, or NC_000007.13. 

#### HGVS_cDNA
"Required for sequence variants if start and stop are not provided. 
Provide the c. or g. portion of the nucleotide HGVS expression(s) for the variants being reported related to condition.  (http://www.hgvs.org/mutnomen/examplesDNA.html, http://www.hgvs.org/mutnomen/standards.html). Can also be used for a haplotype or compound heterozygotes, e.g.  c.[76A>C;105G>A]  or  c.[76A>C];[105G>A] . 

#### BIC Nomenclature
Optional.  A set of alternative,  common, or legacy names for the variant, or an alternate HGVS expression, e.g. for protein. Separate multiple names with a vertical bar (|).

#### Abbrev AA change
Optional. For variants that are assigned official allele names, e.g. CYP3A4*18 for one of the cytochrome P450 gene CYP3A4; or HLA-DRA*0102 for the MHC gene HLA-DRA.

#### Condition ID type
"Required if ""Preferred condition name"" is not supplied. 
Please specify the database for the ID in the Condition ID value column.   Allowed values are: OMIM, MeSH, MedGen, UMLS, Orphanet, HPO.  If you use a different source, please check with us first."

#### Condition ID value
"Required if ""Preferred condition name"" is not supplied. 
Please specify the identifier from the database referenced in the Condition ID type column. If there are multiple values, they all need to be from the same source, and separated by a semicolon. See http://www.ncbi.nlm.nih.gov/clinvar/docs/faq_submitters/#examples for examples of identifiers. Values may include both diseases and phenotypes (but not genes).

#### Condition category
Required but defaults to Disease if not filled in. Allowed spreadsheet values are Disease, DrugResponse, BloodGroup, Finding.  Use Finding for clinical findings/features.

#### Clinical significance
Required. Click cell for drop-down list; allowed values include Pathogenic, Likely pathogenic, Uncertain significance, Likely benign, Benign, association, drug response, confers sensitivity, protective, risk factor, other, not provided. If "other" or "drug response", please expand "Clinical Significance - more" columns to provide additional information.

#### Date last evaluated
Required, if available: date the clinical significance of the variant was last evaluated by the submitter (not the date the patient was evaluated in the clinic).  If not available, leave blank.  If  only month/year is provided we will convert it to the first day of the month.  If only year is provided it will be converted to the first day of the first month of the year. Please use the format  yyyy-mm-dd.

#### Assertion method
Optional; required for expert panel and professional society submissions. Free text name of the document in Assertion method citation, describing the method and criteria that were used to make the assertion of clinical significance. 

#### Assertion method citation
Optional; required for expert panel and professional society submissions; optional for others. Citation or URL describing the method and criteria used to make assertions of clinical significance.  

#### Clinical significance citations
Optional.  Citations documenting the clinical significance. Any of PubMed, PubMedCentral, DOI, NCBI Bookshelf combined with the id in that database (e.g. PMID:123456,  PMCID:PMC3385229, NBK:56955).  Separate multiple citations by a semicolon.

#### Citations or URLs for  clinical significance without database identifiers
Optional.  Citations that require a URL, or that do not have an identifier in one of the resources indicated in the Clinical significance citations column. Separate multiple citations/URLs by a vertical bar (|).

#### Comment on clinical significance
Optional, but highly encouraged.  Free text describing the rationale for the clinical significance.

#### Collection method
Required. Method used to collect the data that supports the assertion of clinical significance. Allowed values: case-control, clinical testing, literature only, reference population, research.

#### Allele origin
Required. The genetic origin of  the variant. Allowed values: biparental, de novo, germline, maternal, paternal, somatic, uniparental, unknown. Note that biparental and uniparental are intended for the context of uniparental disomy. If you'd like to indicate zygosity, please report counts of homozygotes and heterozygotes in columns BU-BW. 

#### Affected status
Required. Indicate whether or not the individual(s) had the condition. Allowed values: yes, no, unknown.

#### ClinVarAccession
Required for updates or if accession numbers were reserved. For new submissions, this accession number will be returned to you after processing of your submission. Provide the SCV for your submission (not the RCV).

#### Novel or Update
 "Optional. If you include ClinVar SCV accessions in the previous column, you must indicate whether the submission is novel (and accessions were reserved prior to submission), an update to an existing SCV record, or to delete an existing SCV record.
"
 

