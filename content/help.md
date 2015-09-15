# Searching

Variants are identified using [HGVS nomenclature](http://www.hgvs.org/mutnomen/), a standard created by the Human Genome Variation Society.  Most genetic test results are reported in HGVS nomenclature. This nomenclature describes a variant by indicating (1) a reference sequence, (s) what kind of sequence it is (genomic, cDNA or protein), (3) the position of the variant in this reference sequence, and (4) the nucloetide or protein differences between the variant and the reference sequence.  For example, _NM_007294.3:c.5053A>G_ indicates a variant relative to reference sequence _NM_007294.3_, which is a cDNA sequence (as indicated by the _c_), and that the variant is in position 5053 and involves changing an _A_ to a _G_.  You can search for a variant by starting to enter its HGVS string.  As you type, the list of variants displayed will update automatically to show the variants that match your search string.

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
 

