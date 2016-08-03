# Searching

To search just start typing a variant in the search box and it will auto complete. Variants are identified using [HGVS nomenclature](http://www.hgvs.org/mutnomen/), a standard created by the Human Genome Variation Society.  Most genetic test results are reported in HGVS nomenclature. This nomenclature describes a variant by indicating (1) a reference sequence, (2) what kind of sequence it is (genomic, cDNA or protein), (3) the position of the variant in relation this reference sequence, and (4) the nucloetide or protein differences between the variant and the reference sequence.  For example, _NM_007294.3:c.5053A>G_ indicates a variant relative to reference sequence _NM_007294.3_, which is a cDNA sequence (as indicated by _c_), that the variant is in position 5053, and that it involves changing an _A_ to a _G_.

# Understanding the data

#### Gene

The Gene column displays the name of the gene on which the variant was found,
as named by [HGNC](http://www.genenames.org/).  This will be either BRCA1 or 
BRCA2.


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

