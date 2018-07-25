The Variant Details page is displayed once you have chosen a single variant from a set of search results. This page is a “bird’s eye view” of the variant in question, based on publicly available data. This page will tell you all the standard aliases of the variant, its most definitive classification, and a brief history of variant interpretations. These details originate from large databases, rather than individual submissions from people who have undergone BRCA testing.

### Variant Details in the Expert Reviewed Data Portal

The Expert Reviewed Data Portal contains an abridged view of the variant data in BRCA Exchange, including a definitive, expertly reviewed ENIGMA interpretation. For all publicly available data on a variant, click the "Show All Public Data on this Variant" button at the bottom of the Variant Details Page.

#### Variant Nomenclature Tile

* #### Gene ((Gene_Symbol))
  The name of the gene on which the variant was found, as named by the HGNC. This will be either _BRCA1_ or _BRCA2_.

* #### HGVS Nomenclature
    A nomenclature system standardized by the Human Genome Variation Society \(HGVS\) which utilizes coordinates to name variants across genomic DNA, coding DNA, RNA, and protein molecules.

    Each coordinate provides an indicator of the biomolecule being referred to \('c.', 'g.', 'p.'\), a location in that biomolecule \('12345'\), and an indication of what changed in the variant \('G&gt;A', 'insA', 'dupT', 'del15', 'Asp...Asn'\). Most variant aliases utilize this standard nomenclature:
  * #### HGVS Nucleotide ((HGVS_cDNA))
    HGVS variant alias which references the nucleotide change based on the location in the coding DNA, not the genomic DNA.
  * #### HGVS RNA ((HGVS_RNA))
    The variant alias that uses RNA location.
  * #### HGVS Protein ((HGVS_Protein))
    The predicted protein-level change \(if any\) that would be introduced by the genomic variant.
  * #### Genomic Nomenclature GRCh38 ((Genomic_Coordinate_hg38))
    HGVS variant alias which references the nucleotide in genomic DNA, per the GRCh38 genome build. Click alias to access the UCSC Genome Browser.
  * #### Genomic Nomenclature GRCh37 ((Genomic_Coordinate_hg37))
    HGVS variant alias which references the nucleotide in genomic DNA, per the GRCh37 genome build. Click alias to access the UCSC Genome Browser.
  * For more information on HGVS nomenclature, visit the [HGVS site](http://varnomen.hgvs.org/bg-material/simple/).
* #### Transcript Identifier ((Reference_Sequence))
    These identifiers, which typically start in ‘NM’, are gene-specific, not variant-specific. They are generally used to keep track of different transcripts submitted to databases such as RefSeq, and are carried over to ClinVar. BRCA Exchange mainly uses two transcript identifiers. There is one main transcript used for _BRCA1_ \(NM 007294.3\) and another identifier used for _BRCA2_ \(NM 000059.3\). Other, less common transcripts are also present in the database.
* #### Abbreviated AA Change ((Protein_Change))
    A shortened abbreviation of the HGVS Protein alias using one-letter abbreviation for amino acids.
* #### BIC designation ((BIC_Nomenclature))
    A variant alias presented in BIC Nomenclature, which predates HGVS nomenclature and thus follows a different format.

#### Clinical Significance Tile

This tile organizes ENIGMA classifications and contextual information for the variant.

* #### ENIGMA interpretation ((Pathogenicity_expert))
    The most definitive, expertly reviewed interpretation for the variant according to [ENIGMA criteria](https://enigmaconsortium.org/wp-content/uploads/2017/12/ENIGMA_Rules_2017-06-29.pdf).
* #### IARC Class ((Clinical_significance_ENIGMA))
    Clinical Classification provided by the [International Agency for Research on Cancer](http://monographs.iarc.fr/ENG/Classification/).

* The rest of the fields in this this tile contain various citations, evidence, and other sources related to ENIGMA’s assessments

#### Previous Versions Tile

Here you will find a list of release dates and the interpretation as of that date. Clicking on the date takes you to the [release notes](http://brcaexchange.org/release/12), which provide you with summaries of changes, as well as a link to download the entire BRCA Exchange output as of that release date. This feature can be used to look at original variant interpretations alongside the data that informed them. These sets can be thought of as ‘variant time capsules,’ containing the data that was available when that release took place.

