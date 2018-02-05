# Help

## How do I search for a variant?

Welcome to the BRCA Exchange, where you can search BRCA1 and BRCA2 variants in our expertly curated and maintained data portal. This site allows you to access important, up-to-date information about a BRCA variant, such as its clinical significance.
To begin, notice the search bar on the home page. This search bar accesses our databases directly, and can handle any identifying names for the variant or variants you are searching for. Here are some examples of what can be typed in this search box:

* BRCA1 or BRCA2 \(Gene Symbol\)
* BRCA1 c.1105G&gt;A \(Gene + HGVS Nucleotide\)
* c.1105G&gt;A \(HGVS Nucleotide\)
* chr17:g.4305831 \(Genomic Nomenclature\)
* IVS19-1179G&gt;T, 2043G&gt;C \(BIC Designation\)
* NM\_007924.3 \(Transcript Identifier\)
* _p.\(Pro1238Leu\)_ \(HGVS Protein\)
* _P1238L_ \(Abbreviated Amino Acid Change\)

Clicking the magnifying glass will execute a search for the variants that fit your search criteria. To learn more about HGVS nomenclature, visit the HGVS Nomenclature portion of the Reference Guide \(found below\). A table of all matching variants is provided once you have used the search box. You can sort search results alphabetically by clicking any column header once. The list will be alphabetized based on the column you clicked. Clicking once more will sort search results in reverse-alphabetical order.
When you have identified your variant of interest in the list, you can click anywhere in the variant’s data row to access the Variant Detail Page. The Variant Detail Page provides an organized summary of the searched variant’s information, including its aliases, its clinical significance, the date it was last evaluated, and other relevant data. Information is grouped into tiles for your convenience.

Please note that searching for a variant using a genomic coordinate will return all variants that match according to any of the hg38, hg19, or hg18 genome builds. For example, if you search for a variant using its hg38 coordinate, and it happens to match the coordinate of some variant  in the hg19 build, both variants will be returned in the search. In this case, make sure to verify that the coordinate\(s\) and genome build are correct once you navigate to the Variant Details Page.

## What do the fields in the Variant Details Page mean?

The Variant Details page is displayed once you have chosen a single variant from a set of search results. This page is a “bird’s eye view” of the variant in question, based on publicly available data. This page will tell you all the standard aliases of the variant, its most definitive classification, and a brief history of variant interpretations. These details originate from large databases, rather than individual submissions from people who have undergone BRCA testing.

### Variant Details in the Expert Reviewed Data Portal

The Expert Reviewed Data Portal contains an abridged view of the variant data in BRCA Exchange, including a definitive, expertly reviewed ENIGMA interpretation. For all publicly available data on a variant, click the "Show All Public Data on this Variant" button at the bottom of the Variant Details Page.

#### Variant Nomenclature Tile

* _Gene_: The name of the gene on which the variant was found, as named by the HGNC. This will be either _BRCA1_ or _BRCA2_.
* _HGVS Nomenclature_: A nomenclature system standardized by the Human Genome Variation Society \(HGVS\) which utilizes coordinates to name variants across genomic DNA, coding DNA, RNA, and protein molecules. Each coordinate provides an indicator of the biomolecule being referred to \('c.', 'g.', 'p.'\), a location in that biomolecule \('12345'\), and an indication of what changed in the variant \('G&gt;A', 'insA', 'dupT', 'del15', 'Asp...Asn'\). Most variant aliases utilize this standard nomenclature:
  * _HGVS Nucleotide:_ HGVS variant alias which references the nucleotide change based on the location in the coding DNA, not the genomic DNA.
  * _HGVS RNA_: The variant alias that uses RNA location.
  * _HGVS Protein_: The predicted protein-level change \(if any\) that would be introduced by the genomic variant.
  * _Genomic Nomenclature GRCh38_: HGVS variant alias which references the nucleotide in genomic DNA, per the GRCh38 genome build. Click alias to access the UCSC Genome Browser.
  * _Genomic Nomenclature GRCh37_: HGVS variant alias which references the nucleotide in genomic DNA, per the GRCh37 genome build. Click alias to access the UCSC Genome Browser.
  * For more information on HGVS nomenclature, visit the [HGVS site](http://varnomen.hgvs.org/bg-material/simple/).
* _Transcript Identifier_: These identifiers, which typically start in ‘NM’, are gene-specific, not variant-specific. They are generally used to keep track of different transcripts submitted to databases such as RefSeq, and are carried over to ClinVar. BRCA Exchange mainly uses contains two transcript identifiers. There is one main transcript used for BRCA1 \(NM 007294.3\) and another identifier used for BRCA2 \(NM 000059.3\). Other, less common transcripts are also present in the database.
* _Abbreviated AA Change_: A shortened abbreviation of the HGVS Protein alias using one-letter abbreviation for amino acids.
* _BIC designation_: A variant alias presented in BIC Nomenclature, which predates HGVS nomenclature and thus follows a different format.

#### Clinical Significance Tile

This tile organizes ENIGMA classifications and contextual information for the variant.

* _ENIGMA interpretation_: The most definitive, expertly reviewed interpretation for the variant according to [ENIGMA criteria](https://enigmaconsortium.org/wp-content/uploads/2017/12/ENIGMA_Rules_2017-06-29.pdf).
* _IARC Class_: Clinical Classification provided by the [International Agency for Research on Cancer](http://monographs.iarc.fr/ENG/Classification/).
* The rest of the fields in this this tile contain various citations, evidence, and other sources related to ENIGMA’s assessments

#### Previous Versions Tile

Here you will find a list of release dates and the interpretation as of that date. Clicking on the date takes you to the [release notes](http://brcaexchange.org/release/12), which provide you with summaries of changes, as well as a link to download the entire BRCA Exchange output as of that release date. This feature can be used to look at original variant interpretations alongside the data that informed them. These sets can be thought of as ‘variant time capsules,’ containing the data that was available when that release took place.

