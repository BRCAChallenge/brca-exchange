# Help
The BRCA Exchange is a curated data platform that provides information on variants of specific genes: _BRCA1_ and _BRCA2_. Rather than contributing new information about a variant, the BRCA Exchange becomes a single source that collects and organizes existing information. The BRCA Exchange retrieves data from a variety of databases, such as ClinVar and LOVD. Using these sources, BRCA Exchange answers three basic questions about a variant:

* What is the most definitive clinical interpretation available for this variant?
* When was this interpretation made?
* What publicly available data informed this interpretation?

There are two data views used by BRCA Exchange to aid researchers and clinicians who work with _BRCA1_ and _BRCA2._ The default, or Expert Reviewed, view of BRCA Exchange will only display variant aliases, ENIGMA variant interpretations, and a brief interpretation history. The ENIGMA Consortium assesses variant pathogenicity using an expertly developed set of _BRCA1_ and _BRCA2_ classification criteria. The _All Public_ view contains all publicly available data on a variant. Thus, this view of the data provides a more complete, but also more complicated, characterization of a variant. Using multiple large-scale sources and two different portals, BRCA Exchange serves as a contextualized resource for curation of BRCA variant data.

## Table of Contents
* How do I search for a variant?
* Variant Details Page
	* Variant Nomenclature
	* Clinical Significance Tiles
	* Transcript Visualization
	* Multifactorial Likelihood Analysis
	* Allele Frequency Reference Sets
* Lollipop Plots
* Downloading Variant Data
* How do I use Filters?
* How do I use Column Selectors?

## How do I search for a variant?

In the BRCA Exchange, you can search _BRCA1_ and _BRCA2_ variants in our expertly curated and maintained data portal. This site allows you to access important, up-to-date information about a BRCA variant, such as its clinical significance. To begin, notice the search bar on the home page. This search bar accesses our databases directly, and can handle any identifying names for the variant or variants you are searching for. Here are some examples of what can be typed in this search box:
* _BRCA1_ or _BRCA2_ \(Gene Symbol\)
* _BRCA1_ c.1105G&gt;A \(Gene + HGVS Nucleotide\)
* _c.1105G&gt;A_ \(HGVS Nucleotide\)
* _chr17:g.4305831_ \(Genomic Nomenclature\)
* _IVS19-1179G&gt;T, 2043G&gt;C_ \(BIC Designation\)
* _NM\_007924.3_ \(Transcript Identifier\)
* _p.\(Pro1238Leu\)_ \(HGVS Protein\)
* _P1238L_ \(Abbreviated Amino Acid Change\)

Clicking the magnifying glass will execute a search for the variants that fit your search criteria. A table of all matching variants is provided once you have used the search box. You can sort search results alphabetically by clicking any column header once. The list will be alphabetized based on the column you clicked. Clicking once more will sort search results in reverse-alphabetical order. When you have identified your variant of interest in the list, you can click anywhere in the variant’s data row to access the Variant Detail Page. The Variant Detail Page provides an organized summary of the searched variant’s information, including its aliases, its clinical significance, the date it was last evaluated, and other relevant data. Information is grouped into tiles for your convenience.

Please note that searching for a variant using a genomic coordinate will return all variants that match according to any of the hg38, hg19, or hg18 genome builds. For example, if you search for a variant using its hg38 coordinate, and it happens to match the coordinate of some variant in the hg19 build, both variants will be returned in the search. In this case, make sure to verify that the coordinate(s) and genome build are correct once you navigate to the Variant Details Page.

# Variant Details Page

## Variant Nomenclature
The Variant Nomenclature fields contain forms of identification by which the variant can be labeled.

* #### Gene Symbol
	The Gene Symbol field displays the name of the gene on which the variant was found, as named by HGNC. This will be either _BRCA1_ or _BRCA2_.

* #### Reference cDNA Sequence
	This field provides the reference sequence on which the variant was observed, corresponding to the HGVS coordinates indicated in the Nucleotide field e.g. NM\_000492.3, NG\_016465.3, or NC\_000007.13.

* #### HGVS Nucleotide
	This is the HGVS variant alias which references the nucleotide change based on the location in the coding DNA, not the genomic DNA.

* #### HGVS Protein
	This alias utilizes predicted protein-level change \(if any\) that would be introduced by the genomic variant.

* #### BIC Variant Identifier
	This variant alias is presented in BIC Nomenclature, which predates HGVS nomenclature and thus follows a different format.

* #### Genome \(GRCh38\)
	This alias utilizes the coordinate of the variant on the GRCh38 reference genome. Click alias to access the UCSC Genome Browser.

* #### Genome \(GRCh37\)
	This alias utilizes the coordinate of the variant on the GRCh37 reference genome. Click alias to access the UCSC Genome Browser.

* #### Genome \(GRCh36\)
	This alias utilizes the coordinate of the variant on the GRCh36 reference genome. Click alias to access the UCSC Genome Browser.

## Clinical Significance Tiles

### ENIGMA
* #### Clinical Significance
	* This field displays the variant’s clinical classification and supporting information as provided by ENIGMA expert panel review.
	##### What do these classifications mean?
	* Pathogenic variants confer a markedly increased risk of disease.
	* Likely pathogenic variants have good evidence to support an association with markedly increased disease risk.
	* Likely benign variants have good evidence to support no association with markedly increased disease risk.
	* Benign variants are not associated with any markedly increased risk of disease.
	* Variants of uncertain significance \(VUS\) are those for which the evidence of disease risk is not clear yet, sometimes because there is not yet enough evidence to classify them as either pathogenic or benign.

	Markedly increased disease risk for a variant in _BRCA1_ and _BRCA2_ is currently defined as risk that is sufficient to alter patient management based on detection of that single genetic variant, after considering other non-genetic risk factors for the individual.

* #### Comment on Clinical Significance
	This field comments on the derivation of the IARC class, provided by ENIGMA.

* #### Assertion Method
	This field provides the citation or URL describing the method and criteria used to make assertions of clinical significance.

* #### Date last evaluated
	This is a display of the date on which the clinical significance of the variant was last evaluated by ENIGMA.

* #### Collection Method
	This is the method used to collect the data that supports the assertion of clinical significance. The following values may be listed: case-control, clinical testing, literature only, reference population, research.

* #### Clinical Significance Citation
	This field provides citations that document the clinical significance. They may be from PubMed, PubMedCentral, DOI, or NCBI Bookshelf.

* #### Allele Origin
	Variants are classified as either germline or somatic, depending on how they are acquired. Germline variants are genetic changes that we inherit from our parents. Somatic variants are DNA changes that we acquire over our lifetime, often through exposure to pollutants, toxins, radiation and other carcinogens.

### ClinVar
##### _Updates Coming Soon_

* #### Clinical Significance
	These fields are extracted from the [ClinVar website](https://www.ncbi.nlm.nih.gov/clinvar/) relating to variant clinical classification, for individual submitters to ClinVar \(includes repeated information from ENIGMA and BIC tiles\).

	##### What do these classifications mean?
	The below classifications will be found in ClinVar's "clinical significance" field. They are implemented per [recommendations by ACMG/AMP](https://www.ncbi.nlm.nih.gov/pubmed/25741868) for variants interpreted by Mendelian disorders.

	* Pathogenic
	* Likely pathogenic
	* Likely benign
	* Benign
	* Uncertain Significance

* #### Submitter
	This field provides the name of the submitting organization. From ClinVar.
* #### SCV Accession
	This is the SCV accession from the ClinVar repository.
* #### Date last updated \(ClinVar\)
	The field that shows the date the variant was last submitted to ClinVar. Note that this date is not equivalent to the date last evaluated.
* ##### _Coming Soon:_
	Updated Classification descriptions, Information about tile functionality, resources for star review, definition of summary evidence.

For more information, please visit the [ClinVar website](https://www.ncbi.nlm.nih.gov/clinvar/) or the [ClinVar Data Dictionary](https://www.ncbi.nlm.nih.gov/projects/clinvar/ClinVarDataDictionary.pdf).

### LOVD

##### _Updates Coming Soon_

* #### Clinical Significance
	These fields are extracted from [LOVD](http://www.lovd.nl/3.0/home) and are relevant to variant clinical classification.
* ##### _Coming Soon:_
	Submitters, Variant Effect, Genetic Origin, Submission ID, Variant Haplotype

	For more information, please visit the [LOVD website](http://www.lovd.nl/3.0/home).

### BIC
* #### Patient Nationality
	The field that indicates nationality of the patient associated with the variant.
* #### Ethnicity
	The field that provides the ethnicities of the observed patients.
* #### Family members carrying this variant
	Number of family members carrying this variant.
* #### Literature Reference
	Literature Reference\(s\). From BIC.
* #### Allele Origin
	Origin of the allele, Germline \(G\) or Somatic \(S\). Variants are classified as either germline or somatic, depending on how they are acquired. Germline variants are genetic changes that we inherit from our parents. Somatic variants are DNA changes that we acquire over our lifetime, often through exposure to pollutants, toxins, radiation and other carcinogens. From BIC

## Transcript Visualization

The transcript visualization depicts the reference transcript of the gene \(either _BRCA1_ or _BRCA2_\) on which the selected variant is found. The gene is drawn as a series of exons with intervening introns. Note that exons are drawn scaled \(nonlinearly\) to their actual size, to make the smaller exons easier to see, while introns are drawn at a fixed width.

The gene is shown in full at the top of the visualization, with the exons/introns in which the variant lies \(i.e., the area of interest\) highlighted in yellow. The area of interest is shown zoomed-in below the full depiction, but contains the same information.

The variant's modifications to the reference transcript are depicted in three colored regions:
* Substitution variants \(depicted in lime green\)
* Deletions \(depicted in red with cross-hatched black lines\)
* Insertions \(depicted in blue\)

The visualization also displays other regions of interest for the gene, such as donor and acceptor sites and clinically important functional domains. The visibility of these regions can be toggled individually via the checkboxes below the transcript.

Clinically important functional domains are regions of the gene that ENIGMA has assessed to be relevant to risk. More specifically, these domains are regions of the gene that ENIGMA has determined to be relevant to risk when modified by missense alterations or in-frame deletions that could encode a stable \(non-truncated\) protein. This determination is made using existing knowledge of pathogenicity for individual variants in these domains. Clicking the label for a given clinically important domain will display the corresponding domain in an animated outline. Clicking again on the label will disable the outline, as will clicking another domain designated as clinically important.

## [Multifactorial Likelihood Analysis Classification](https://www.ncbi.nlm.nih.gov/pubmed/21990134)
Component likelihoods and final classification are from multifactorial likelihood analysis, which is a quantitative integrated evaluation of variant pathogenicity.
* #### Posterior Probability of pathogenicity \(ExUV\)
	This value is the posterior probability of pathogenicity from ExUV.
* #### Prior probability of pathogenicity \(ExUV\)
	This value is the combined prior probability in favor of pathogenicity. It is a combination of the missense analysis prior probability and the splicing analysis prior probability. Generally, it is the higher of these two prior probabilities. From ExUV.
* #### Missense analysis probability of pathogenicity \(ExUV\)
	This prior probability estimate combines position in the protein with an evaluation of missense substitutions that fall in the proteins key functional domains. From ExUV.
* #### Co-occurrence likelihood \(ExUV\)
	This value is the likelihood ratio based on the frequency of co-occurrence between the variant of interest and clearly pathogenic variants in the same gene. From ExUV.
* #### Segregation Likelihood Ratio \(ExUV\)
	This value is the likelihood ratio based on segregation analysis
* #### Summary Family History Likelihood Ratio \(ExUV\)
	This value is the likelihood ratio based on an analysis of the severity of summary family histories of breast and/ or ovarian cancer.

## Allele Frequency Reference Sets
The Allele Frequency fields describe how often the variant was observed overall, or in specific populations.
* #### Allele Frequency \(ExAC minus TCGA\)
	Minor allele frequency from ExAC data that does not include TCGA data
* #### Allele Frequency \(1000 Genomes\)
	Overall allele frequency, from 1000 Genomes
* #### AFR Allele Frequency \(1000 Genomes\)
	African-American minor allele frequency, from 1000 Genomes
* #### AMR Allele Frequency \(1000 Genomes\)
	Allele frequency in the Admixed American population, from 1000 Genomes
* #### EAS Allele Frequency \(1000 Genomes\)
	Allele frequency in the East Asian population, from 1000 Genomes
* #### EUR Allele Frequency \(1000 Genomes\)
	Allele frequency in the European population, from 1000 Genomes
* #### SAS Allele Frequency \(1000 Genomes\)
	Allele frequency in the Southeast Asian population, from 1000 Genomes
* #### Allele Frequencies: EA|AA|All \(ESP\)
	Allele frequencies from ESP, expressed as EA (European)|AA (African American)|All
* #### Bar Charts
	* ##### _Coming Soon_:
		Supplementary information about bar chart details and functionality


## Lollipop Plots
Lollipop plots are a tool to visualize the chromosomal position and pathogenicity classification for each variant in a gene. Here, each circle-capped ‘lollipop’ indicates whether a BRCA1 and BRCA2 variant is pathogenic (labeled in red), benign (labeled in light-blue), or the variants clinical significance is uncertain (labeled in dark-blue). The y-axis represents the pathogenicity classification of a variant. The x-axis represents the genomic coordinates of those variants in GRCh38 human reference space. The bottom box that runs along the x-axis of the diagram displays the position of each exon in the selected gene.

To open the lollipop chart click on the ‘Show Lollipop Chart’ button located at the top of the Variants page. To alternate between the BRCA1 and BRCA2 lollipop charts, click on the tab that contains the relevant gene name located in the upper-left hand corner next to the lollipop chart. To select a specific area of the gene to zoom in on a region, simply click and drag on the lower box to create a shaded box which displays which region you have selected. You can also select for specific exon ranges by clicking on the relevant colored exon box on the lower box-chart. Mouse scrolling on the main chart itself will also zoom in and out of a region. To de-select or zoom all the way back out, you can either click on an unshaded space in the lower box or click the ‘Hide Lollipop Chart’ button at the top twice to refresh and reset the chart. On the legend located within the chart you can see the relative percentages of pathogenicity classifications of variants that match the search and filter criteria for a particular gene. To hide/dim a class of variants set of lollipops, click on the colored box on the legend. You can reset this by clicking on the colored box again.

<div style="max-width: 1000px; margin-left: auto; margin-right: auto; height: 0px; position: relative; padding-bottom: 59.552%;"><iframe src="https://streamable.com/e/zx9c?muted=1&amp;autoplay=1&amp;hd=1" frameborder="0" allowfullscreen webkitallowfullscreen mozallowfullscreen scrolling="no" style="width: 100%; height: 100%; position: absolute;"></iframe></div>

* * * * * * * * * * *


## Downloading Variant Data
To download the variant data click on the ‘Download’ button located above the search bar within the Variants Table page. When you click on it the data that is downloaded is a comma delimited .csv file containing the set of variant details for the variants that matched the search and/or filtering criteria. The first row in the file represent the column labels.

## How do I use Column Selectors?
Columns allow you to filter the variant table (in search results) according to fields found in the tiles on the Variant Details Page. Column selectors correspond to each tile’s fields, whose descriptions are found above. These selectors allow you to include any tile’s field as its own column in your search results. You can then sort results using added columns.

##### _Updates Coming Soon_

## How do I use Filters?
##### _Updates Coming Soon_
