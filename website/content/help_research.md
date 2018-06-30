Help
============================

The BRCA Exchange is a curated data platform that provides information on variants of specific genes: _BRCA1_ and _BRCA2_. Rather than contributing new information about a variant, the BRCA Exchange becomes a single source that collects and organizes existing information. The BRCA Exchange retrieves data from a variety of databases, such as ClinVar and LOVD. Using these sources, BRCA Exchange answers three basic questions about a variant:

* What is the most definitive clinical interpretation available for this variant?
* When was this interpretation made?
* What publicly available data informed this interpretation?

There are two data views used by BRCA Exchange to aid researchers and clinicians who work with _BRCA1_ and _BRCA2._ The default, or Expert Reviewed, view of BRCA Exchange will only display variant aliases, ENIGMA variant interpretations, and a brief interpretation history. The ENIGMA Consortium assesses variant pathogenicity using an expertly developed set of _BRCA1_ and _BRCA2_ classification criteria. The _All Public_ view contains all publicly available data on a variant. Thus, this view of the data provides a more complete, but also more complicated, characterization of a variant. Using multiple large-scale sources and two different portals, BRCA Exchange serves as a contextualized resource for curation of BRCA variant data.


Table of Contents
----------------------------

* [How do I search for a variant?](#how-do-i-search-for-a-variant)
* [Variant Details Page](#variant-details-page)
	* [Variant Nomenclature](#variant-nomenclature)
	* [Clinical Significance Tiles](#clinical-significance-tiles)
	* [Transcript Visualization](#transcript-visualization)
	* [Multifactorial Likelihood Analysis](#multifactorial-likelihood-analysis)
	* [Allele Frequency Reference Sets](#allele-frequency-reference-sets)
		* [ExAC](#exac)
		* [1000 Genomes](#1000-genomes)
		* [Exome Sequencing Project](#exome-sequencing-project)
	* [Cravat/MuPIT Interactive Protein Structure Viewer](#cravatmupit-interactive-protein-structure-viewer)
* [Lollipop Plots](#lollipop-plots)
* [Downloading Variant Data](#downloading-variant-data)
* [How do I use Filters?](#how-do-i-use-filters)
* [How do I use Column Selectors?](#how-do-i-use-column-selectors)

How do I search for a variant?
----------------------------

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



Variant Details Page
============================

## Variant Nomenclature
The Variant Nomenclature fields contain forms of identification by which the variant can be labeled.

* #### Gene Symbol ((Gene_Symbol))
	The Gene Symbol field displays the name of the gene on which the variant was found, as named by HGNC. This will be either _BRCA1_ or _BRCA2_.

* #### Reference cDNA Sequence ((Reference_Sequence))
	This field provides the reference sequence on which the variant was observed, corresponding to the HGVS coordinates indicated in the Nucleotide field e.g. NM\_000492.3, NG\_016465.3, or NC\_000007.13.

* #### HGVS Nucleotide ((HGVS_cDNA))
	This is the HGVS variant alias which references the nucleotide change based on the location in the coding DNA, not the genomic DNA.

* #### HGVS Protein ((HGVS_Protein))
	This alias utilizes predicted protein-level change \(if any\) that would be introduced by the genomic variant.

* #### BIC Variant Identifier ((BIC_Nomenclature))
	This variant alias is presented in BIC Nomenclature, which predates HGVS nomenclature and thus follows a different format.

* #### Genome \(GRCh38\) ((Genomic_Coordinate_hg38))
	This alias utilizes the coordinate of the variant on the GRCh38 reference genome. Click alias to access the UCSC Genome Browser.

* #### Genome \(GRCh37\) ((Genomic_Coordinate_hg37))
	This alias utilizes the coordinate of the variant on the GRCh37 reference genome. Click alias to access the UCSC Genome Browser.

* #### Genome \(GRCh36\) ((Genomic_Coordinate_hg36))
	This alias utilizes the coordinate of the variant on the GRCh36 reference genome. Click alias to access the UCSC Genome Browser.


Clinical Significance Tiles
----------------------------

### ENIGMA ((clinical-significance-enigma))
* #### Clinical Significance ((Clinical_significance_ENIGMA))

    This field displays the variant’s clinical classification and supporting information as provided by ENIGMA expert panel review.

	##### What do these classifications mean?
	* Pathogenic variants confer a markedly increased risk of disease.
	* Likely pathogenic variants have good evidence to support an association with markedly increased disease risk.
	* Likely benign variants have good evidence to support no association with markedly increased disease risk.
	* Benign variants are not associated with any markedly increased risk of disease.
	* Variants of uncertain significance \(VUS\) are those for which the evidence of disease risk is not clear yet, sometimes because there is not yet enough evidence to classify them as either pathogenic or benign.

	Markedly increased disease risk for a variant in _BRCA1_ and _BRCA2_ is currently defined as risk that is sufficient to alter patient management based on detection of that single genetic variant, after considering other non-genetic risk factors for the individual.

* #### Comment on Clinical Significance ((Comment_on_clinical_significance_ENIGMA))
	This field comments on the derivation of the IARC class, provided by ENIGMA.

* #### Assertion Method ((Assertion_method_ENIGMA))
	This field provides the citation or URL describing the method and criteria used to make assertions of clinical significance.

* #### Date last evaluated ((Date_last_evaluated_ENIGMA))
	This is a display of the date on which the clinical significance of the variant was last evaluated by ENIGMA.

* #### Collection Method ((Collection_method_ENIGMA))
	This is the method used to collect the data that supports the assertion of clinical significance. The following values may be listed: case-control, clinical testing, literature only, reference population, research.

* #### Clinical Significance Citation ((Clinical_significance_citations_ENIGMA))
	This field provides citations that document the clinical significance. They may be from PubMed, PubMedCentral, DOI, or NCBI Bookshelf.

* #### Allele Origin ((Allele_origin_ENIGMA))
	Variants are classified as either germline or somatic, depending on how they are acquired. Germline variants are genetic changes that we inherit from our parents. Somatic variants are DNA changes that we acquire over our lifetime, often through exposure to pollutants, toxins, radiation and other carcinogens.


### ClinVar

* #### Clinical Significance ((report-Clinical_Significance_ClinVar))

	These fields, extracted from [ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/), relate to variant clinical classification from individual submitters to ClinVar \(includes repeated information from ENIGMA and BIC tiles\).

	##### What do these classifications mean?
	The below classifications will be found in ClinVar's "clinical significance" field. They are implemented per [recommendations by ACMG/AMP](https://www.ncbi.nlm.nih.gov/pubmed/25741868) for variants interpreted by Mendelian disorders.

	* Pathogenic
	* Likely pathogenic
	* Likely benign
	* Benign
	* Uncertain Significance

* #### Submitter ((report-Submitter_ClinVar))
	This field provides the name of the submitting organization. From ClinVar.
* #### SCV Accession ((report-SCV_ClinVar))
	This is the SCV accession from the ClinVar repository.
* #### Date last updated \(ClinVar\) ((report-Date_Last_Updated_ClinVar))
	The field that shows the date the variant was last submitted to ClinVar. Note that this date is not equivalent to the date last evaluated.
* ##### _Coming Soon:_
	Updated Classification descriptions, Information about tile functionality, resources for star review, definition of summary evidence.

For more information, please visit the [ClinVar website](https://www.ncbi.nlm.nih.gov/clinvar/) or the [ClinVar Data Dictionary](https://www.ncbi.nlm.nih.gov/projects/clinvar/ClinVarDataDictionary.pdf).


### Leiden Open Variation Database (LOVD)

According to LOVD’s homepage, the Leiden Open Variation Database seeks to “provide a flexible, freely available tool for Gene-centered collection and display of DNA variations.” This data is similar to the data in ClinVar, but uses different standards and nomenclatures. Below, you can find a list of the LOVD data fields displayed on BRCA Exchange.

* #### Submitter ((report-Submitters_LOVD))
	The name of the researcher, research group, or individual who submitted the variant data to LOVD.

* #### Variant Effect ((report-Variant_effect_LOVD))
	Variant effect shows the effect of the variant on the function of the protein. It is comparable, but not identical, to clinical significance.

	It uses a symbolic nomenclature:
	* **R/C format**
		* There are two assessments made in LOVD data, which are combined using the R/C format.
		* R is a value reported by the data source/submitter.
		* C is a value concluded by the curator (a LOVD administrator) who reviews data submissions to LOVD.
	* **\+**

	    The variant affects function.
	* **\+?**

	    The variant probably affects function.
	* **\+∗**

	    The variant affects function, but is not associated with individual’s disease phenotype.
	* **\#**

	    The variant affects function, but is not associated with any known disease phenotype.
	* **\-**

	    The variant does not affect function.
	* **\-?**

	    The variant probably does not affect function.
	* **?**

	    The variant effect is unknown.
	* **\.**

	    The variant effect is not classified.

    A symbol is assigned for both the reported classification (from source, “R”) and the concluded classification (from curator, “C”).

    * ##### Examples:
        * \-/\-: denotes a variant that has been determined to not affect function both by a curator and a submitter
        * \+/\+?: denotes a variant that has been determined to affect function by a submitter, but determined as probably affecting function by a curator

* #### Variant Data Type (Origin) ((report-Genetic_origin_LOVD))
	The Data Type field indicates what kind of data was submitted to LOVD.  This field often indicates the source of the variant.

	For the BRCA genes, the following types of data are most common.

	* ##### Germline
		Germline variants are those that the individual is typically born with due to inheritance from one or more parent(s).
	* ##### De novo
		De novo variants are those that an individual is born with, but are not inherited through either parent. These variants might occur via processes such as random mutation.
	* ##### Somatic
		Somatic variants are those that an individual acquires during their lifetime
	* ##### SUMMARY
		Summary records do not describe a variant found in an individual, but summarize multiple records that were provided by the submitter.

	The data reported by LOVD is partially  public and partially private. Because some private data is used, data entries in these fields may look identical from a public view, while perhaps differing in the private information provided to LOVD. Data privacy is often required to secure patient’s personal information.

* #### Individuals ((report-Individuals_LOVD))
	The number of individuals the variant has been observed in, per submission. The individuals field does not account for all LOVD data on a variant, but accounts for the number of individuals observed by the submitter.

* #### Variant ID (Database ID) ((report-DBID_LOVD))
	The ID given to the variant by LOVD. This ID is variant-specific, but not submitter-specific. Please note that this is an important difference from ClinVar’s data specifications.

* #### Created Date ((report-Created_date_LOVD))
	The date that variant data was first submitted to LOVD.

* #### Edited Date ((report-Edited_date_LOVD))
	The date that the variant data was last updated or changed in LOVD to reflect new evidence.


### BIC ((clinical-significance-bic))
* #### Patient Nationality ((Patient_nationality_BIC))
	The field that indicates nationality of the patient associated with the variant.
* #### Ethnicity ((Ethnicity_BIC))
	The field that provides the ethnicities of the observed patients.
* #### Family members carrying this variant ((Number_of_family_member_carrying_mutation_BIC))
	Number of family members carrying this variant.
* #### Literature Reference ((Literature_citation_BIC))
	Literature Reference\(s\). From BIC.
* #### Allele Origin ((Germline_or_Somatic_BIC))
	Origin of the allele, Germline \(G\) or Somatic \(S\).

	<p>Variants are classified as either germline or somatic, depending on how they are acquired. Germline variants are genetic changes that we inherit from our parents. Somatic variants are DNA changes that we acquire over our lifetime, often through exposure to pollutants, toxins, radiation and other carcinogens.</p>



Transcript Visualization
----------------------------

The transcript visualization depicts the reference transcript of the gene \(either _BRCA1_ or _BRCA2_\) on which the selected variant is found. The gene is drawn as a series of exons with intervening introns. Note that exons are drawn scaled \(nonlinearly\) to their actual size, to make the smaller exons easier to see, while introns are drawn at a fixed width.

The gene is shown in full at the top of the visualization, with the exons/introns in which the variant lies \(i.e., the area of interest\) highlighted in yellow. The area of interest is shown zoomed-in below the full depiction, but contains the same information.

The variant's modifications to the reference transcript are depicted in three colored regions:
* Substitution variants \(depicted in lime green\)
* Deletions \(depicted in red with cross-hatched black lines\)
* Insertions \(depicted in blue\)

The visualization also displays other regions of interest for the gene, such as donor and acceptor sites and clinically important functional domains. The visibility of these regions can be toggled individually via the checkboxes below the transcript.

Clinically important functional domains are regions of the gene that ENIGMA has assessed to be relevant to risk. More specifically, these domains are regions of the gene that ENIGMA has determined to be relevant to risk when modified by missense alterations or in-frame deletions that could encode a stable \(non-truncated\) protein. This determination is made using existing knowledge of pathogenicity for individual variants in these domains. Clicking the label for a given clinically important domain will display the corresponding domain in an animated outline. Clicking again on the label will disable the outline, as will clicking another domain designated as clinically important.


## [Multifactorial Likelihood Analysis Classification](https://www.ncbi.nlm.nih.gov/pubmed/21990134) ((multifactorial-likelihood-analysis))
Component likelihoods and final classification are from multifactorial likelihood analysis, which is a quantitative integrated evaluation of variant pathogenicity.
* #### Posterior Probability of pathogenicity \(ExUV\) ((Posterior_probability_exLOVD))
	This value is the posterior probability of pathogenicity from ExUV.
* #### Prior probability of pathogenicity \(ExUV\) ((Combined_prior_probablility_exLOVD))
	This value is the combined prior probability in favor of pathogenicity. It is a combination of the missense analysis prior probability and the splicing analysis prior probability. Generally, it is the higher of these two prior probabilities. From ExUV.
* #### Missense analysis probability of pathogenicity \(ExUV\) ((Missense_analysis_prior_probability_exLOVD))
	This prior probability estimate combines position in the protein with an evaluation of missense substitutions that fall in the proteins key functional domains. From ExUV.
* #### Co-occurrence likelihood \(ExUV\) ((Co_occurrence_LR_exLOVD))
	This value is the likelihood ratio based on the frequency of co-occurrence between the variant of interest and clearly pathogenic variants in the same gene. From ExUV.
* #### Segregation Likelihood Ratio \(ExUV\) ((Segregation_LR_exLOVD))
	This value is the likelihood ratio based on segregation analysis
* #### Summary Family History Likelihood Ratio \(ExUV\) ((Sum_family_LR_exLOVD))
	This value is the likelihood ratio based on an analysis of the severity of summary family histories of breast and/ or ovarian cancer.


Allele Frequency Reference Sets
----------------------------

The allele frequency reference sets show the frequency of a BRCA1 or BRCA2 variant in a reference population. To view or collapse all nested tiles, click the arrows available at the top right of this tile. Though the two sets of populations in ExAC and 1000 Genomes closely resemble each other, they are not identical. ESP also uses different population categories.

### ExAC

ExAC is a data source that provides BRCA1 and BRCA2 allele frequencies for the BRCA Exchange. The goal of ExAC is to “aggregate and harmonize exome sequencing data from a variety of large-scale sequencing projects, and to make summary data available for the wider scientific community” \([About ExAC](http://exac.broadinstitute.org/about)\). The ExAC data set used by BRCA exchange excludes data from [TCGA](https://tcga-data.nci.nih.gov/docs/publications/tcga/about.html), to ensure that frequencies used to assess pathogenicity are not skewed by sampling errors. ExAC data will soon be updated to the GnomAD data set.

For more information about ExAC, please refer to the [ExAC browser](http://exac.broadinstitute.org/) or their [flagship publication](https://www.nature.com/articles/nature19057).


#### Graphical ExAC Data
Graphical ExAC data can be viewed by expanding the ExAC \(Graphical\) nested tile.  Two Graphs are available; one of the graphs is custom scaled to the allele frequencies by default \(right side\). Hovering over each bar will give you the numerical value represented in the population subset. You can click anywhere on the ExAC \(scaled\) graph to change the scale between 1.0% \(.01\), 0.1% \(.001\), and the custom, default scale. Because some Allele Frequencies can be very small, a variety of scales will allow you to view all possible Allele Frequencies graphically.

Each group found on the x-axis of the bar chart can be found in the list of fields described in the ExAC \(Numerical\) section.


#### Numerical ExAC Data
Numerical ExAC data fields show numerical minor allele frequency data associated with each population, as well as the overall allele frequency. All of the minor allele frequencies are consistent with the graphs shown in the ExAC \(Graphical\) nested tile.
* #### Allele Frequency \(ExAC minus TCGA\) ((Allele_frequency_ExAC))
	Minor allele frequency, per ExAC \(excluding TCGA data\)
* #### African/African American \(AFR\) ((Allele_frequency_AFR_ExAC))
	Allele frequency in African/African American populations, per ExAC
* #### Admixed American/Latino \(AMR\) ((Allele_frequency_AMR_ExAC))
	Allele frequency in Admixed American/Latino populations, per ExAC
* #### East Asian \(EAS\) ((Allele_frequency_EAS_ExAC))
	Allele frequency in East Asian populations, per ExAC
* #### Finnish \(FIN\) ((Allele_frequency_FIN_ExAC))
	Allele frequency in Finnish populations, per ExAC and separated from European because of an enriched data set
* #### Non-Finnish European \(NFE\) ((Allele_frequency_NFE_ExAC))
	Allele frequency in Non-Finnish European populations, per ExAC
* #### South Asian \(SAS\) ((Allele_frequency_SAS_ExAC))
	Allele frequency in South Asian populations, per ExAC
* #### Other \(OTH\) ((Allele_frequency_OTH_ExAC))
	Allele frequency in populations other than those listed above, per ExAC


### 1000 Genomes

The 1000 Genomes Project also contributes allele frequency data to BRCA Exchange. The goal of the 1000 Genomes Project is “to provide a comprehensive description of common human genetic variation by applying whole-genome sequencing to a diverse set of individuals from multiple populations” \([The 1000 Genomes Project Consortium](https://www.nature.com/articles/nature15393)\).  Whereas ExAC aggregates sequence data, 1000 Genomes works on adding new, more diverse genomes to existing data sets.

For more information about 1000 Genomes, please visit the [International Genome Sample Resource](http://www.internationalgenome.org/) website or their [publications](http://www.internationalgenome.org/1000-genomes-project-publications).

#### Graphical 1000 Genomes Data
Graphical 1000 Genomes data can be viewed by expanding the 1000 Genomes \(Numerical\) nested tile. Two Graphs are available; one of the graphs is custom scaled to the allele frequencies by default \(right side\). Hovering over each bar will give you the numerical value represented in the population subset. You can click anywhere on the 1000 Genomes \(scaled\) graph to change the scale between 1.0% (.01), 0.1% (.001), and the custom, default scale. Because some Allele Frequencies can be very small, a variety of scales will allow you to view all possible Allele Frequencies graphically.

Each group found on the x-axis of the bar chart can be found in the list of fields described in the 1000 Genomes \(Numerical\) section.

#### Numerical 1000 Genomes Data
Numerical 1000 Genomes data fields show numerical minor allele frequency data associated with each population, as well as the overall allele frequency. All of the minor allele frequencies are consistent with the graphs shown in the 1000 Genomes \(Graphical\) nested tile.
* #### Allele Frequency ((Allele_frequency_1000_Genomes))
	Overall allele frequency, per 1000 Genomes
* #### AFR Allele Frequency ((AFR_Allele_frequency_1000_Genomes))
	Allele frequency in African-American populations, per 1000 Genomes
* #### AMR Allele Frequency ((AMR_Allele_frequency_1000_Genomes))
	Allele frequency in Admixed American populations, per 1000 Genomes
* #### EAS Allele Frequency ((EAS_Allele_frequency_1000_Genomes))
	Allele frequency in East Asian populations, per 1000 Genomes
* #### EUR Allele Frequency ((EUR_Allele_frequency_1000_Genomes))
	Allele frequency in European populations, per 1000 Genomes
* #### SAS Allele Frequency ((SAS_Allele_frequency_1000_Genomes))
	Allele frequency in Southeast Asian Populations, per 1000 Genomes


### Exome Sequencing Project

#### Numerical ESP Data
The NHLBI GO Exome Sequencing Project is yet another contributor of allele frequency data to the BRCA Exchange. The Exome Sequencing Project \(ESP\) aims to “to discover novel genes and mechanisms contributing to heart, lung and blood disorders by pioneering the application of next-generation sequencing of the protein coding regions of the human genome across diverse, richly-phenotyped populations and to share these datasets and findings with the scientific community to extend and enrich the diagnosis, management and treatment of heart, lung and blood disorders” \([NHLBI Exome Sequencing Project](http://evs.gs.washington.edu/EVS/)\).

* #### Allele Frequency (ESP) ((Allele_Frequency_ESP))
	Allele frequency in entire data set, per ESP
* #### EA Allele Frequency (ESP) ((EA_Allele_Frequency_ESP))
	Allele frequency in European-American populations, per ESP
* #### AA Allele Frequency (ESP) ((AA_Allele_Frequency_ESP))
	Allele frequency in African-American populations, per ESP


CRAVAT/MuPIT Interactive Protein Structure Viewer ((cravat-mupit-3d-protein-view))
----------------------------

For missense variants that occur within a region of the protein with a well-defined three-dimensional structure, bioinformatics and protein structure analysis can help suggest the impact of the variation (1-3). Certain regions in protein structures tend to be more sensitive to variation, such as positions buried within the core of the protein (where variation could destabilize the protein’s structure) and positions near binding sites (where variation could impact the protein’s function). The MuPIT interactive viewer from the CRAVAT project (4) facilitates such analysis by showing the position of the variant in the context of its three-dimensional protein structure.

By clicking on the CRAVAT/MuPIT thumbnail image, the user opens a new browser tab running the interactive MuPIT viewer.  The display is divided into three panels.

The center panel displays the three-dimensional protein structure, and indicates the variant’s location with spacefill spheres. The positions in the structure are also color-coded, per predictions about the impact of variation in those regions. By default, colors at a position are assigned according to the most severe in silico prediction of pathogenicity (6). For example, if the reference base was a ‘C’, and the possible missense variants had respective in silico pathogenicity predictions of 0.67 for ‘A’, 0.33 for ‘G’, and 0.88 for ‘T’, then the position would be color-coded according to the highest prediction score, which is 0.88 for ‘T’. The user can rotate the protein structure by holding down and dragging the mouse.

The rightmost panel provides controls to allow the user to select a different color map, and to change the color of the variant. The available color maps describe:
* In silico prediction of the most sever variant at each position (6)
* Multifactorial analysis, incorporating in silico prediction with probabilities estimated from case-level patient data (7)

Coming soon: color maps currently under development will describe known pathogenic and benign variants, augmented with in silico prediction of whether the pathogenic variants are pathogenic due to impact on splicing or impact on protein function (6).

The leftmost panel displays additional controls, including a link for further help.

CRAVAT and MuPIT are available for missense variants that map to positions within curated three-dimensional protein structures.  It is not available for insertions, deletions, duplications, or positions in regions of the protein where no three-dimensional protein structure has been determined.

##### References

1. [Carvalho et al. 2009. PMID 18992264](https://www.ncbi.nlm.nih.gov/pubmed/18992264)
2. [Karchin et al. 2008. PMID 19043619](https://www.ncbi.nlm.nih.gov/pubmed/?term=19043619)
3. [Karchin et al. 2007. PMID 17305420.](https://www.ncbi.nlm.nih.gov/pubmed/?term=17305420)
4. [Masica et al. 2017. PMID 29092935.](https://www.ncbi.nlm.nih.gov/pubmed/29092935)
5. [Vall&#x00e9;e et al. 2016. PMID 26913838.](https://www.ncbi.nlm.nih.gov/pubmed/26913838)
6. [Vall&#x00e9;e et al. 2012. PMID 21990165.](https://www.ncbi.nlm.nih.gov/pubmed/21990165)


Lollipop Plots
============================

Lollipop plots are a tool to visualize the chromosomal position and pathogenicity classification for each variant in a gene. Here, each circle-capped ‘lollipop’ indicates whether a BRCA1 and BRCA2 variant is pathogenic (labeled in red), benign (labeled in light-blue), or the variants clinical significance is uncertain (labeled in dark-blue). The y-axis represents the pathogenicity classification of a variant. The x-axis represents the genomic coordinates of those variants in GRCh38 human reference space. The bottom box that runs along the x-axis of the diagram displays the position of each exon in the selected gene.

To open the lollipop chart click on the ‘Show Lollipop Chart’ button located at the top of the Variants page. To alternate between the BRCA1 and BRCA2 lollipop charts, click on the tab that contains the relevant gene name located in the upper-left hand corner next to the lollipop chart. To select a specific area of the gene to zoom in on a region, simply click and drag on the lower box to create a shaded box which displays which region you have selected. You can also select for specific exon ranges by clicking on the relevant colored exon box on the lower box-chart. Mouse scrolling on the main chart itself will also zoom in and out of a region. To de-select or zoom all the way back out, you can either click on an unshaded space in the lower box or click the ‘Hide Lollipop Chart’ button at the top twice to refresh and reset the chart. On the legend located within the chart you can see the relative percentages of pathogenicity classifications of variants that match the search and filter criteria for a particular gene. To hide/dim a class of variants set of lollipops, click on the colored box on the legend. You can reset this by clicking on the colored box again.

<div style="max-width: 1000px; margin-left: auto; margin-right: auto; height: 0px; position: relative; padding-bottom: 59.552%;"><iframe src="https://streamable.com/e/zx9c?muted=1&amp;autoplay=1&amp;hd=1" frameborder="0" allowfullscreen webkitallowfullscreen mozallowfullscreen scrolling="no" style="width: 100%; height: 100%; position: absolute;"></iframe></div>


Downloading Variant Data
============================

To download the variant data click on the ‘Download’ button located above the search bar within the Variants Table page. When you click on it the data that is downloaded is a comma delimited .csv file containing the set of variant details for the variants that matched the search and/or filtering criteria. The first row in the file represent the column labels.


How do I use Column Selectors?
============================

Columns allow you to filter the variant table (in search results) according to fields found in the tiles on the Variant Details Page. Column selectors correspond to each tile’s fields, whose descriptions are found above. These selectors allow you to include any tile’s field as its own column in your search results. You can then sort results using added columns.

##### _Updates Coming Soon_


How do I use Filters?
============================

##### _Updates Coming Soon_
