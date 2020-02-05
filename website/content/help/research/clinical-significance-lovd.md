According to LOVD’s homepage, the Leiden Open Variation Database seeks to “provide a flexible, freely available tool for Gene-centered collection and display of DNA variations.” LOVD is a part of the standardized technology infrastructure of the [Human Variome Project](https://www.humanvariomeproject.org/), and utilizes a large network of submitters to collect variant data. LOVD is available both as a downloadable, privately managed software package and as a shared online database. Accordingly, this tool is designed to meet the diverse and interdisciplinary needs of the variant curation community. Any consortia can utilize and privately manage their own data. Additionally, LOVD also aggregates shared data on a gene-by-gene basis, and often utilizes the knowledge of an expert curator to manage variants for a specific gene. LOVD is also compliant with HGVS guidelines, and collaborates with other tools to ensure data quality and compatibility. For more information about LOVD, read their [documentation](https://databases.lovd.nl/shared/docs/), visit their [home page](https://www.lovd.nl/3.0/home), or read the [Wikipedia entry](https://en.wikipedia.org/wiki/Leiden_Open_Variation_Database) on their work. 

BRCA Exchange aggregates data from the shared LOVD database for _BRCA1_ and _BRCA2_. LOVD data is similar to the data in ClinVar, in that it contains clinical interpretations, but uses somewhat different standards and nomenclatures. Below, you can find a list of the LOVD data fields displayed on BRCA Exchange. 

* #### Submitter ((report-Submitters_LOVD))
	The name of the researcher, research group, or individual who submitted the variant data to LOVD.

* #### Clinical Classification ((report-Classification_LOVD))
	The classification of variant based on clinical consequences, generally using standardised criteria.

	Please note that there may be variability in this field, since different submitters might use different systems of classification. If you have questions or concerns about a classification, visit LOVD for more information.  

* #### Variant Data Type (Origin) ((report-Genetic_origin_LOVD))
	The Data Type field indicates what kind of data was submitted to LOVD. This field often indicates the source of the variant. For the BRCA genes, the following types of data are most common.
	* ##### Germline
		Germline variants are those that the individual is born with, due to inheritance from one or more parent(s).
	* ##### De novo
		De novo variants are those that an individual is born with, but are not inherited through either parent. These variants might occur via processes such as random mutation.
	* ##### Somatic
		Somatic variants are those that an individual acquires during their lifetime
	* ##### SUMMARY
		Summary records do not describe a variant found in an individual, but summarize multiple records that were provided by the submitter.

* #### Frequency ((report-Variant_frequency_LOVD)) 
	The frequency with which the variant was observed in the cohort studied.

* #### Individuals ((report-Individuals_LOVD))
	The number of individuals the variant has been observed in, per submission. The individuals field does not account for all LOVD data on a variant, but accounts for the number of individuals observed by the submitter.

* #### Variant ID (Database ID) ((report-DBID_LOVD))
	The ID given to the variant by LOVD. This ID is variant-specific, but not submitter-specific. Please note that this is an important difference from ClinVar’s data specifications.

* #### Variant Haplotype ((report-Variant_haplotype_LOVD))
	The haplotype on which the variant was found.

* #### Created Date ((report-Created_date_LOVD))
	The date that variant data was first submitted to LOVD.

* #### Edited Date ((report-Edited_date_LOVD))
	The date that the variant data was last updated or changed in LOVD to reflect new evidence.

* #### Variant Remarks ((report-Remarks_LOVD))
	Remarks regarding the variant described, provided by submitter. 


