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
