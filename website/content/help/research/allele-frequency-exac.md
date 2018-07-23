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
