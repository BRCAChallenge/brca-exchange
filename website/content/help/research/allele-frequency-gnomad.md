
gnomAD provides _BRCA1_ and _BRCA2_ allele frequencies for the BRCA Exchange. The goal of gnomAD is to “aggregate and harmonize exome and genome sequencing data from a variety of large-scale sequencing projects, and to make summary data available for the wider scientific community.” Its predecessor, the Exome Aggregation Consortium (ExAC) contained only exome data. The gnomAD data set used by BRCA Exchange is the “non-cancer” subset, and excludes data from TCGA and other cancer cohorts. This ensures that frequencies used to represent normal, unaffected populations are not skewed by the inclusion of cancer patients. 

As an aggregation database, gnomAD data subsumes much of ExAC and 1000 Genomes data. 

#### Exomes and Genomes 
gnomAD provides frequency from both exome and genome sequencing projects. On the BRCA Exchange, these frequencies are presented  in different subtiles. A small "E" indicates Exomes data, while a small "G" icon indicates Genomes data. If data is available on both Exomes and Genomes, a total of four subtiles will be available:

* gnomAD Exomes (Graphical) 
* gnomAD Exomes (Numerical) 
* gnomAD Genomes (Graphical) 
* gnomAD Genomes (Numerical) 

This data is updated during the monthly release cycle, and will contain any new information provided by gnomAD. 

#### gnomAD Populations

The BRCA Exchange displays the primary subsets used by ExAC (i.e. “East  Asian”, but does not display detailed population data (i.e. “Japanese,” “Korean”). Detailed populations, as well as female and male distributions, can be found in gnomAD. 

* ##### Allele Frequency (gnomAD non-cancer cohort)
	* Minor allele frequency, per ExAC (excluding cancer sequencing data)
* ##### African/African American (AFR)
	* Allele frequency in African/African American populations, per gnomAD
* ##### Latino (AMR)
	* Allele frequency in Admixed American/Latino populations, per gnomAD
* ##### East Asian (EAS)
	* Allele frequency in East Asian populations, per gnomAD
	* Subsets available on gnomAD:
		* Japanese 
		* Korean 
		* Other East Asian
* ##### Finnish (FIN)
	* Allele frequency in Finnish populations, per gnomAD and separated from European because of an enriched data set
* ##### Non-Finnish European (NFE)
	* Allele frequency in Non-Finnish European populations, per gnomAD
	* Subsets available on gnomAD:
		* Bulgarian
		* Estonian 
		* North-western European
		* Other non-Finnish European
		* Southern European
		* Swedish
* ##### South Asian (SAS)
	* Allele frequency in South Asian populations, per gnomAD
* ##### Other (OTH)
	* Allele frequency in populations that “did not unambiguously cluster with the major populations in a principal component analysis,” per gnomAD
	* South Asian genomes only contain 31 samples; thus for genome data, the SAS population is grouped with Other.

#### Graphical gnomAD Data
Graphical gnomAD data can be viewed by expanding the gnomAD (Graphical) subtile. Two Graphs are available; one of the graphs is custom scaled to the allele frequencies by default (right side). Hovering over each bar will give you the numerical value represented in the population subset. You can click anywhere on the gnomAD (scaled) graph to change the scale between 1.0% (.01), 0.1% (.001), and the custom, default scale. Because some Allele Frequencies can be very small, a variety of scales will allow you to view all possible Allele Frequencies graphically.
Each group found on the x-axis of the bar chart can be found in the list of fields described in the gnomAD Populations section.

#### Numerical gnomAD Data
Numerical gnomAD data fields show numerical minor allele frequency data associated with each population, as well as the overall allele frequency. All of the minor allele frequencies are consistent with the graphs shown in the gnomAD (Graphical) subtile.

#### Exomic and Genomic Allele Frequencies for Variant Interpretation
There is not yet scientific consensus about whether exomic, genomic, or a combined allele frequency is most reliable for variant interpretation. The data display currently provides them separately. If you have expertise or opinions about how to use exomic and genomic allele frequencies together, please feel free to [contact us](mailto:brca-exchange-contact@genomicsandhealth.org?subject=BRCA%20Exchange%20Literature%20Search). 

#### gnomAD Flags
A flag on a variant can serve one of two purposes: either it indicates that there is something unusual about the data (such as low quality), or it indicates an annotation for the variant. For further information about why a variant might be flagged, use the gnomAD link in the tile to view the variant in gnomAD, along with detailed quality metrics. For more information about flags, visit the MacArthur Lab’s [blog post](https://macarthurlab.org/2018/10/17/gnomad-v2-1/) about the latest gnomAD release. 

#### Transcripts

The functional impact of a variant depends in part on where it appears in the RNA transcript.  The choice of transcript can affect the computation of the flags in gnomAD.  While gnomAD computes information relative to many transcripts, by default it shows information relative to the ENSEMBL canonical transcripts.  For _BRCA1_, this transcript is not clinically-relevant, and BRCA Exchange shows data from a different transcript.  For _BRCA2_, it shows data relative to gnomAD’s default transcript.  

The BRCA Exchange displays data associated with ENSEMBL transcript ENST00000357654 for BRCA1, and ENSMBL transcript ENST00000544455 for BRCA 2. These correspond to RefSeq transcript [NM_007294.3](http://ftp.ebi.ac.uk/pub/databases/lrgex/LRG_292.xml) and [NM_000059.3](http://ftp.ebi.ac.uk/pub/databases/lrgex/LRG_293.xml) respectively (per LRG).  Additional data, including detailed populations, quality scores, and flags relative to other transcripts, are available at gnomAD.

For more information about efforts to unify NCBI and EMBL-EBI datasets, please visit [MANE](https://ncbiinsights.ncbi.nlm.nih.gov/2019/03/12/mane-select-v0-5/) (Matched Annotation from the NCBI and EMBL-EBI). Though they have not yet covered _BRCA1_ and _BRCA2_, the MANE Select set will eventually comprise of “a matched representative transcript for every human protein-coding gene.”

