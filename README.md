# brca
A demonstration of the GA4GH api system whose goals are to showcase the capabilities of the system in distributing and categorizing large quantities of BRCA variant data in a manner that is consistant with the needs of the genomics and biomedical community. So far there is support and implementation of a reference server for BRCA1 and BRCA2 variants contained within 1000-genomes, ClinVAR, LOVD, ex-lovd, ExAC, UMD, and BIC databases.

##Authors
[Benedict Paten](https://github.com/benedictpaten/), [Charles Markello](https://github.com/cmarkello), [Molly Zhang](https://github.com/MollyZhang), [Max Haeussler](https://github.com/maximilianh), [Melissa Cline](https://github.com/melissacline), [Mark Diekhans](https://github.com/diekhans)

## Repository directory
  bic: Contains relevant scripts pertaining to converting data originating from the [Breast Cancer Information Core](https://research.nhgri.nih.gov/projects/bic/index.shtml) database to GA4GH reference server friendly vcf format.
  
  umd: Contains relevant scripts and web-scrapped data pertaining to converting data originating from the [Universal Mutation Database](http://www.umd.be/BRCA1/) to GA4GH reference server friendly vcf format.
  
  doc: Directory contains miscellaneous logfiles of various runs of scripts that are in development.

## Misc Instructions
### Convert refseq .psl file to .gp (genepred) format (required format for hgvs conversion)
  1. Add '/cluster/bin/x86_64/mrnaToGene' to your PATH environment variable
  2. mrnaToGene [options] psl genePredFile
  3. Add the open reading frame coordinates in the genepred file (column 6 is start position and column 7 is the stop position)

  e.g. mrnaToGene -insertMergeSize=-1 -noCds refseq_annotation.hg19.psl refseq_annotation.hg19.gp
