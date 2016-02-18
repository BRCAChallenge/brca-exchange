# brca
A demonstration of the GA4GH api system whose goals are to showcase the capabilities of the system in distributing and categorizing large quantities of BRCA variant data in a manner that is consistant with the needs of the genomics and biomedical community. So far there is support and implementation of a reference server for BRCA1 and BRCA2 variants contained within 1000-genomes, ClinVAR, LOVD, ex-lovd, ExAC, UMD, and BIC databases.

##Authors
[Benedict Paten](https://github.com/benedictpaten/), [Charles Markello](https://github.com/cmarkello), [Molly Zhang](https://github.com/MollyZhang), [Max Haeussler](https://github.com/maximilianh), [Melissa Cline](https://github.com/melissacline), [Mark Diekhans](https://github.com/diekhans)

## Repository directory
  bic: Contains relevant scripts pertaining to converting data originating from the [Breast Cancer Information Core](https://research.nhgri.nih.gov/projects/bic/index.shtml) database to GA4GH reference server friendly vcf format.
  
  umd: Contains relevant scripts and web-scrapped data pertaining to converting data originating from the [Universal Mutation Database](http://www.umd.be/BRCA1/) to GA4GH reference server friendly vcf format.
  
  doc: Directory contains miscellaneous logfiles of various runs of scripts that are in development.

## Misc Instructions
### Restart GA4GH reference server
  1. Log in to the AWS instance at `ubuntu@ec2-54-148-207-224.us-west-2.compute.amazonaws.com` using ssh
  2. type "screen -r"
  3. ctrl + c to stop the server
  4. restart the server with command 'python server_dev.py -f /srv/ga4gh/brca_config.py'
  5. ctrl + AD to detach from screen
  6. logout of the AWS instance

### Convert refseq .psl file to .gp (genepred) format (required format for hgvs conversion)
  This is the gene feature coordinate file that coincides with the '-r' option of umd2vcf and bic2vcf scripts.
  
  1. Add '/cluster/bin/x86_64/mrnaToGene' to your PATH environment variable
  2. mrnaToGene [options] psl genePredFile
  3. Insert an extra column on the left-hand most side for each row in the genepred file and put any number there designating the id of the refseq annotation. The exact number doesn't matter as long as it is unique. This is needed for proper formatting so that the hgvs python package can properly interpret the genepred file.
  4. Add the open reading frame coordinates in the genepred file (column 7 is start-codon position and column 8 is the position at the end of the stop-codon)

  e.g. mrnaToGene -insertMergeSize=-1 -noCds refseq_annotation.hg19.psl refseq_annotation.hg19.gp

### Data extraction and conversion to vcf:
  Requires installation of hgvs python package as can be found here: https://github.com/counsyl/hgvs
  
#### Generate ClinVar VCF files (`ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/xml/`)
See `https://github.com/BD2KGenomics/brca/blob/master/clinvar/README.txt`

#### Generate umd vcf files from webscrapped data and upload to server
See `https://raw.githubusercontent.com/BD2KGenomics/brca/master/umd/README.txt`

#### Generate bic vcf files from webscrapped data and upload to server
See `https://raw.githubusercontent.com/BD2KGenomics/brca/master/bic/README.txt`
 
#### Webscrap data from LOVD, generate LOVD vcf files and upload to server
  1. Install leiden package and its dependencies via command `python setup.py install` in the `leidenv1.0_package` directory.

##### Webscrap exLOVD (`http://hci-exlovd.hci.utah.edu/`) 
See `https://raw.githubusercontent.com/BD2KGenomics/brca/master/lovd/README.md`

##### Webscrap sharedLOVD (`http://databases.lovd.nl/shared/`)
See `https://raw.githubusercontent.com/BD2KGenomics/brca/master/lovd/README.md`

#### ExAC

#### ESP (`http://evs.gs.washington.edu/EVS/`)
See `https://github.com/BD2KGenomics/brca/edit/master/esp/README.txt`


