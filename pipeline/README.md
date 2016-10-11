# brca-pipeine

Contains source code used to wrangle the BRCA variation data that is presented in the brca-exchange site.
##### Overall Dependencies
 * create a directory for genomic resources.  Point to it with the environment variable $BCRA_RESOURCES
 * download http://hgwdev.soe.ucsc.edu/~cline/BRCA/resources/refseq_annotation.hg38.gp to $BRCA_RESOURCES
 * download http://hgwdev.soe.ucsc.edu/~cline/BRCA/resources/refseq_annotation.hg19.gp to $BRCA_RESOURCES
 * download http://hgwdev.soe.ucsc.edu/~cline/BRCA/resources/refseq_annotation.hg18.gp to $BRCA_RESOURCES
 * download http://hgwdev.soe.ucsc.edu/~cline/BRCA/resources/hg19.fa.gz to $BRCA_RESOURCES and uncompress
 * download http://hgwdev.soe.ucsc.edu/~cline/BRCA/resources/hg18.fa.gz to $BRCA_RESOURCES and uncompress
 * download http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz to $BRCA_RESOURCES and uncompress
 * download http://hgdownload.cse.ucsc.edu/gbdb/hg19/liftOver/hg19ToHg38.over.chain.gz to $BRCA_RESOURCES.  There is no need to uncompress it.
 * Download brca1_hg38.txt, brca2_hg38.txt, brca1_hg19.txt and brca2_hg19.txt from http://hgwdev.soe.ucsc.edu/~cline/BRCA/resources/ to $BRCA_RESOURCES
 * download the latest available version of BRCA2.txt from the lovd subdirectory and move it into $BRCA_RESOURCES (can be found on the production VM)
 * install postgresql
 * install the leiden cleanup package at http://leiden-open-variation-database-lovd-cleanup.readthedocs.io/en/latest/index.html (NOTE: this can be problematic, please submit an issue if you encounter errors with the leiden package)
 * pip install git+https://github.com/counsyl/hgvs.git
 * pip install psycopg2
 * pip install parsley
 * pip install configparser (Enigma requires version 3.5.0, will fail with latest version)
 * pip install bioutils
 * pip install biopython
 * install lzo
 * install bx-python from https://pypi.python.org/packages/55/db/fa76af59a03c88ad80494fc0df2948740bbd58cd3b3ed5c31319624687cc/bx-python-0.7.3.tar.gz#md5=d8c50c01c9e421bae0bbdbfa00fef6e4 (note: pip install bx-python left an outstanding dependency to lzo1x.h, which proved to be hard to resolve) NOTE: This may lead to errors in the bx.bbi directory, please submit an issue if you encounter this problem. Sometimes `pip install bx-python` is also necessary.
 * pip install pysam
 * install lzo (sudo apt-get install liblzo2-dev)
 * pip install git+https://github.com/bxlab/bx-python
 * install Cython
 * pip install CrossMap (8/25/16: install from
 * install vcftools from https://vcftools.github.io/index.html
 * pip install pyvcf
 * install tabix (see http://genometoolbox.blogspot.com/2013/11/installing-tabix-on-unix.html)
 * pip install luigi
 * pip install retrying
 * pip install synapseclient[pandas,pysftp]
 * cd into enigma and `pip install -r requirements.txt`

## Misc Instructions
### Convert refseq .psl file to .gp (genepred) format (required format for hgvs conversion)
  This is the gene feature coordinate file that coincides with the '-r' option of umd2vcf and bic2vcf scripts.
  
  1. Add '/cluster/bin/x86_64/mrnaToGene' to your PATH environment variable
  2. mrnaToGene [options] psl genePredFile
  3. Insert an extra column on the left-hand most side for each row in the genepred file and put any number there designating the id of the refseq annotation. The exact number doesn't matter as long as it is unique. This is needed for proper formatting so that the hgvs python package can properly interpret the genepred file.
  4. Add the open reading frame coordinates in the genepred file (column 7 is start-codon position and column 8 is the position at the end of the stop-codon)

  e.g. mrnaToGene -insertMergeSize=-1 -noCds refseq_annotation.hg19.psl refseq_annotation.hg19.gp

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

## Luigi

Luigi is intended to automate many of the tasks in the pipeline directory. To get started, `pip install luigi`.
