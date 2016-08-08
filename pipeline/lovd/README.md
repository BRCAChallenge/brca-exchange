
##### Webscrap exLOVD (`http://hci-exlovd.hci.utah.edu/`)
  1. create an output directory.  Point to it with the environment variable $EXLOVD
  2. Extract variant data to .txt flat file format: e.g. `extract_data.py -u http://hci-exlovd.hci.utah.edu/ -l BRCA1 BRCA2 -o $EXLOVD`
  3. Convert extracted flat file to vcf format `./lovd2vcf -i $EXLOVD/BRCA1.txt -o $EXLOVD/exLOVD_brca1.hg19.vcf -a exLOVDAnnotation -b 1  -r $BRCA_RESOURCES/refseq_annotation.hg19.gp -g $BRCA_RESOURCES/hg19.fa`
  3. `./lovd2vcf -i output_directory/BRCA2.txt -o exLOVD_brca2.vcf -a $EXLOVD/exLOVDAnnotation -b 2 -r $BRCA_RESOURCES/refseq_annotation.hg19.gp -g $BRCA_RESOURCES/hg19.fa`                     
  4. `vcf-concat $EXLOVD/exLOVD_brca1.hg19.vcf $EXLOVD/exLOVD_brca2.hg19.vcf > $EXLOVD/exLOVD_brca12.hg19.vcf`
  5. `CrossMap.py vcf $BRCA_RESOURCES/hg19ToHg38.over.chain.gz $EXLOVD/exLOVD_brca12.hg19.vcf $BRCA_RESOURCES/hg38.fa $EXLOVD/exLOVD_brca12.hg38.vcf
  6. `vcf-sort $EXLOVD/exLOVD_brca12.hg38.vcf > $EXLOVD/exLOVD_brca12.sorted.hg38.vcf`
  7. `cp $EXLOVD/exLOVD_brca12.sorted.hg38.vcf $PIPELINE_INPUT

##### Webscrap sharedLOVD (`http://databases.lovd.nl/shared/`)
  1. create an output directory.  Point to it with the environment variable $LOVD
  2. `extract_data.py -u http://databases.lovd.nl/shared/ -l BRCA1 BRCA2 -o $LOVD`
  3. `./lovd2vcf -i $LOVD/BRCA1.txt -o $LOVD/sharedLOVD_brca1.hg19.vcf -a sharedLOVDAnnotation -b 1 -r $BRCA_RESOURCES/refseq_annotation.hg19.gp -g $BRCA_RESOURCES/hg19.fa`
  4. `./lovd2vcf -i $LOVD/BRCA2.txt -o $LOVD/sharedLOVD_brca2.hg19.vcf -a sharedLOVDAnnotation -b 2 -r $BRCA_RESOURCES/refseq_annotation.hg19.gp -g $BRCA_RESOURCES/hg19.fa`
  # Comment: BRCA2.txt might or might not be created.  It seems like there's an error condition that kills BRCA2.  It's not worth fixing the error condition at this point.
  5. `vcf-concat $LOVD/sharedLOVD_brca1.hg19.vcf $LOVD/sharedLOVD_brca2.hg19.vcf > $LOVD/sharedLOVD_brca12.hg19.vcf`
  6. `CrossMap.py vcf $BRCA_RESOURCES/hg19ToHg38.over.chain.gz $LOVD/sharedLOVD_brca12.hg19.vcf $BRCA_RESOURCES/hg38.fa $LOVD/sharedLOVD_brca12.hg38.vcf
  7. `vcf-sort $LOVD/sharedLOVD_brca12.38.vcf > $LOVD/sharedLOVD_brca12.sorted.38.vcf`
  8. `cp $LOVD/sharedLOVD_brca12.sorted.38.vcf $PIPELINE_INPUT
  
