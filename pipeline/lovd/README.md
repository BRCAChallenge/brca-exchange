##### Webscrap exLOVD (`http://hci-exlovd.hci.utah.edu/`)
  1. Extract variant data to .txt flat file format: e.g. `extract_data.py -u http://hci-exlovd.hci.utah.edu/ -l BRCA1 BRCA2 -o output_directory`
  2. Convert extracted flat file to vcf format `./lovd2vcf -i output_directory/BRCA1.txt -o exLOVD_brca1.vcf -a exLOVDAnnotation -b 1 > stdout.brca1.txt`
  3. `./lovd2vcf -i output_directory/BRCA2.txt -o exLOVD_brca2.vcf -a exLOVDAnnotation -b 2 > stdout.brca2.txt`
  4. `vcf-concat exLOVD_brca1.vcf exLOVD_brca2.vcf > exLOVD_brca12.vcf`
  5. `vcf-sort exLOVD_brca12.vcf > exLOVD_brca12.sorted.vcf`
  6. `bgzip exLOVD_brca12.sorted.vcf`
  7. `bcftools index exLOVD_brca12.sorted.vcf.gz`
  8. `scp -i ~/brca_server.pem exLOVD_brca12.sorted.vcf.gz ubuntu@ec2-54-148-207-224.us-west-2.compute.amazonaws.com:/srv/ga4gh/brca_data/variants/ex-lovd`
  9. `scp -i ~/brca_server.pem exLOVD_brca12.sorted.vcf.gz.csi ubuntu@ec2-54-148-207-224.us-west-2.compute.amazonaws.com:/srv/ga4gh/brca_data/variants/ex-lovd`

##### Webscrap sharedLOVD (`http://databases.lovd.nl/shared/`)
  1. `extract_data.py -u http://databases.lovd.nl/shared/ -l BRCA1 BRCA2 -o output_directory2`
  2. `./lovd2vcf -i output_directory2/BRCA1.txt -o sharedLOVD_brca1.vcf -a sharedLOVDAnnotation -b 1`
  3. `./lovd2vcf -i output_directory2/BRCA2.txt -o sharedLOVD_brca2.vcf -a sharedLOVDAnnotation -b 2`
  4. `vcf-concat sharedLOVD_brca1.vcf sharedLOVD_brca2.vcf > sharedLOVD_brca12.vcf`
  5. `vcf-sort sharedLOVD_brca12.vcf > sharedLOVD_brca12.sorted.vcf`
  6. `bgzip sharedLOVD_brca12.sorted.vcf`
  7. `bcftools index sharedLOVD_brca12.sorted.vcf.gz`
  8. `scp -i ~/brca_server.pem sharedLOVD_brca12.sorted.vcf.gz ubuntu@ec2-54-148-207-224.us-west-2.compute.amazonaws.com:/srv/ga4gh/brca_data/variants/lovd`
  9. `scp -i ~/brca_server.pem sharedLOVD_brca12.sorted.vcf.gz.csi ubuntu@ec2-54-148-207-224.us-west-2.compute.amazonaws.com:/srv/ga4gh/brca_data/variants/lovd`
  
