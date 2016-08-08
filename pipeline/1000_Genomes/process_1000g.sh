FOLDER="/cluster/home/mollyzhang/common_brca_folder/release1.0/data/1000_genomes"

echo downloading 1000 genome variant data from ftp....takes about 5 minutes....
echo -e "\ndata location: $FOLDER"
ftp -n < ftp_data_download.txt

echo unzipping those data ...takes about 10 mintues...
gzip -d $FOLDER/ftp_download/chr13_1000g_GRCh37_vcf.gz
gzip -d $FOLDER/ftp_download/chr17_1000g_GRCh37_vcf.gz

echo remove sample columns ... another 5 minutes ...
cut -f1-9 $FOLDER/ftp_download/chr13_1000g_GRCh37_vcf > $FOLDER/ftp_download/chr13_1000g_GRCh37_no_sample.vcf
cut -f1-9 $FOLDER/ftp_download/chr17_1000g_GRCh37_vcf > $FOLDER/ftp_download/chr17_1000g_GRCh37_no_sample.vcf

echo extract BRCA gene region from chr13 and chr17
head -440648 $FOLDER/ftp_download/chr13_1000g_GRCh37_no_sample.vcf | tail -2286 > $FOLDER/ftp_download/chr13_brca2_1000g_GRCh37.vcf
head -1146432 $FOLDER/ftp_download/chr17_1000g_GRCh37_no_sample.vcf | tail -2052 > $FOLDER/ftp_download/chr17_brca1_1000g_GRCh37.vcf

echo cacatenate header, brca1 and brca2 files to generate final output
grep "^#" $FOLDER/ftp_download/chr17_1000g_GRCh37_no_sample.vcf > $FOLDER/ftp_download/header
cat $FOLDER/ftp_download/header $FOLDER/ftp_download/chr13_brca2_1000g_GRCh37.vcf $FOLDER/ftp_download/chr17_brca1_1000g_GRCh37.vcf > $FOLDER/1000g_brca12_GRCh37.vcf 

echo Done!
