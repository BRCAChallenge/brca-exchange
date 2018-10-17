This directory contains code to create a .bigBed file for a UCSC genome browser
track hub from the brcaexchange tarball.

First, get the release archive you would like to use. For example, if you're working with the latest production release:
```
wget https://brcaexchange.org/backend/downloads/releases/current_release.tar.gz
tar -zxvf current_release.tar.gz
```

Convert the release to a .bed file:
```
python brcaToBed.py
```

Sort the data for both hg19 and hg38:
```
sort -k1,1 -k2,2 brcaExchange.hg19.bed -o brcaExchange.hg19.bed
sort -k1,1 -k2,2 brcaExchange.hg38.bed -o brcaExchange.hg38.bed
```

Make sure you have the correct bedToBigBed script as found in the UCSC kent tool collection: http://hgdownload.cse.ucsc.edu/admin/exe/, then run:
```
chmod a+x bedToBigBed
./bedToBigBed -type=bed9+ -as=brcaExchange.as -tab brcaExchange.hg19.bed hg19.chrom.sizes hg19/brcaExchange.bb
./bedToBigBed -type=bed9+ -as=brcaExchange.as -tab brcaExchange.hg38.bed hg38.chrom.sizes hg38/brcaExchange.bb
```

Move the output files to a single directory:
```
mkdir trackhubs
cp -R hub.txt genomes.txt hg19 hg38 ./trackhubs/
```

Move the folder to the correct location on the brcaexchange server:
```
scp -prq trackhubs/. brca@brcaexchange.org:/var/www/html/production/trackhubs/
```
