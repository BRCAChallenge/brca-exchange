# get the indexes of mutations
DATA=/hive/groups/cgl/brca/phase1/data/umd/
wget http://www.umd.be/BRCA1/4DACTION/W_DMDT1/1 -O $DATA/umdBrca1.html
wget http://www.umd.be/BRCA2/4DACTION/W_DMDT1/1 -O $DATA/umdBrca2.html

# get the IDs
grep ../../4DACTION/WV/[0-9]* $DATA/umdBrca1.html  | egrep -o "/[0-9]{1,9}'" | tr -d '/' | tr -d "'" > brca1Ids.txt
grep ../../4DACTION/WV/[0-9]* $DATA/umdBrca2.html  | egrep -o "/[0-9]{1,9}'" | tr -d '/' | tr -d "'" > brca2Ids.txt

# get all pages

mkdir -p $DATA/brca1Html $DATA/brca2Html
for i in `cat brca1Ids.txt`; do wget http://www.umd.be/BRCA2/4DACTION/WV/$i -O /hive/groups/cgl/brca/phase1/data/umdHtml/brca1Html/; done
for i in `cat brca2Ids.txt`; do wget http://www.umd.be/BRCA2/4DACTION/WV/$i -O /hive/groups/cgl/brca/phase1/data/umdHtml/brca2Html/; done

# some failed
grep "A runtime error" brca1Html/*.html | cut -d/ -f2 | cut -d. -f1 > brc1Failed.txt
grep "A runtime error" brca2Html/*.html | cut -d/ -f2 | cut -d. -f1 > brc2Failed.txt 

# extract data from pages as a tsv file
python scrapeUmdDetails.py $DATA/umdBrca1.html  $DATA/brca1Html/*.html > brca1.tsv
python scrapeUmdDetails.py $DATA/umdBrca2.html  $DATA/brca2Html/*.html > brca2.tsv
