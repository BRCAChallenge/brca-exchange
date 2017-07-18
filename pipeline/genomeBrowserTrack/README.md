This directory contains code to create a .bigBed file for a UCSC genome browser
track hub from the brcaexchange tarball.

On the hgwdev server at UCSC, Max runs the code like this:

    wget http://brcaexchange.org/backend/downloads/releases/current_release.tar.gz
    python brcaToBed.py 
    sort -k1,1 -k2,2 brcaExchange.hg19.bed -o brcaExchange.hg19.bed
    sort -k1,1 -k2,2 brcaExchange.hg38.bed -o brcaExchange.hg38.bed

    bedToBigBed -type=bed9+ -as=brcaExchange.as -tab brcaExchange.hg19.bed hg19.chrom.sizes hg19/brcaExchange.bb
    bedToBigBed -type=bed9+ -as=brcaExchange.as -tab brcaExchange.hg38.bed hg38.chrom.sizes hg38/brcaExchange.bb
    cp brcaExchange.bb ~/public_html/immuno/track/hub/hg19/

bedSort and bedToBigBed are tools from the UCSC kent tool collection: http://hgdownload.cse.ucsc.edu/admin/exe/
