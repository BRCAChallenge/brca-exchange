This directory contains code to create a .bigBed file for a UCSC genome browser
track hub from the brcaexchange tarball.

On the hgwdev server at UCSC, Max runs the code like this:

    wget http://brcaexchange.org/backend/downloads/releases/current_release.tar.gz
    python brcaToBed.py 
    bedSort brcaExchange.bed brcaExchange.bed
    bedToBigBed -type=bed9+ -as=brcaExchange.as -tab brcaExchange.bed /scratch/data/hg19/chrom.sizes brcaExchange.bb
    cp brcaExchange.bb ~/public_html/immuno/track/hub/hg19/

bedSort and bedToBigBed are tools from the UCSC kent tool collection: http://hgdownload.cse.ucsc.edu/admin/exe/
