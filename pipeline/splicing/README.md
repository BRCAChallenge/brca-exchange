# Splice Site Analysis Pipeline
Generates splice site affect from a list of variants

# Docker
The pipeline has been dockerized with facilities for downloading required reference files, unit and functional tests.

# Notes
The file brcaPase.py contains the class for parsing through the .tsv file from the BRCA Exchange website. The class takes the mutations in the Ref column and creates a new variant sequence for MaxEntScan to parse through and generate MaxEntScores for variants. 
