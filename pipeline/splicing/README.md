# Splice Site Analysis Pipeline
Calculates the prior probability of pathogenicity or a prior ENGIMA classification based on variant type and variant location

# Running
The pipeline has been dockerized with integrated support for downloading reference files and running unit and functional tests.

Download reference files:

	docker run -it --rm \
		--user=`id -u`:`id -g` \
		-v <host path to references>:/references \
		brcachallenge/splicing-pipeline references

Output:

	Downloading references. This may take 1-2 hours...
	--2018-06-25 16:26:32--  http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
	...
	Validating MD5 of references...
	/references/hg38.fa.gz: OK
	/references/homo_sapiens_vep_92_GRCh38.tar.gz: OK
	Unpacking VEP references...
	Calculating one variant to force installation of references...
	connect
	DEBUG seqdb._create_seqLenDict: Building sequence length index...
	[fai_load] build FASTA index.
	variant 1 complete
	Reference download and installation complete.

Downloads genome and Ensembl VEP references, verifies their MD5, unpacks and calculates a single variant to cause final VEP and reference installation. After this step the references directory can be mapped read-only. The download can take up to 2 hours and the resulting references folder will occupy 18GBs.

Run short test:

	docker run -it --rm \
		--user=`id -u`:`id -g` \
		-v <host path to references>:/references:ro \
		brcachallenge/splicing-pipeline test short

Runs integrated unit tests via pytest and then calculates priors for 8 variants and verifies the MD5 of the output. You can run a longer test of 88 variants via 'test long'.

Output:

	================ test session starts ================
	...
	=========== 267 passed in 1.14 seconds =============
	variant 1 complete
	...
	variant 8 complete
	/tmp/priors_short.tsv: OK

Calculate:

	docker run --rm -it \
		--user=`id -u`:`id -g` \
		-v <host path to references>:/references:ro \
		-v <host path to input and output data>:/data \
		brcachallenge/splicing-pipeline calc <input variants tsv> <output priors tsv>

See the tests/ folder in this repo for examples of input and output files.

# Developing
To develop the pipeline you can either install all of its requirements in a local virtualenv or debug within the running docker. The later approach is highly recommended:

Clone this repo, change into the pipelines/splicing directory and then:

	docker run --rm -it \
		--entrypoint /bin/bash \
		--user=`id -u`:`id -g` \
		-v <host path to references>:/references \
		-v `pwd`:/app \
		brcachallenge/splicing-pipeline

You will end up with a shell inside the container but running the code mapped from your local directory. You can then edit the source outside of the docker and run/test the code using the shell inside the docker.

# Notes
The file brcaPase.py contains the class for parsing through the .tsv file from the BRCA Exchange website. The class takes the mutations in the Ref column and creates a new variant sequence for MaxEntScan to parse through and generate MaxEntScores for variants.

Run the docker with -h to see the various other options:

	docker run --rm -it brcachallenge/splicing-pipeline -h
