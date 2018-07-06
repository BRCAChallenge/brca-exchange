build:
	# Build and tag the image prefixed by user to avoid conflicts on shared machines
	docker build -t $(USER)-splicing-pipeline .

push:
	# Tag and push our current docker build to dockerhub
	docker tag $(USER)-splicing-pipeline:latest brcachallenge/splicing-pipeline:latest
	docker push brcachallenge/splicing-pipeline:latest

debug:
	# Run the docker mapping out of local directory for development
	docker run --rm -it \
		--entrypoint /bin/bash \
		--user=`id -u`:`id -g` \
		-v `readlink -f ~/references/splicing`:/references \
		-v `pwd`:/app \
		$(USER)-splicing-pipeline

test:
	# Run a short test from within the docker container
	python calcVarPriors.py test short

references:
	# Download and install references
	docker run --rm -it \
		--user=`id -u`:`id -g` \
		-v `readlink -f ~/references/splicing`:/references \
		$(USER)-splicing-pipeline references

short:
	# Run a short test using the docker container
	docker run --rm -it \
		--user=`id -u`:`id -g` \
		-v `readlink -f ~/references/splicing`:/references:ro \
		$(USER)-splicing-pipeline test short

long:
	# Run a long test using the docker container with 22 cores
	docker pull brcachallenge/splicing-pipeline:latest
	docker run --rm -it \
		--user=`id -u`:`id -g` \
		-v `readlink -f ~/references/splicing`:/references:ro \
		brcachallenge/splicing-pipeline:latest test long -p 22 

vep:
	# Run a single request through local vep
	vep \
		--cache --dir_cache /references/vep/ --no_stats --offline \
		--fasta /references/hg38.fa \
		--output_file STDOUT \
		--input_data 13:32346895-32346895:1/G

profile:
	# Run a short set of variants with profiling
	python -m cProfile -o profile calcVarPriors.py \
		-b "enigma" \
		-v mod_res_dn_brca20160525.txt \
		-g /references/hg38.fa \
		-t /references/refseq_annotation.hg38.gp \
		-i tests/variants_short.tsv \
		-o tests/priors_short.tsv

upstream:
	# Merge upstream into this fork
	git pull https://github.com/BRCAChallenge/brca-exchange master
