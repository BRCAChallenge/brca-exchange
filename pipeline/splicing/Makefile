# Host path where references should be stored - symlink tolerant
REFERENCES = $(shell readlink -f ~/references/splicing)
IMAGE_NAME = $(USER)-splicing-pipeline

default: build references short

build:
	# Build and tag the image prefixed by user to avoid conflicts on shared machines
	docker build -t $(IMAGE_NAME) .

push:
	# Tag and push our current docker build to dockerhub
	docker tag $(USER)-splicing-pipeline:latest brcachallenge/splicing-pipeline:latest
	docker push brcachallenge/splicing-pipeline:latest

debug:
	# Run the docker mapping out of local directory for development
	docker run --rm -it \
		--entrypoint /bin/bash \
		--user=`id -u`:`id -g` \
		-v `readlink -f $(REFERENCES)`:/references \
		-v `pwd`:/app \
	  $(IMAGE_NAME)

test:
	# Run a short test from within the docker container
	python calcVarPriors.py test short

unit_tests_docker:
	docker run --rm -it \
		--user=`id -u`:`id -g` \
    --entrypoint pytest \
		$(IMAGE_NAME)


references:
	# Download and install references
	mkdir -p $(REFERENCES)
	docker run --rm -it \
		--user=`id -u`:`id -g` \
		-v `readlink -f $(REFERENCES)`:/references \
		$(IMAGE_NAME) references

short:
	# Run a short test using the docker container
	docker run --rm -it \
		--user=`id -u`:`id -g` \
		-v `readlink -f $(REFERENCES)`:/references:ro \
		$(IMAGE_NAME) test short

long:
	# Run a long test using the docker container with 22 cores
	docker pull brcachallenge/splicing-pipeline:latest
	docker run --rm -it \
		--user=`id -u`:`id -g` \
		-v `readlink -f $(REFERENCES)`:/references:ro \
		brcachallenge/splicing-pipeline:latest --processes 22 test long

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
