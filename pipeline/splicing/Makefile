build:
	# Build and tag the image prefixed by user to avoid conflicts on shared machines
	docker build -t $(USER)-splicing-pipeline .

push:
	# Tag and push our current docker build to dockerhub
	docker tag $(USER)-splicing-pipeline brcachallenge/splicing-pipeline

debug:
	# Run the docker mapping out of local directory for development
	docker run --rm -it \
		--entrypoint /bin/bash \
		--user=`id -u`:`id -g` \
		-v `readlink -f ~/references/splicing`:/references \
		-v `pwd`:/app \
		$(USER)-splicing-pipeline

references:
	# Download and verify references - to be called from within the docker
	echo "Downloading references. This may take 1-2 hours..."
	wget -N -P /references http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
	wget -N -P /references ftp://ftp.ensembl.org/pub/release-92/variation/VEP/homo_sapiens_vep_92_GRCh38.tar.gz
	echo "Validating MD5 of references..."
	md5sum -c references.md5
	gunzip -k -f /references/hg38.fa.gz
	mkdir -p /references/vep
	echo "Unpacking VEP references..."
	tar -xzf /references/homo_sapiens_vep_92_GRCh38.tar.gz -C /references/vep
	echo "Running vep to create fasta index and verify a single request..."
	vep \
		--cache --dir_cache /references/vep/ --no_stats --offline \
		--fasta /references/hg38.fa \
		--output_file STDOUT \
		--input_data 13:32346895-32346895:1/G | tail -2 | md5sum -c vep.md5

vep:
	# Run a single request through local vep
	vep \
		--cache --dir_cache /references/vep/ --no_stats --offline \
		--fasta /references/hg38.fa \
		--output_file STDOUT \
		--input_data 13:32346895-32346895:1/G

test:
	# Run automated unit tests
	pytest

profile:
	# Run a short set of variants with profiling
	python -m cProfile -o profile calcVarPriors.py \
		-b "enigma" \
		-v mod_res_dn_brca20160525.txt \
		-g /references/hg38.fa \
		-t /references/refseq_annotation.hg38.gp \
		-i tests/variants_short.tsv \
		-o tests/priors_short.tsv

short:
	# Run a short set of variants
	python calcVarPriors.py \
		-b "enigma" \
		-v mod_res_dn_brca20160525.txt \
		-g /references/hg38.fa \
		-t /references/refseq_annotation.hg38.gp \
		-i tests/variants_short.tsv \
		-i /tmp/short_output

		# -o tests/priors_short.tsv

long:
	# Run a longer set of variants
	python calcVarPriors.py \
		-b "enigma" \
		-v mod_res_dn_brca20160525.txt \
		-g /references/hg38.fa \
		-t /references/refseq_annotation.hg38.gp \
		-i tests/variants_long.tsv \
		-o tests/priors_long.tsv

upstream:
	# Merge upstream into this fork
	git pull https://github.com/ccdlott/brca-exchange.git splicingUpdates
