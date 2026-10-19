.DEFAULT_GOAL := help

# If the .ONESHELL special target appears anywhere in the makefile then all
# recipe lines for each target will be provided to a single invocation of the shell.
.ONESHELL:

CONFIG_PATH := brca_pipeline_cfg.mk

ifneq ("$(wildcard $(CONFIG_PATH))","")
	include $(CONFIG_PATH)
endif

OS := $(shell uname -s)

MKFILE_PATH := $(abspath $(lastword $(MAKEFILE_LIST)))

GIT_BRANCH ?= master

# local (non-Docker) pipeline configuration -- these are the paths/settings
# used to run the real, current Luigi pipeline directly on this host.
# PYTHONPATH_PIPELINE is resolved to an absolute path (relative to this
# Makefile's directory) so it stays correct across the chained `cd`s in the
# test targets below, regardless of a relative default.
PYTHONPATH_PIPELINE ?= ../../code/pipeline
PYTHONPATH_PIPELINE := $(abspath $(PYTHONPATH_PIPELINE))

VCF_ASSEMBLY_OUT_DIR ?= ../../data_out
VCF_ASSEMBLY_SEQ_REPO_DIR ?= /data/seqrepo/latest
VCF_ASSEMBLY_SEQ_REPO_ROOT ?= $(patsubst %/,%,$(dir $(VCF_ASSEMBLY_SEQ_REPO_DIR)))
GENE_CONFIG_PIPELINE ?= $(PYTHONPATH_PIPELINE)/workflow/gene_config_brca_only.txt
VCF_ASSEMBLY_UTA_DB_URL ?= postgresql://anonymous@localhost:50828/uta/uta_20241220
# PostgreSQL schema for pipeline tables; matches PipelineParams.db_schema's
# own default ('pipeline') unless overridden -- e.g. by brca_pipeline_cfg.mk,
# which sets this to the work directory's name for a given release (a schema
# generate_release.py creates in Postgres before invoking this Makefile).
PIPELINE_DB_SCHEMA ?= pipeline

VCF_ANALYSIS_MODULE ?= workflow.variant_analysis
VCF_ANALYSIS_ROOT_TASK ?= workflow.variant_analysis.VariantAnalysis
VCF_ANALYSIS_VEP_SERVER_URL ?= http://localhost:8888
VCF_ANALYSIS_GENOME_FA ?= ../../resources/hg38.fa
# Optional: path to the SpliceAI-scored VCF from the previous release, to
# avoid rescoring already-scored variants. Leave unset for a from-scratch
# run (e.g. no previous release exists yet) to score every variant.
VCF_ANALYSIS_PREVIOUS_SPLICEAI_VCF ?=
ifneq ($(strip $(VCF_ANALYSIS_PREVIOUS_SPLICEAI_VCF)),)
VCF_ANALYSIS_PREVIOUS_SPLICEAI_VCF_FLAG = --previous-spliceai-vcf $(VCF_ANALYSIS_PREVIOUS_SPLICEAI_VCF)
endif

# Whether workflow.variant_assembly.DownloadGnomADCoverage should recompute the gnomAD
# coverage parquets from raw gnomAD coverage data instead of downloading the pre-built
# parquets from brcaexchange.org (the default -- fast, and normally a no-op since
# download-resources already staged the same files into RESOURCES_DIR). Override from
# generate_release.py (--run-gnomad-coverage-estimation) or directly as
# DOWNLOAD_GNOMAD_COVERAGE_COMPUTE_NEW=true.
#
# DownloadGnomADCoverage is a plain luigi.Task (not a luigi.Config like
# PipelineParams), and it lives in workflow.variant_assembly while the root
# task runs out of workflow.variant_analysis -- so luigi only recognizes its
# CLI flag in the fully module-qualified form below, and since compute_new is
# a BoolParameter with implicit (store_true) parsing, the flag must be passed
# bare (no value) and only when true.
DOWNLOAD_GNOMAD_COVERAGE_COMPUTE_NEW ?= false
ifeq ($(strip $(DOWNLOAD_GNOMAD_COVERAGE_COMPUTE_NEW)),true)
DOWNLOAD_GNOMAD_COVERAGE_COMPUTE_NEW_FLAG = --workflow.variant-assembly.DownloadGnomADCoverage-compute-new
endif

LUIGI_ROOT_TASK = $(VCF_ANALYSIS_ROOT_TASK)
LUIGI_WORKERS = 8
LUIGI_LOG_LEVEL ?= WARNING

# run-pipeline's log destination -- a timestamped file directly under this
# release's root (LOG_DIR defaults to WORK_DIR; see brca_pipeline_cfg.mk.j2).
# Override PIPELINE_LOG explicitly (e.g. from generate_release.py) to agree on
# the exact path up front rather than parsing it back out of `make` output --
# its directory is created from wherever it actually points, so overriding
# just this one variable is always enough.
PIPELINE_LOG ?= $(LOG_DIR)/pipeline_$(shell date +%Y%m%d_%H%M%S).log

# build-release's own log destination -- covers everything in the chain
# (checkout, downloads, docker service startup, ...) except run-pipeline's
# Luigi output, which goes to PIPELINE_LOG instead. Same override convention
# as PIPELINE_LOG: set BUILD_RELEASE_LOG explicitly to agree on the path
# up front (e.g. from generate_release.py).
BUILD_RELEASE_LOG ?= $(LOG_DIR)/tmp/build_release_$(shell date +%Y%m%d_%H%M%S).log

LOCAL_PIPELINE_ENV = \
	HGVS_SEQREPO_DIR=$(VCF_ASSEMBLY_SEQ_REPO_DIR) \
	UTA_DB_URL=$(VCF_ASSEMBLY_UTA_DB_URL) \
	PYTHONPATH=$(PYTHONPATH_PIPELINE)

# --gene-config below is DownloadGnomADCoverage's own gene_config parameter (distinct
# from --PipelineParams-gene-config-path above), inherited via @requires into
# ProvideGnomADv41JointVCF/ProvideGnomADv41ExomeVCF/VRSAnnotateGnomADv41 -- without it
# those tasks silently default to gene_config_brca_only.txt and never pick up a
# per-release gene set (verified: this silently dropped PALB2 from gnomAD v4.1 data on
# the palb2 release, since the joint/exome VCFs came back with zero chr16 rows).
# Unlike DOWNLOAD_GNOMAD_COVERAGE_COMPUTE_NEW_FLAG above, this must be passed at the ROOT
# task (unprefixed, not module-qualified as --DownloadGnomADCoverage-gene-config) --
# only the root-level flag actually threads a value down through the whole @requires
# chain, since each level's auto-generated requires() clones its OWN parameter value
# (defaulting per-class otherwise) rather than reading the leaf class's CLI flag directly.
PIPELINE_ARGS := \
	--PipelineParams-gene-config-path $(GENE_CONFIG_PIPELINE) \
	--PipelineParams-seq-repo-dir $(VCF_ASSEMBLY_SEQ_REPO_DIR) \
	--PipelineParams-file-parent-dir $(VCF_ASSEMBLY_OUT_DIR) \
	--PipelineParams-output-dir $(VCF_ASSEMBLY_OUT_DIR) \
	--PipelineParams-db-schema $(PIPELINE_DB_SCHEMA) \
	--gene-config $(GENE_CONFIG_PIPELINE) \
	$(DOWNLOAD_GNOMAD_COVERAGE_COMPUTE_NEW_FLAG) \
	--workers $(LUIGI_WORKERS)

RUN_TASK_CMD = $(LOCAL_PIPELINE_ENV) python -m luigi --module $(VCF_ANALYSIS_MODULE)

init: ## setup config file
	echo "creating file in $(CONFIG_PATH). Please edit to your needs in particular the WORK_DIR!"
	jinja2 -D DATA_DATE=2099-12-31 \
		-D CODE_BASE=$(shell realpath "$(shell dirname $(MKFILE_PATH))/..") \
		-D WORK_DIR='/tmp/brca_workdir' \
		-D RELEASE_TAG='dev' \
		pipeline_running/brca_pipeline_cfg.mk.j2 > $(CONFIG_PATH)

checkout: ## checking out code base; set GIT_COMMIT to pin an exact commit, or GIT_BRANCH (default master) to track a branch's latest commit
	[ -d $(CODE_BASE) ] || git clone $(BRCA_GIT_REPO) $(CODE_BASE)
	cd $(CODE_BASE) && git fetch origin
ifdef GIT_COMMIT
	cd $(CODE_BASE) && git checkout $(GIT_COMMIT)
else
	cd $(CODE_BASE) && git checkout $(GIT_BRANCH) && git merge --ff-only origin/$(GIT_BRANCH)
endif

build-docker: ## build main pipeline docker image
	cd $(CODE_BASE)/pipeline/docker && ./build_docker_image.sh $(RELEASE_TAG) || { echo "$(date) Something wrong with docker image creation"; exit 1; }

download-resources: ## download resources files
	$(CODE_BASE)/pipeline/download_resources_files.sh $(RESOURCES_DIR) >> "$(LOG_DIR)/download_resources_$(RELEASE_TAG).log" 2>&1 || { echo "Downloading resource files failed"; exit 1; }

VEP_CONTAINER    ?= vep-server
VEP_PORT         ?= 8888
VEP_CACHE_DIR    ?= /data/vep_cache
VEP_FASTA        ?= /data/references/hg38.fa
VEP_IMAGE        ?= brca-vep-server

build-vep-server: ## Build the VEP REST server Docker image (skipped if already built)
	docker image inspect $(VEP_IMAGE) >/dev/null 2>&1 || docker build -t $(VEP_IMAGE) vep_server/

setup-vep-cache: ## Download Ensembl VEP homo_sapiens GRCh38 cache (~15 GB)
	mkdir -p $(VEP_CACHE_DIR)
	docker run --rm --user root -v $(VEP_CACHE_DIR):/vep_cache \
	  ensemblorg/ensembl-vep:latest \
	  INSTALL.pl --AUTO c --SPECIES homo_sapiens --ASSEMBLY GRCh38 \
	             --CACHEDIR /vep_cache --NO_HTSLIB

start-vep-server: build-vep-server ## Start VEP REST server on VEP_PORT (default 8888) if not already running
	@echo "Starting VEP REST server..."
	[ $$(docker ps -f name=$(VEP_CONTAINER) | wc -l) -gt 1 ] || \
	docker run -d --name $(VEP_CONTAINER) \
	  --restart unless-stopped \
	  -p $(VEP_PORT):8080 \
	  -v $(VEP_CACHE_DIR):/vep_cache \
	  -v $(VEP_FASTA):/references/hg38.fa:ro \
	  -v $(VEP_FASTA).fai:/references/hg38.fa.fai:ro \
	  -e VEP_FASTA=/references/hg38.fa \
	  $(VEP_IMAGE)

stop-vep-server: ## Stop and remove the VEP REST server container
	docker stop $(VEP_CONTAINER) && docker rm $(VEP_CONTAINER)

start-local-uta: ## starting local uta docker container
	@echo "Starting local UTA docker..."
	[ `docker ps -f name=$(UTA_CONTAINER) | wc -l` -gt 1 ] || \
	docker run -dit --name $(UTA_CONTAINER) \
          -v /tmp:/tmp -v $(UTA_VOLUME):/var/lib/postgresql/data \
          -e $(UTA_PG_PW) -p $(UTA_PORT):5432 \
          $(UTA_DOCKER_IMAGE)

start-seqrepo-rest-service:
	@echo "Starting seqrepo REST service..."
	[ `docker ps -f name="seqrepo-rest-service" | wc -l` -gt 1 ] || \
	utilities/launch_seqrepo.sh

setup-files: ## setup various directories to run pipeline
	[ -f $(CREDENTIALS_PATH) ] || { touch $(CREDENTIALS_PATH); echo "WARNING: files $(CREDENTIALS_PATH) doesn't exist.Creating it"; }
	[ -f $(RELEASE_NOTES_PATH) ] || { touch $(RELEASE_NOTES_PATH) ; echo "WARNING: files $(RELEASE_NOTES_PATH) doesn't exist. Creating it"; }
	[ -d $(RESOURCES_DIR) ] || mkdir -p $(RESOURCES_DIR)
	[ -d $(PRIORS_REFERENCES) ] || mkdir -p $(PRIORS_REFERENCES)
	[ -d $(OUT_DIR) ] || mkdir -p $(OUT_DIR)
	[ -f $(PREVIOUS_RELEASE_DIR)/latest_release.tar.gz ] || wget --progress=dot:giga https://brcaexchange.org/backend/downloads/releases/current_release.tar.gz -O $(PREVIOUS_RELEASE_DIR)/latest_release.tar.gz

download-seqrepo: ## Download seq repo data locally (no Docker; seqrepo pull skips downloading the latest release if we already have it)
	@echo "Download seqrepo data..."
	mkdir -p $(VCF_ASSEMBLY_SEQ_REPO_ROOT)
	[ -d $(LOG_DIR)/tmp ] || mkdir -p $(LOG_DIR)/tmp
	[ -d $(LOG_DIR)/tmp/seqrepo_venv ] || python -m venv $(LOG_DIR)/tmp/seqrepo_venv
	. $(LOG_DIR)/tmp/seqrepo_venv/bin/activate
	pip install -q biocommons.seqrepo
	seqrepo -v --root-directory $(VCF_ASSEMBLY_SEQ_REPO_ROOT) pull --update-latest

run-pipeline: ## Run the current luigi pipeline locally (no Docker); output logged to PIPELINE_LOG, not the console
	@echo "Running Luigi pipeline..."
	mkdir -p $(dir $(PIPELINE_LOG))
	@echo "Pipeline log: $(PIPELINE_LOG)"
	[ -d $(LOG_DIR)/tmp/luigi_venv ] || python -m venv $(LOG_DIR)/tmp/luigi_venv
	. $(LOG_DIR)/tmp/luigi_venv/bin/activate
	pip install -q -r requirements.txt
	# spliceai 1.3.1 (latest) uses pkg_resources which was removed from setuptools 71+.
	# Same patch as pipeline/docker/Dockerfile: fall back to importlib when pkg_resources
	# is missing. find_spec (not import) locates the package without executing its
	# pkg_resources-dependent __init__.py, so this works whether or not the venv already
	# has pkg_resources; the ^-anchored regex only matches the original unindented import
	# line, so re-running this on an already-patched file is a no-op instead of corrupting it.
	SPLICEAI_DIR=$$(python -c "import importlib.util, os; print(os.path.dirname(importlib.util.find_spec('spliceai').origin))")
	python3 -c "import pathlib, re; p = pathlib.Path('$$SPLICEAI_DIR/__init__.py'); src = p.read_text(); src = re.sub(r'^from pkg_resources import get_distribution', 'try:\n    from pkg_resources import get_distribution\nexcept ImportError:\n    from importlib.metadata import distribution as get_distribution', src, count=1, flags=re.MULTILINE); p.write_text(src)"
	python3 -c "import pathlib, re; p = pathlib.Path('$$SPLICEAI_DIR/utils.py'); src = p.read_text(); src = re.sub(r'^from pkg_resources import resource_filename', 'try:\n    from pkg_resources import resource_filename\nexcept ImportError:\n    import importlib.resources\n    def resource_filename(package, resource):\n        return str(importlib.resources.files(package).joinpath(resource))', src, count=1, flags=re.MULTILINE); p.write_text(src)"
	PATH=$(LOG_DIR)/tmp/luigi_venv/bin:$(HOME)/bin:$(PATH) \
	$(LOCAL_PIPELINE_ENV) \
	python -m luigi --module $(VCF_ANALYSIS_MODULE) $(VCF_ANALYSIS_ROOT_TASK) \
		$(PIPELINE_ARGS) \
		--vep-server-url $(VCF_ANALYSIS_VEP_SERVER_URL) \
		--genome-fa $(VCF_ANALYSIS_GENOME_FA) \
		$(VCF_ANALYSIS_PREVIOUS_SPLICEAI_VCF_FLAG) \
		--local-scheduler \
		--log-level $(LUIGI_LOG_LEVEL) \
		> $(PIPELINE_LOG) 2>&1

# If the first make target has arguments (run-task and others)...
# https://stackoverflow.com/questions/2214575/passing-arguments-to-make-run
ifeq ($(firstword $(MAKECMDGOALS)),$(filter $(firstword $(MAKECMDGOALS)), run-task force-run-task clean-files-from))
  # use the rest as arguments for "run"
  RUN_ARGS := $(wordlist 2,$(words $(MAKECMDGOALS)),$(MAKECMDGOALS))
  # ...and turn them into do-nothing targets
  $(eval $(RUN_ARGS):;@:)
endif

run-task: ## Running a specific luigi task locally (no Docker)
	$(RUN_TASK_CMD) $(RUN_ARGS) $(PIPELINE_ARGS) --local-scheduler

force-run-task: ## Running a specific task locally (deleting its dependencies first)
	$(LOCAL_PIPELINE_ENV) python $(PYTHONPATH_PIPELINE)/workflow/prune_output_files.py $(VCF_ANALYSIS_MODULE) $(LUIGI_ROOT_TASK) $(RUN_ARGS) $(PIPELINE_ARGS)
	$(RUN_TASK_CMD) $(RUN_ARGS) $(PIPELINE_ARGS) --local-scheduler

clean-files-from: ## Delete files generated by tasks on the task path between the given task and $(LUIGI_ROOT_TASK)
	$(LOCAL_PIPELINE_ENV) python $(PYTHONPATH_PIPELINE)/workflow/prune_output_files.py $(VCF_ANALYSIS_MODULE) $(LUIGI_ROOT_TASK) $(RUN_ARGS) $(PIPELINE_ARGS)

show-luigi-graph: ## print some representation of the luigi compute graph on the console
	HGVS_SEQREPO_DIR=$(VCF_ASSEMBLY_SEQ_REPO_DIR) \
	UTA_DB_URL=$(VCF_ASSEMBLY_UTA_DB_URL) \
	PYTHONPATH="$(PYTHONPATH_PIPELINE):$(PYTHONPATH_PIPELINE)/workflow" \
	luigi-deps-tree --module $(VCF_ANALYSIS_MODULE) $(PIPELINE_ARGS) $(LUIGI_ROOT_TASK) | sed -E "s/([^{])\{.*\}([^}])/\1/g"

test: ## Running pipeline unit tests locally (no Docker)
	cd $(PYTHONPATH_PIPELINE)/data && bash getdata && cd $(PYTHONPATH_PIPELINE) && PYTHONPATH=$(PYTHONPATH_PIPELINE) pytest --ignore=splicing/

test-coverage: ## Running pipeline unit tests with coverage information locally (no Docker)
	cd $(PYTHONPATH_PIPELINE)/data && bash getdata && cd $(PYTHONPATH_PIPELINE) && PYTHONPATH=$(PYTHONPATH_PIPELINE) pytest --cov --ignore=splicing/ && coverage html --include="$(PYTHONPATH_PIPELINE)/*" --omit="*/test_*"


build-release: ## create new data release; output logged to BUILD_RELEASE_LOG, not the console (run-pipeline's Luigi output still goes to its own PIPELINE_LOG)
	mkdir -p $(dir $(BUILD_RELEASE_LOG))
	@echo "Build-release log: $(BUILD_RELEASE_LOG)"
	$(MAKE) PIPELINE_LOG=$(PIPELINE_LOG) checkout setup-files download-resources download-seqrepo start-local-uta start-seqrepo-rest-service start-vep-server run-pipeline variants-by-source > $(BUILD_RELEASE_LOG) 2>&1

variants-by-source: ## postprocessing: compute statistics for changes with respect to the last release
	PYTHONPATH=$(PYTHONPATH_PIPELINE) python $(PYTHONPATH_PIPELINE)/utilities/variantsBySource.py -i $(VCF_ASSEMBLY_OUT_DIR)/output/release/built_with_change_types.tsv -c true

# TODO: GenerateReleaseNotes has no equivalent task yet in workflow.variant_analysis;
# re-enable once release-notes/archive generation is ported to the new schema pipeline.
# prune-release-notes-output: ## postprocessing: removes some files, s.t. luigi tasks can be triggered to include updated release notes into release archive
# 	make clean-files-from GenerateReleaseNotes

include-release-notes: run-pipeline ## postprocessing: include updated release notes into archive

cleanup-failed: ## postprocessing: cleaning up some dangling FAILED files, which may be confusing
	find $(OUT_DIR) -name 'FAILED*' -exec rm -i {} \;

tag-release: ## postprocessing: tags and pushes git commit used for the data release
	cd $(CODE_BASE)
	git tag -m $(RELEASE_TAG) -a $(RELEASE_TAG) `git rev-parse HEAD`
	git push origin $(RELEASE_TAG)

push-docker: ## postprocessing: pushes docker image to docker hub
	docker login
	docker push brcachallenge/brca-exchange-pipeline:$(RELEASE_TAG)
	docker logout

copy-to-previous: ## postprocessing: copy new release tar to previous release folder
# 	take latest, in case there should be several tar.gz from some reason
	cp `ls -t $(OUT_DIR)/release-*.tar.gz | head -n 1` $(PREVIOUS_RELEASE_DIR)
	rm $(PREVIOUS_RELEASE_DIR)/latest_release.tar.gz
	ln -s `ls -t $(PREVIOUS_RELEASE_DIR)/release-*.tar.gz | head -n 1` $(PREVIOUS_RELEASE_DIR)/latest_release.tar.gz

post-release-cmds: cleanup-failed include-release-notes push-docker tag-release copy-to-previous ## postprocessing: do all the necessary postprocessing to wrap up a data release

setup-data-from-latest-release-tar: setup-files ## sets up brca output dir with data contained in release archive from last release (only data from variant merging onwards)
	tar -C $(OUT_DIR) -zxf $(PREVIOUS_RELEASE_PATH)

setup-dev-env: build-docker setup-files download-resources ## setup development environment

run-variant-assembly: ## Run the VariantAssembly luigi task locally (no Docker)
	$(LOCAL_PIPELINE_ENV) \
	python -m luigi --module workflow.variant_assembly workflow.variant_assembly.VCFAssembly \
		$(PIPELINE_ARGS) \
		--local-scheduler

run-variant-analysis: ## Run the full pipeline (VariantAssembly + VariantAnalysis: VEP/BayesDel/SpliceAI/priors/popfreq) locally (no Docker)
	PATH=$(HOME)/bin:$(PATH) \
	$(LOCAL_PIPELINE_ENV) \
	python -m luigi --module $(VCF_ANALYSIS_MODULE) $(VCF_ANALYSIS_ROOT_TASK) \
		$(PIPELINE_ARGS) \
		--vep-server-url $(VCF_ANALYSIS_VEP_SERVER_URL) \
		--genome-fa $(VCF_ANALYSIS_GENOME_FA) \
		$(VCF_ANALYSIS_PREVIOUS_SPLICEAI_VCF_FLAG) \
		--local-scheduler

clean-pyc: ## remove Python file artifacts
	find . -name '*.pyc' -exec rm -f {} +
	find . -name '*.pyo' -exec rm -f {} +
	find . -name '*~' -exec rm -f {} +
	find . -name '__pycache__' -exec rm -fr {} +


define PRINT_HELP_PYSCRIPT
import re, sys

for line in sys.stdin:
        match = re.match(r'^([a-zA-Z_-]+):.*?## (.*)$$', line)
        if match:
                target, help = match.groups()
                print("%-20s %s" % (target, help))
endef
export PRINT_HELP_PYSCRIPT

help:
	@python -c "$$PRINT_HELP_PYSCRIPT" < $(MAKEFILE_LIST)
