.DEFAULT_GOAL := help

CONFIG_PATH := brca_pipeline_cfg.mk

ifneq ("$(wildcard $(CONFIG_PATH))","")
	include $(CONFIG_PATH)
endif

OS := $(shell uname -s)

MKFILE_PATH := $(abspath $(lastword $(MAKEFILE_LIST)))

LUIGI_ROOT_TASK = GenerateReleaseArchive
LUIGI_WORKERS = 8

SEQ_REPO_DIR_DOCKER := /files/resources/seq_repo

GENE_CONFIG_FILE := /opt/brca-exchange/pipeline/luigi/gene_config_brca_only.txt

init: ## setup config file
	echo "creating file in $(CONFIG_PATH). Please edit to your needs in particular the WORK_DIR!"
	jinja2 -D DATA_DATE=2099-12-31 \
		-D CODE_BASE=$(shell realpath "$(shell dirname $(MKFILE_PATH))/..") \
		-D WORK_DIR='/tmp/brca_workdir' \
		-D RELEASE_TAG='dev' \
		pipeline_running/brca_pipeline_cfg.mk.j2 > $(CONFIG_PATH)

checkout: ## checking out code base
	[ -d $(CODE_BASE) ] || git clone $(BRCA_GIT_REPO) $(CODE_BASE)
	cd $(CODE_BASE) && git checkout $(GIT_COMMIT)

build-docker: ## build main pipeline docker image
	cd $(CODE_BASE)/pipeline/docker && ./build_docker_image.sh $(RELEASE_TAG) || { echo "$(date) Something wrong with docker image creation"; exit 1; }

download-resources: ## download resources files
	$(CODE_BASE)/pipeline/download_resources_files.sh $(RESOURCES_DIR) >> "$(LOG_DIR)/download_resources_$(RELEASE_TAG).log" 2>&1 || { echo "Downloading resource files failed"; exit 1; }

start-local-uta: ## starting local uta docker container
	[ `docker ps -f name=$(UTA_CONTAINER) | wc -l` -gt 1 ] || docker run -dit --name $(UTA_CONTAINER) -p $(UTA_PORT):5432 $(UTA_DOCKER_IMAGE)

.ONESHELL:
setup-files: ## setup various directories to run pipeline
	[ -f $(CREDENTIALS_PATH) ] || { touch $(CREDENTIALS_PATH); echo "WARNING: files $(CREDENTIALS_PATH) doesn't exist.Creating it"; }
	[ -f $(RELEASE_NOTES_PATH) ] || { touch $(RELEASE_NOTES_PATH) ; echo "WARNING: files $(RELEASE_NOTES_PATH) doesn't exist. Creating it"; }
	[ -d $(RESOURCES_DIR) ] || mkdir -p $(RESOURCES_DIR)
	[ -d $(PRIORS_REFERENCES) ] || mkdir -p $(PRIORS_REFERENCES)
	[ -d $(OUT_DIR) ] || mkdir -p $(OUT_DIR)
	[ -f $(PREVIOUS_RELEASE_DIR)/latest_release.tar.gz ] || wget https://brcaexchange.org/backend/downloads/releases/current_release.tar.gz -O $(PREVIOUS_RELEASE_DIR)/latest_release.tar.gz

.ONESHELL:
setup-lovd: ## setting up LOVD data via pipeline machine, if not running on pipeline machine, as the access is IP restricted. Requires access to main pipeline machine
	if [ `hostname -s` != "brcaexchange-dev" ]; then
		mkdir -p $(OUT_DIR)/LOVD
		ssh brca@brcaexchange-dev.gi.ucsc.edu curl https://databases.lovd.nl/shared/export/BRCA > $(OUT_DIR)/LOVD/BRCA.txt
	fi


ifeq ($(OS), Darwin)
DOCKER_GRP = 0
else
DOCKER_GRP = `stat -c '%g' /var/run/docker.sock`
endif

COMMON_DOCKER_ARGS = --rm -u `id -u ${USER}`:$(DOCKER_GRP) \
	-e "DATA_DATE=$(DATA_DATE)" \
        -e "UTA_DB_URL=postgresql://anonymous@0.0.0.0:$(UTA_PORT)/uta/uta_$(UTA_RELEASE_DATE)" \
	-e "HGVS_SEQREPO_DIR=$(SEQ_REPO_DIR_DOCKER)/latest" \
  -e "PYTHONPATH=/opt/brca-exchange/pipeline" \
	--network host \
	-v $(RESOURCES_DIR):/files/resources \
	-v $(OUT_DIR):/files/data \
	-v $(CREDENTIALS_PATH):/opt/luigi_pipeline_credentials.cfg \
	-v $(PREVIOUS_RELEASE_PATH):/files/previous_release.tar.gz \
	-v $(RELEASE_NOTES_PATH):/files/release_notes.txt \
	-v $(CODE_BASE):/opt/brca-exchange \
	-v /var/run/docker.sock:/var/run/docker.sock

download-seqrepo: ## Download seq repo data
	docker run $(COMMON_DOCKER_ARGS) $(PIPELINE_IMAGE) bash -c "mkdir -p $(SEQ_REPO_DIR_DOCKER) && seqrepo --root-directory $(SEQ_REPO_DIR_DOCKER) pull && seqrepo --root-directory $(SEQ_REPO_DIR_DOCKER) update-latest"

run-pipeline: ## running entire pipeline
	TS=`date +%Y%m%d_%H%M%S`
	LOG_FILE=$(LOG_DIR)/pipeline_run_$${TS}.log
	echo "Log files are in $${LOG_FILE}"
	docker run $(COMMON_DOCKER_ARGS) $(PIPELINE_IMAGE) /opt/brca-exchange/pipeline/docker/run_luigi.sh $(PRIORS_REFERENCES) $(OUT_DIR) $(PRIORS_IMAGE) $(LUIGI_ROOT_TASK) $(LUIGI_WORKERS) $(GENE_CONFIG_FILE) > $${LOG_FILE} 2>&1

run-interactive: ## starting docker container in interactive mode
	docker run -it $(COMMON_DOCKER_ARGS) $(PIPELINE_IMAGE) bash

# If the first argument is "run-task"...
# https://stackoverflow.com/questions/2214575/passing-arguments-to-make-run
ifeq ($(firstword),$(filter $(firstword), run-task, force-run-task, clean-files-from))
  # use the rest as arguments for "run"
  RUN_ARGS := $(wordlist 2,$(words $(MAKECMDGOALS)),$(MAKECMDGOALS))
  # ...and turn them into do-nothing targets
  $(eval $(RUN_ARGS):;@:)
endif

RUN_TASK_CMD := docker run $(COMMON_DOCKER_ARGS) $(PIPELINE_IMAGE) /opt/brca-exchange/pipeline/docker/run_luigi.sh $(PRIORS_REFERENCES) $(OUT_DIR) $(PRIORS_IMAGE)

run-task: ## Running a specific task
	$(RUN_TASK_CMD) $(RUN_ARGS) $(LUIGI_WORKERS) $(GENE_CONFIG_FILE)

PIPELINE_ARGS := --PipelineParams-output-dir /files/data/output --PipelineParams-file-parent-dir /files/data --PipelineParams-date $(DATA_DATE) --PipelineParams-gene-config-path $(GENE_CONFIG_FILE)

force-run-task: ## Running a specific task (deleting its dependencies first)
	docker run $(COMMON_DOCKER_ARGS) $(PIPELINE_IMAGE) bash -c 'cd /opt/brca-exchange/pipeline/luigi && python prune_output_files.py CompileVCFFiles $(RUN_ARGS) $(RUN_ARGS) $(PIPELINE_ARGS)'
	$(RUN_TASK_CMD) $(RUN_ARGS) $(LUIGI_WORKERS) $(GENE_CONFIG_FILE)

clean-files-from: ## Delete files generated by tasks on the task path between the given task and $(LUIGI_ROOT_TASK)
	docker run $(COMMON_DOCKER_ARGS) $(PIPELINE_IMAGE) bash -c 'cd /opt/brca-exchange/pipeline/luigi && python prune_output_files.py CompileVCFFiles $(LUIGI_ROOT_TASK) $(RUN_ARGS) $(PIPELINE_ARGS)'

show-luigi-graph: ## print some representation of the luigi compute graph on the console
	docker run $(COMMON_DOCKER_ARGS) $(PIPELINE_IMAGE) bash -c 'cd /opt/brca-exchange/pipeline/luigi && PYTHONPATH="/opt/brca-exchange/pipeline:/opt/brca-exchange/pipeline/luigi" luigi-deps-tree --module CompileVCFFiles $(PIPELINE_ARGS) $(LUIGI_ROOT_TASK) | sed -E "s/([^{])\{.*\}([^}])/\1/g" '


test: ## Running pipeline unit tests
	docker run $(COMMON_DOCKER_ARGS) $(PIPELINE_IMAGE) bash -c 'cd /opt/brca-exchange/pipeline/data && bash getdata && cd /opt/brca-exchange/pipeline && pytest --ignore=splicing/'

test-coverage: ## Running pipeline unit tests with coverage information
	docker run $(COMMON_DOCKER_ARGS) $(PIPELINE_IMAGE) bash -c 'cd /opt/brca-exchange/pipeline/data && bash getdata && cd /opt/brca-exchange/pipeline && pytest --cov --ignore=splicing/ && coverage html --include="/opt/brca-exchange/pipeline/*" --omit="*/test_*"'


build-release: start-local-uta checkout build-docker setup-files setup-lovd download-resources download-seqrepo run-pipeline variants-by-source ## create new data release

variants-by-source: ## postprocessing: compute statistics for changes with respect to the last release
	docker run $(COMMON_DOCKER_ARGS) $(PIPELINE_IMAGE) python /opt/brca-exchange/pipeline/utilities/variantsBySource.py  -i /files/data/output/release/built_with_change_types.tsv -c true

prune-release-notes-output: ## postprocessing: removes some files, s.t. luigi tasks can be triggered to include updated release notes into release archive
	make clean-files-from GenerateReleaseNotes

include-release-notes: prune-release-notes-output run-pipeline ## postprocessing: include updated release notes into archive

cleanup-failed: ## postprocessing: cleaning up some dangling FAILED files, which may be confusing
	find $(OUT_DIR) -name 'FAILED*.tsv' -exec rm -i {} \;

.ONESHELL:
tag-release: ## postprocessing: tags and pushes git commit used for the data release
	cd $(CODE_BASE)
	git tag -m $(RELEASE_TAG) -a $(RELEASE_TAG) `git rev-parse HEAD`
	git push origin $(RELEASE_TAG)

.ONESHELL:
push-docker: ## postprocessing: pushes docker image to docker hub
	docker login
	docker push brcachallenge/brca-exchange-pipeline:$(RELEASE_TAG)
	docker logout

.ONEHSELL:
copy-to-previous: ## postprocessing: copy new release tar to previous release folder
# 	take latest, in case there should be several tar.gz from some reason
	cp `ls -t $(OUT_DIR)/release-*.tar.gz | head -n 1` $(PREVIOUS_RELEASE_DIR)
	rm $(PREVIOUS_RELEASE_DIR)/latest_release.tar.gz
	ln -s `ls -t $(PREVIOUS_RELEASE_DIR)/release-*.tar.gz | head -n 1` $(PREVIOUS_RELEASE_DIR)/latest_release.tar.gz

post-release-cmds: cleanup-failed include-release-notes push-docker tag-release copy-to-previous ## postprocessing: do all the necessary postprocessing to wrap up a data release

setup-data-from-latest-release-tar: setup-files ## sets up brca output dir with data contained in release archive from last release (only data from variant merging onwards)
	tar -C $(OUT_DIR) -zxf $(PREVIOUS_RELEASE_PATH)

setup-dev-env: start-local-uta build-docker setup-files download-resources ## setup development environment

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
