.DEFAULT_GOAL := help

CONFIG_PATH := brca_pipeline_cfg.mk

ifneq ("$(wildcard $(CONFIG_PATH))","")
	include $(CONFIG_PATH)
endif

OS := $(shell uname -s)

MKFILE_PATH := $(abspath $(lastword $(MAKEFILE_LIST)))


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
	[ -d $(SYNAPSE_CACHE) ] || mkdir -p $(SYNAPSE_CACHE)
	[ -d $(RESOURCES_DIR) ] || mkdir -p $(RESOURCES_DIR)
	[ -d $(PRIORS_REFERENCES) ] || mkdir -p $(PRIORS_REFERENCES)
	[ -d $(OUT_DIR) ] || mkdir -p $(OUT_DIR)
	[ -f $(PREVIOUS_RELEASE_DIR)/latest_release.tar.gz ] || wget http://brcaexchange.org/backend/downloads/releases/current_release.tar.gz -O $(PREVIOUS_RELEASE_DIR)/latest_release.tar.gz

.ONESHELL:
setup-lovd: ## setting up LOVD data via pipeline machine, if not running on pipeline machine, as the access is IP restricted. Requires access to main pipeline machine	
	if [ `hostname -s` != "brca-pipeline" ]; then
		mkdir -p $(OUT_DIR)/LOVD
		ssh pipeline@52.151.62.136 curl https://databases.lovd.nl/shared/export/BRCA > $(OUT_DIR)/LOVD/BRCA.txt
	fi


ifeq ($(OS), Darwin)
DOCKER_GRP = 0
else
DOCKER_GRP = `stat -c '%g' /var/run/docker.sock`
endif

DOCKER_RUN_CMD = docker run -it --rm -u `id -u ${USER}`:$(DOCKER_GRP) \
	-e "DATA_DATE=$(DATA_DATE)" \
        -e"UTA_DB_URL=postgresql://anonymous@0.0.0.0:$(UTA_PORT)/uta/uta_$(UTA_RELEASE_DATE)" \
	--network host \
	-v $(RESOURCES_DIR):/files/resources \
	-v $(OUT_DIR):/files/data \
	-v $(CREDENTIALS_PATH):/opt/luigi_pipeline_credentials.cfg \
	-v $(PREVIOUS_RELEASE_PATH):/files/previous_release.tar.gz \
	-v $(RELEASE_NOTES_PATH):/files/release_notes.txt \
	-v ${SYNAPSE_CACHE}:/.synapseCache \
  -v $(CODE_BASE):/opt/brca-exchange \
	-v /var/run/docker.sock:/var/run/docker.sock \
	$(PIPELINE_IMAGE)

run-pipeline: ## running entire pipeline
	TS=`date +%Y%m%d_%H%M%S`
	LOG_FILE=$(LOG_DIR)/pipeline_run_$${TS}.log
	echo "Log files are in $${LOG_FILE}"
	$(DOCKER_RUN_CMD) /opt/brca-exchange/pipeline/docker/run_luigi.sh $(PRIORS_REFERENCES) $(OUT_DIR) $(PRIORS_IMAGE) > $${LOG_FILE} 2>&1

run-interactive: ## starting docker container in interactive mode
	$(DOCKER_RUN_CMD) bash

# If the first argument is "run-task"...
# https://stackoverflow.com/questions/2214575/passing-arguments-to-make-run
ifeq (run-task,$(firstword $(MAKECMDGOALS)))
  # use the rest as arguments for "run"
  RUN_ARGS := $(wordlist 2,$(words $(MAKECMDGOALS)),$(MAKECMDGOALS))
  # ...and turn them into do-nothing targets
  $(eval $(RUN_ARGS):;@:)
endif

run-task: ## Running a specific task
	$(DOCKER_RUN_CMD) /opt/brca-exchange/pipeline/docker/run_luigi.sh $(PRIORS_REFERENCES) $(OUT_DIR) $(PRIORS_IMAGE) $(RUN_ARGS)

build-release: start-local-uta checkout build-docker setup-files setup-lovd download-resources run-pipeline variants-by-source ## create new data release

variants-by-source: ## postprocessing: compute statistics for changes with respect to the last release
	$(DOCKER_RUN_CMD) python /opt/brca-exchange/pipeline/utilities/variantsBySource.py  -i /files/data/output/release/built_with_change_types.tsv -c true

prune-release-notes-output: ## postprocessing: removes some files, s.t. luigi tasks can be triggered to include updated release notes into release archive
	rm $(OUT_DIR)/release-*.tar.gz
	rm $(OUT_DIR)/output/md5sums.txt
	rm $(OUT_DIR)/output/README.txt
	rm $(OUT_DIR)/output/release/metadata/version.json

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
	docker login hub.docker.com
	docker push brcachallenge/brca-exchange-pipeline:$(RELEASE_TAG)
	docker logout hub.docker.com

.ONEHSELL:
copy-to-previous: ## postprocessing: copy new release tar to previous release folder 
# 	take latest, in case there should be several tar.gz from some reason
	cp `ls -t $(OUT_DIR)/release-*.tar.gz | head -n 1` $(PREVIOUS_RELEASE_DIR)
	rm $(PREVIOUS_RELEASE_DIR)/latest_release.tar.gz
	ln -s `ls -t $(PREVIOUS_RELEASE_DIR)/release-*.tar.gz | head -n 1` $(PREVIOUS_RELEASE_DIR)/latest_release.tar.gz

post-release-cmds: cleanup-failed include-release-notes push-docker tag-release copy-to-previous ## postprocessing: do all the necessary postprocessing to wrap up a data release

setup-data-from-latest-release-tar: setup-files ## sets up brca output dir with data contained in release archive from last release (only data from variant merging onwards)
	tar -C $(OUT_DIR) -zxf $(PREVIOUS_RELEASE_PATH)


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
