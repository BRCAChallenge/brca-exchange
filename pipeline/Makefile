
CONFIG_PATH := brca_pipeline_cfg.mk

include $(CONFIG_PATH)

checkout:
	[ -d $(CODE_BASE) ] || git clone $(BRCA_GIT_REPO) $(CODE_BASE)
	cd $(CODE_BASE) && git checkout $(GIT_COMMIT)

build-docker:
	cd $(CODE_BASE)/pipeline/docker && ./build_docker_image.sh $(GIT_COMMIT) $(RELEASE_TAG) || { echo "$(date) Something wrong with docker image creation"; exit 1; }

download-resources:
	$(CODE_BASE)/pipeline/download_resources_files.sh $(RESOURCES_DIR) >> "$(LOG_DIR)/download_resources_$(RELEASE_TAG).log" 2>&1 || { echo "Downloading resource files failed"; exit 1; }


start-local-uta:
	[ `docker ps -f name=$(UTA_CONTAINER) | wc -l` -gt 1 ] || docker run -dit --name $(UTA_CONTAINER) -p $(UTA_PORT):5432 $(UTA_DOCKER_IMAGE)

.ONESHELL:
setup-files:
	[ -f $(CREDENTIALS_PATH) ] || { touch $(CREDENTIALS_PATH); echo "WARNING: files $(CREDENTIALS_PATH) doesn't exist.Creating it"; }
	[ -f $(RELEASE_NOTES_PATH) ] || { touch $(RELEASE_NOTES_PATH) ; echo "WARNING: files $(RELEASE_NOTES_PATH) doesn't exist. Creating it"; }
	[ -d $(SYNAPSE_CACHE) ] || mkdir -p $(SYNAPSE_CACHE)
	[ -d $(RESOURCES_DIR) ] || mkdir -p $(RESOURCES_DIR)
	[ -d $(PRIORS_REFERENCES) ] || mkdir -p $(PRIORS_REFERENCES)
	[ -d $(OUT_DIR) ] || mkdir -p $(OUT_DIR)
	[ -f $(PREVIOUS_RELEASE_DIR)/latest_release.tar.gz ] || wget http://brcaexchange.org/backend/downloads/releases/current_release.tar.gz -O $(PREVIOUS_RELEASE_DIR)/latest_release.tar.gz

.ONESHELL:
setup-lovd:
#	setting up LOVD data via pipeline machine, if not running on pipeline machine, as the access is IP restricted
	if [ `hostname -s` != "brca-pipeline" ]; then
		mkdir -p $(OUT_DIR)/brca_out/LOVD
		ssh pipeline@52.151.62.136 curl https://databases.lovd.nl/shared/export/BRCA > $(OUT_DIR)/brca_out/LOVD/BRCA.txt
	fi

# TODO: docker access within docker doesn't work yet on Mac OS X
DOCKER_RUN_CMD = docker run -it --rm -u `id -u ${USER}`:`stat -c '%g' /var/run/docker.sock` \
	-e "DATA_DATE=$(DATA_DATE)" \
        -e"UTA_DB_URL=postgresql://anonymous@0.0.0.0:$(UTA_PORT)/uta/uta_$(UTA_RELEASE_DATE)" \
       --network host \
       -v $(RESOURCES_DIR):/files/resources \
       -v $(OUT_DIR):/files/data \
       -v $(CREDENTIALS_PATH):/opt/luigi_pipeline_credentials.cfg \
       -v $(PREVIOUS_RELEASE_PATH):/files/previous_release.tar.gz \
       -v $(RELEASE_NOTES_PATH):/files/release_notes.txt \
       -v ${SYNAPSE_CACHE}:/.synapseCache \
       -v /var/run/docker.sock:/var/run/docker.sock \
       $(PIPELINE_IMAGE) 

run-pipeline:
	$(DOCKER_RUN_CMD) /opt/run_luigi.sh $(PRIORS_REFERENCES) $(OUT_DIR) $(PRIORS_IMAGE)

run-interactive:
	$(DOCKER_RUN_CMD) bash

variants-by-source:
	$(DOCKER_RUN_CMD) 'python /opt/brca-exchange/pipeline/utilities/variantsBySource.py  -i /files/data/output/release/built_with_change_types.tsv -c true'

build-release: start-local-uta checkout build-docker setup-files setup-lovd download-resources run-pipeline variants-by-source

.ONESHELL:
tag-release:
	cd $(CODE_BASE)
	git tag -m $(RELEASE_TAG) -a $(RELEASE_TAG) `git rev-parse HEAD`
	git push origin $(RELEASE_TAG)

push-docker:
	docker push brcachallenge/brca-exchange-pipeline:$(RELEASE_TAG)

.ONEHSELL:
copy-to-previous:
	cp `ls -t /home/pipeline/monthly_releases/$(RELEASE_TAG)/brca_out/release-*.tar.gz | head -n 1` $(PREVIOUS_RELEASE_DIR)
	rm $(PREVIOUS_RELEASE_DIR)/latest_release.tar.gz
	ln -s `ls -t $(PREVIOUS_RELEASE_DIR)/release-*.tar.gz | head -n 1` $(PREVIOUS_RELEASE_DIR)/latest_release.tar.gz

post-release-cmds: push-docker tag-release copy-to-previous

