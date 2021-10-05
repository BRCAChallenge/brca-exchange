# Variant Scoring

Prototype implementation of the variant scoring algorithm

## Running the Pipeline

1. download the necessary data from gnomad: `download_data.sh`
2. meanwhile, create docker image by changing to 'docker' directory and running `./build_docker.sh`
3. adapt `run_pipeline.sh` to your needs and run it.

The script `run_pipeline.sh` runs two steps: the first (`process_gnomad_data.py`) extracts all the relevant data from the downloaded gnomad data
(that is, currently only the BRCA1 and BRCA2 variants) and saves it to parquet. This step is fairly resource intensive (45min using 16 cores), the
more cores available, the faster. It also needs internet access, as some code artifacts are downloaded on the fly.

The second step (`variant_scoring.py`) processes the relevant data to generate a dataframe with the evidence code. This requires little resources.

Note, that these two scripts don't necessarily need to be run within docker, a conda environment built with `environment.yml` should also do.



