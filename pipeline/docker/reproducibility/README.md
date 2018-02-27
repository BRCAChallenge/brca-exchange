## BRCA Exchange Reproducibility

Reproducing the results of the entire pipeline is difficult since raw data at the sources may have evolved meanwhile and also because the access to raw data used for some sources is restricted. For the latter case, you'd have to request access to those.

The pipeline consists of 2 stages where the first one is to normalize the data from each source and the second to merge all variants and do further processing. Since the output of the first stage is included in a release archive, the second stage can be relatively easily reproduced. Note, however, that also for this part some discrepancies may arise due to the fact, that the second stage relies on external webservices.

## Reproducing Merging Stage
  
1. Make sure you have docker installed on your system 
2. Copy, edit accordingly the variables annotated with `PLEASE EDIT` and run [this script](https://github.com/BRCAChallenge/brca-exchange/blob/master/pipeline/docker/reproducibility/reproduce_merging.sh): 
   It will download necessary release archives, some auxiliary data and run the docker image of the merging part.
3. If everything went well, the script will output a path were the newly generated release archive can be found.
