# brca-exchange
Overall management and deployment of the BRCA Exchange web portal and pipeline scripts

The goals of the BRCA Exchange are to aggregate, curate and share BRCA variation data, which generally exists in disparate silos.  

##brca-website
### deployment
Deployment scripts for brca-website are located in the root of this repository. There are two scripts:

`deploy-beta`   deploy to brcaexchange.cloudapp.net, the staging site.

`deploy-prod`   deploy to brcaexchange.org, the production site.

A collection of configuration files for the frontend and backend are located in `site_settings`. These config files have been scrubbed of secrets and passwords, so those fields need to be populated.

The deployment scripts and the site\_settings folder expect to live in `~brca/`, and they expect the code to live in `~brca/brca-exchange/brca-website/`.

To deploy the latest commit to the staging site:
```
cd ~brca/brca-exchange/brca-website
git pull
cd ~brca
./deploy-beta
```

To deploy to the production site:
```
cd ~brca
./deploy-prod
```

##Authors
[Benedict Paten](https://github.com/benedictpaten/), [Charles Markello](https://github.com/cmarkello), [Molly Zhang](https://github.com/MollyZhang), [Max Haeussler](https://github.com/maximilianh), [Melissa Cline](https://github.com/melissacline), [Mark Diekhans](https://github.com/diekhans)
