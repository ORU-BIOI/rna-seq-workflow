# NGBI RNASeq pipeline in Snakemake #

Snakemake requires python3, which is usually not the default python interpreter  in Linux distros. So to get a working python3 environment, and a python2 environment, with all packages that we require installed, one can execute the [INSTALL.sh](INSTALL.sh) script (for new installs, change the prefix in the file). This downloads and installs miniconda, a virtual environment from ContinuumIO, and sets up two virtual envs, sci2 and sci3 which are python2 and python3 respectively. 

### Execute pipeline ###
