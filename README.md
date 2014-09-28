# NGBI RNASeq pipeline in Snakemake #

Snakemake requires python3, which is usually not the default python interpreter  in Linux distros. So to get a working python3 environment, and a python2 environment, with all packages that we require installed, one can execute the [INSTALL.sh](https://bitbucket.org/binnisb/ngbi-rna-pipeline/src/master/INSTALL.sh) script (for new installs, change the prefix in the file). This downloads and installs miniconda, a virtual environment from ContinuumIO, and sets up two virtual envs, sci2 and sci3 which are python2 and python3 respectively. 

### Execute pipeline on sbatch queue system ###
First we need to activate the base conda environment by adding it to path and then we activate the sci3 environment:
```
#!bash
export PATH="/proj/b2014206/miniconda/bin:$PATH"
source activate sci3
```
Then we also load the required modules on uppmax:
```
#!bash
module load bioinfo-tools FastQC/0.11.2 samtools/0.1.19 bowtie2/2.2.3 tophat/2.0.12
```
To execute the pipeline on uppmax with the sbatch system we can execute the [Snakefile-sbatch.sh](https://bitbucket.org/binnisb/ngbi-rna-pipeline/src/master/Snakefile-sbatch.sh) which in turn executes the python script [Snakefile-sbatch.py](https://bitbucket.org/binnisb/ngbi-rna-pipeline/src/master/Snakefile-sbatch.sh) with sensible parameters. This submits all unfinished jobs to the sbatch queue with dependencies, so jobs that require the output from another job(s) as input do not start until all the inputs are ready.