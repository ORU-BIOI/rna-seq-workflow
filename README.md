# NGBI RNASeq pipeline in Snakemake #

Snakemake requires python3, which is usually not the default python interpreter  in Linux distros. So to get a working python3 environment, and a python2 environment, with all packages that we require installed, one can execute the [INSTALL.sh](https://bitbucket.org/binnisb/ngbi-rna-pipeline/src/master/INSTALL.sh) script (for new installs, change the prefix in the file). This downloads and installs miniconda, a virtual environment from ContinuumIO, and sets up two virtual envs, sci2 and sci3 which are python2 and python3 respectively. 

### Execute pipeline ###
First we need to activate the base conda environment by adding it to path and then we activate the sci3 environment:
```
#!bash
export PATH="/proj/b2014206/miniconda/bin:$PATH"
source activate sci3
```
Then we also load the required modules on uppmax (make sure these are installed otherwise):
```
#!bash
module load bioinfo-tools FastQC/0.11.2 samtools/0.1.19 bowtie2/2.2.3 tophat/2.0.12
```
To execute the pipeline on uppmax with the sbatch queueing system we can execute the [Snakefile-sbatch.sh](https://bitbucket.org/binnisb/ngbi-rna-pipeline/src/master/Snakefile-sbatch.sh) which in turn executes the python script [Snakefile-sbatch.py](https://bitbucket.org/binnisb/ngbi-rna-pipeline/src/master/Snakefile-sbatch.py) with sensible parameters. This submits all unfinished jobs to the sbatch queue with dependencies, so jobs that require the output from another job(s) as input do not start until all the inputs are ready.

If we are not on a sbatch queue system we can exectue the pipeline without the previous scripts. then each rule runs when the previous rule finishes.

To visualize the workflow, one can execute the [dag.sh](https://bitbucket.org/binnisb/ngbi-rna-pipeline/src/master/dag.sh) script, which calls snakemake --dag with some parameres:
```
#!bash
# To visualize the workflow executed when we run the rule htseq_tophat_cutadapt 
# we can either run the dag.sh script or call snakemake in the following way
snakemake htseq_tophat_cutadapt --dag | dot -Tpdf > HTSeqTophatCutadapt.pdf 

# To run the pipeline on a regular server, execute the following snakemake rule:
snakemake htseq_tophat_cutadapt

# To run the pipeline on sbatch server, either run the Snakefile_sbatch.sh or
# call snakemake in the following way:
 snakemake -j 99 --debug --immediate-submit --cluster './Snakefile-sbatch.py {dependencies}' htseq_tophat_cutadapt
```