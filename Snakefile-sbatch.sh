 snakemake -j 99 --debug --immediate-submit --cluster './Snakefile-sbatch.py {dependencies}' rnaseqc_tophat_cutadapt
