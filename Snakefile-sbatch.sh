 snakemake -j 99 --debug --immediate-submit --cluster './Snakefile-sbatch.py {dependencies}' count_tophat_cutadapt
