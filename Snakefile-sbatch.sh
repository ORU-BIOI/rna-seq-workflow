 snakemake -j 99 --debug --immediate-submit --cluster './Snakefile-sbatch.py {dependencies}' htseq_tophat_cutadapt
