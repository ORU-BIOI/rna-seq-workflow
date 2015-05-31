# If Snakefile in same directory, run so:
# $ ./scheduler.sh rule1 rule2
# If snakemake in another directory:
# $ ./scheduler.sh -s ../Snakefile rule1 rule2
# you can even pass in parameters:
# $ ./scheduler.sh -n -s ../Snakefile rule1 rule2
snakemake -j 350000 --debug -p --immediate-submit --cluster './scheduler.py {dependencies}' $*
