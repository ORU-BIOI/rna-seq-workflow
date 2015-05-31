configfile: "config.json"

import os

# Knockout and Wilde type
KO_DS = config["samplegroups"]["1"]
WT_DS = config["samplegroups"]["2"]

# It is possible to give different samples as environmental variables instead
# of using the KO_DS and WT_DS samples. Can be useful for testing with small
# dataset. The useage for providing samples:
# $ SAMPLES="S1 S2 S3" snakemake ...
SAMPLES = os.environ.get("SAMPLES"," ".join([KO_DS,WT_DS])).split()

# paths to Transcriptome annotation, Whole Genome Sequence, WGS bowtie2 index
ANNOTATION   = config["annotation"]
WGS          = config["wgs"]
WGS_BOWTIE2  = config["wgs_bowtie2"]

# Path to tophat transcriptome index, all tophat runs generate this index so 
# it is nice to create the index just once.
TRANS_INDEX  = config["trans_index"]
# Python2 path
## if env variable PY2 exists, this is used, otherwise the 2nd value:
PY2 = os.environ.get("PY2", config["py2"])
# Other tools
PICARD = os.environ.get("PICARD",config["picard"])
RNASEQC = os.environ.get("RNASEQC",config["rnaseqc"])

## function to create a dir if it does not exist (even dir hierarchy)
def mkdir(path_to_dir):
    """Create a directory if it does not exist, otherwise do nothing"""
    if not os.path.isdir(path_to_dir):
        os.makedirs(path_to_dir)

# Rule to do qc of original reads, qc of the cutadapt reads, and
# to generate merged counts of reads mapping to annotated transcripts.
# It trims the reads with cutadapt, runs tophat with annotation to map
# reads to transcripts and uses htseq to generate counts of reads to
# transcripts.
rule rnaseqc_tophat_cutadapt:
    input: 
        expand("./{sample}_fastqc.html",sample = SAMPLES),
        expand("cutadapt/{sample}_fastqc.html",sample = SAMPLES),
        expand("cutadapt/tophat/{sample}/accepted_hits_rg_reordered.bam.bai",sample = SAMPLES),
        "cutadapt/tophat/rnaseqc/samples.info.txt"

# Rule to do qc of original reads, qc of the cutadapt reads, and
# to generate merged counts of reads mapping to annotated transcripts.
# It trims the reads with cutadapt, runs tophat with annotation to map
# reads to transcripts and uses htseq to generate counts of reads to
# transcripts.
rule htseq_tophat_cutadapt:
    input: 
        expand("./{sample}_fastqc.html",sample = SAMPLES),
        expand("cutadapt/{sample}_fastqc.html",sample = SAMPLES),
        expand("cutadapt/tophat/{sample}/accepted_hits.bam.bai",sample = SAMPLES),
      #  expand("cutadapt/tophat/{sample}/accepted_hits.XS.bam.bai",sample = SAMPLES),
      #  expand("cutadapt/tophat/{sample}/accepted_hits.no_XS.bam.bai",sample = SAMPLES),
        "cutadapt/tophat/htseq/map_count.txt"
    output:
        "all-finished.txt"
    shell:"""
	  sleep 6
	  touch {output}
	  sleep 6
          """

# Rule to do qc of original reads, qc of the cutadapt reads, and
# to generate merged counts of reads mapping to annotated transcripts.
# It trims the reads with cutadapt, runs tophat with annotation to map
# reads to transcripts onlly and uses htseq to generate counts of reads to
# transcripts.
rule htseq_tophat_transcriptome_only_cutadapt:
    input: 
        expand("./{sample}_fastqc.html",sample = SAMPLES),
        expand("cutadapt/{sample}_fastqc.html",sample = SAMPLES),
        expand("cutadapt/tophat_transcriptome_only/{sample}/accepted_hits.bam.bai",sample = SAMPLES),
        expand("cutadapt/tophat_transcriptome_only/{sample}/accepted_hits.XS.bam.bai",sample = SAMPLES),
        expand("cutadapt/tophat_transcriptome_only/{sample}/accepted_hits.no_XS.bam.bai",sample = SAMPLES),
        "cutadapt/tophat_transcriptome_only/htseq/map_count.txt"

rule qc:
    input: 
        "{sample}.fastq"
    output: 
        qcout  = "{sample}_fastqc.html",
        zipout = "{sample}_fastqc.html"
    params:
        modules   = "module load bioinfo-tools FastQC",
        threads   = config["qc"]["threads"],
        noextract = "--noextract"
    run:
        abspath = os.path.abspath(output.qcout)
        outdir  = os.path.dirname(abspath)
        mkdir(outdir)
        shell("""
		 {params.modules}
		 fastqc --outdir {outdir} {params.noextract} --threads {params.threads} {input}
	      """)

rule cutadapt:
    input: 
        "{sample}.fastq"
    output:
        cutout = "cutadapt/{sample}.fastq"
    params:
        cutadatp_exec = PY2 + "cutadapt",
        cut           = "5",
        qual_cutoff   = "20",
        adapter       = "AGATCGGAAGAGC",
        min_len       = "20",
    run:
        mkdir(os.path.dirname(output.cutout))
        shell("""{params.cutadatp_exec} --adapter {params.adapter} --minimum-length {params.min_len} --cut {params.cut} --quality-cutoff {params.qual_cutoff} --output {output.cutout} {input}""")

## tophat is python2-code, hence calling with specific py2
rule tophat_index:
    output:
        od = "{transcript}.gff"
    params:
        py2     = PY2 + "python2.7",
        modules = "module load bioinfo-tools bowtie2/2.2.3 samtools tophat/2.0.12",
        path    = "/sw/apps/bioinfo/tophat/2.0.12/milou/bin/tophat",
        ref     = WGS_BOWTIE2, 
        ann     = ANNOTATION
    run:
        outdir = os.path.dirname(output.od)
        mkdir(outdir)

        shell("""
		 {params.modules}
		 {params.py2} {params.path} -G {params.ann} --transcriptome-index={wildcards.transcript} -o {outdir} {params.ref}
	      """)
        
rule tophat:
    input:
        WGS_BOWTIE2+".1.bt2", 
        fastq = "{dir}/{sample}.fastq", 
        ann   = ANNOTATION, 
        ti    = TRANS_INDEX+".gff"
    output:
        topout = "{dir}/tophat/{sample}/accepted_hits.bam"
    params:
        py2     = PY2 + "python2.7",
        modules = "module load bioinfo-tools bowtie2/2.2.3 samtools tophat/2.0.12",
        path    = "/sw/apps/bioinfo/tophat/2.0.12/milou/bin/tophat",
        p       = "16"
    run:
        outdir = os.path.dirname(output.topout)
        mkdir(outdir)
        shell("""
                 {params.modules}
                 {params.py2} {params.path} -p {params.p} -o {outdir} --transcriptome-index={TRANS_INDEX} {WGS_BOWTIE2} {input.fastq}
              """)

# tophat for transcriptome only:
rule tophat_to:
    input:
        WGS_BOWTIE2+".1.bt2", 
        fastq = "{dir}/{sample}.fastq", 
        ann   = ANNOTATION, 
        ti    = TRANS_INDEX+".gff"
    output:
        topout = "{dir}/tophat_transcriptome_only/{sample}/accepted_hits.bam"
    params:
        py2     = PY2 + "python2.7",
        modules = "module load bioinfo-tools bowtie2/2.2.3 samtools tophat/2.0.12",
        path    = "/sw/apps/bioinfo/tophat/2.0.12/milou/bin/tophat",
        p       = "16"
    run:
        outdir = os.path.dirname(output.topout)
        mkdir(outdir)
        shell("""
		 {params.modules}
		 {params.py2} {params.path} -p {params.p} -o {outdir} --transcriptome-index={TRANS_INDEX} -T {WGS_BOWTIE2} {input.fastq}
              """)

rule bam_index:
    input:
        bam = "{sample}.bam"
    output:
        bai = "{sample}.bam.bai"
    params:
        modules = "module load bioinfo-tools samtools",
    shell:
        """
	{params.modules}
        samtools index {input.bam}
        """

rule bam_XS_index:
    input:
        bam = "{sample}.bam"
    output:
        bam = "{sample}.XS.bam"
    params:
        modules = "module load bioinfo-tools samtools",
    shell:
        """
        {params.modules}
	samtools view -h {input.bam} | grep -e "^@" -e "XS:A" | samtools view -bS -o {output.bam} -
        """

rule bam_no_XS_index:
    input:
        bam = "{sample}.bam"
    output:
        bam = "{sample}.no_XS.bam"
    params:
        modules = "module load bioinfo-tools samtools",
    shell:
        """
        {params.modules}
        samtools view -h {input.bam} | grep -ve "XS:A" | samtools view -bS -o {output.bam} -
        """

rule count:
    input:    
        bam = "{dir}/{sample}/accepted_hits.bam",
        ti  = TRANS_INDEX+".gff"
    output:
        countout = "{dir}/htseq/map_count_{sample}.txt"
    params:
        modules = "module load bioinfo-tools samtools",
        path    = PY2 + "htseq-count",
    run:
        outdir = os.path.dirname(output.countout)
        mkdir(outdir)
        shell("""
	    {params.modules}
            echo -e "id\t{wildcards.sample}" > {output.countout}
            samtools view {input.bam} | {params.path} -m union - {input.ti} | cat >> {output.countout}
            """)

rule merge_count:
    input:    
        bam = expand("{{dir}}/htseq/map_count_{sample}.txt",sample=SAMPLES)
    output:
        countout = "{dir}/htseq/map_count.txt"
    shell:
        """
        paste {input.bam} | awk '{{for (i=1; i<=NF; i++) if (i == 1 || i % 2 == 0) printf $i "\t"; print""}}' | cat > {output.countout}
        """

rule read_group:
    input:    
        bam = "{dir}/{sample}/accepted_hits.bam",
    output:
        bamout = "{dir}/{sample}/accepted_hits_rg.bam"
    params:
        picard  = PICARD,
    shell:
        """
        java -jar {params.picard}/AddOrReplaceReadGroups.jar I={input.bam} O={output.bamout} LB=lb PL=illumina PU=pu SM={wildcards.sample}
        """

rule reorder:
    input:    
        bam = "{dir}/{sample}/accepted_hits_rg.bam",
    output:
        bamout = "{dir}/{sample}/accepted_hits_rg_reordered.bam"
    params:
        picard  = PICARD,
        wgs     = WGS,
    shell:
        """
        java -jar {params.picard}/ReorderSam.jar I={input.bam} O={output.bamout} R={params.wgs}
        """

rule rnaseqc:
    input:    
        bam = expand("{{dir}}/{sample}/accepted_hits_rg_reordered.bam",sample=SAMPLES),
    output:
        rnaseqc = "{dir}/rnaseqc",
        info = "{dir}/rnaseqc/samples.info.txt"
    params:
        rnaseqc = RNASEQC,
        wgs     = WGS,
        gtf     = ANNOTATION,
    run:
        mkdir(output.rnaseqc)
        samples = ["\t".join(["ID","Input File","Description"])]
        for i in input.bam:
            s = i.split("/")
            samples.append("\t".join([s[-2],i,"No description"]))
        with open(output.info,"w") as fh:
            fh.write("\n".join(samples))

        shell("""
              java -jar {params.rnaseqc} -o {output.rnaseqc} -r {params.wgs} -s {output.info} -singleEnd -t {params.gtf}
              """)

