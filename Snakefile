import os

# Knockout and Wilde type
KO_DS = "AAGACG_read1 AAGGGA_read1 GGACCC_read1"
WT_DS = "CCTTCA_read1 TGTTGC_read1 TTCAGC_read1"

# It is possible to give different samples as environmental variables instead
# of using the KO_DS and WT_DS samples. Can be useful for testing with small
# dataset. The useage for providing samples:
# $ SAMPLES="S1 S2 S3" snakemake ...
SAMPLES = os.environ.get("SAMPLES"," ".join([KO_DS,WT_DS])).split()

# paths to Transcriptome annotation, Whole Genome Sequence, WGS bowtie2 index
ANNOTATION   = "../../pipeline-Hoefig-Aug2014/genome-mouse-NCBI/Mus_musculus/NCBI/build37.2/Annotation/Genes/genes.gtf"
WGS          = "../../pipeline-Hoefig-Aug2014/genome-mouse-NCBI/Mus_musculus/NCBI/build37.2/Sequence/WholeGenomeFasta/genome.fa"
WGS_BOWTIE2  = "../../pipeline-Hoefig-Aug2014/genome-mouse-NCBI/Mus_musculus/NCBI/build37.2/Sequence/Bowtie2Index/genome"

# Path to tophat transcriptome index, all tophat runs generate this index so 
# it is nice to create the index just ones.
TRANS_INDEX  = "../transcriptome-index/mus_known"
# Python2 path
PY2 = os.environ.get("PY2","/proj/b2014206/miniconda/envs/sci2/bin/python2.7")
# Other tools
PICARD = os.environ.get("PICARD","/sw/apps/bioinfo/picard/1.118/milou")
RNASEQC = os.environ.get("RNASEQC","/proj/b2014206/opt/RNA-SeQC_v1.1.8.jar")

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
        "cutadapt/tophat/rnaseqc"

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
        expand("cutadapt/tophat/{sample}/accepted_hits.XS.bam.bai",sample = SAMPLES),
        expand("cutadapt/tophat/{sample}/accepted_hits.no_XS.bam.bai",sample = SAMPLES),
        "cutadapt/tophat/htseq/map_count.txt"

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
        threads   = "16",
        noextract = "--noextract"
    run:
        abspath = os.path.abspath(output.qcout)
        outdir  = os.path.dirname(abspath)
        mkdir(outdir)
        shell("""fastqc --outdir {outdir} {params.noextract} --threads {params.threads} {input}""")

rule cutadapt:
    input: 
        "{sample}.fastq"
    output:
        cutout = "cutadapt/{sample}.fastq"
    params:
        cutadatp_exec = "/pica/h1/brynjar/miniconda/envs/sci2/bin/cutadapt",
        cut           = "5",
        qual_cutoff   = "20",
        adapter       = "AGATCGGAAGAGC",
        min_len       = "20",
    run:
        mkdir(os.path.dirname(output.cutout))
        shell("""{params.cutadatp_exec} --adapter {params.adapter} --minimum-length {params.min_len} --cut {params.cut} --quality-cutoff {params.qual_cutoff} --output {output.cutout} {input}""")

rule tophat_index:
    output:
        od = "{transcript}.gff"
    params:
        py2     = PY2,
        modules = "module load bioinfo-tools bowtie2/2.2.3 samtools tophat/2.0.12",
        path    = "/sw/apps/bioinfo/tophat/2.0.12/milou/bin/tophat",
        ref     = WGS_BOWTIE2, 
        ann     = ANNOTATION
    run:
        outdir = os.path.dirname(output.od)
        mkdir(outdir)

        shell("""{params.py2} {params.path} -G {params.ann} --transcriptome-index={wildcards.transcript} -o {outdir} {params.ref}""")
        
rule tophat:
    input:
        WGS_BOWTIE2+".1.bt2", 
        fastq = "{dir}/{sample}.fastq", 
        ann   = ANNOTATION, 
        ti    = TRANS_INDEX+".gff"
    output:
        topout = "{dir}/tophat/{sample}/accepted_hits.bam"
    params:
        py2     = PY2,
        modules = "module load bioinfo-tools bowtie2/2.2.3 samtools tophat/2.0.12",
        path    = "/sw/apps/bioinfo/tophat/2.0.12/milou/bin/tophat",
        p       = "16"
    run:
        outdir = os.path.dirname(output.topout)
        mkdir(outdir)
        shell("""{params.py2} {params.path} -p {params.p} -o {outdir} --transcriptome-index={TRANS_INDEX} {WGS_BOWTIE2} {input.fastq}""")

rule tophat_to:
    input:
        WGS_BOWTIE2+".1.bt2", 
        fastq = "{dir}/{sample}.fastq", 
        ann   = ANNOTATION, 
        ti    = TRANS_INDEX+".gff"
    output:
        topout = "{dir}/tophat_transcriptome_only/{sample}/accepted_hits.bam"
    params:
        py2     = PY2,
        modules = "module load bioinfo-tools bowtie2/2.2.3 samtools tophat/2.0.12",
        path    = "/sw/apps/bioinfo/tophat/2.0.12/milou/bin/tophat",
        p       = "16"
    run:
        outdir = os.path.dirname(output.topout)
        mkdir(outdir)
        shell("""{params.py2} {params.path} -p {params.p} -o {outdir} --transcriptome-index={TRANS_INDEX} -T {WGS_BOWTIE2} {input.fastq}""")

rule bam_index:
    input:
        bam = "{sample}.bam"
    output:
        bai = "{sample}.bam.bai"
    shell:
        """
        samtools index {input.bam}
        """

rule bam_XS_index:
    input:
        bam = "{sample}.bam"
    output:
        bam = "{sample}.XS.bam"
    shell:
        """
        samtools view -h {input.bam} | grep -e "^@" -e "XS:A" | samtools view -bS -o {output.bam} -
        """

rule bam_no_XS_index:
    input:
        bam = "{sample}.bam"
    output:
        bam = "{sample}.no_XS.bam"
    shell:
        """
        samtools view -h {input.bam} | grep -ve "XS:A" | samtools view -bS -o {output.bam} -
        """

rule count:
    input:    
        bam = "{dir}/{sample}/accepted_hits.bam",
        ti  = TRANS_INDEX+".gff"
    output:
        countout = "{dir}/htseq/map_count_{sample}.txt"
    params:
        py2     = PY2,
        modules = "module load bioinfo-tools samtools htseq/0.6.1",
        path    = "/pica/h1/brynjar/miniconda/envs/sci2/bin/htseq-count",
    run:
        outdir = os.path.dirname(output.countout)
        mkdir(outdir)
        shell("""
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
        modules = "module load bioinfo-tools samtools picard/1.118",
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
        modules = "module load bioinfo-tools samtools picard/1.118",
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
        modules = "module load bioinfo-tools samtools picard/1.118",
    run:
        mkdir(output.rnaseqc)
        samples = []
        for i in input.bam:
            s = i.split("/")
            samples.append("\t".join([s[-2],i,"No description"]))
        with open(output.info,"w") as fh:
            fh.write("\n".join(i))

        shell("""
              java -jar {params.rnaseqc} -o {output.rnaseqc} -r {params.wgs} -s {output.info} -singleEnd -t {params.gtf}
              """)

