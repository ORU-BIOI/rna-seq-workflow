import os

# Knockout and Wilde type
#KO_DS = "AAGACG_read1 AAGGGA_read1 GGACCC_read1".split()
#WT_DS = "CCTTCA_read1 TGTTGC_read1 TTCAGC_read1".split()
#SAMPLES = KO_DS + WT_DS

SAMPLES = ["AAGGGA_read1_head_1G"]

# paths to Transcriptome annotation, Whole Genome Sequence, WGS bowtie2 index
ANNOTATION   = "../../pipeline-Hoefig-Aug2014/genome-mouse-NCBI/Mus_musculus/NCBI/build37.2/Annotation/Genes/genes.gtf"
WGS          = "../../pipeline-Hoefig-Aug2014/genome-mouse-NCBI/Mus_musculus/NCBI/build37.2/Sequence/WholeGenomeFasta/genome.fa"
WGS_BOWTIE2  = "../../pipeline-Hoefig-Aug2014/genome-mouse-NCBI/Mus_musculus/NCBI/build37.2/Sequence/Bowtie2Index/genome"

# Path to tophat transcriptome index
TRANS_INDEX  = "../transcriptome-index/mus_known"
# Python2 path
PY2 = "/pica/h1/brynjar/miniconda/envs/sci2/bin/python2.7"

MODULES={"htseq":["bioinfo-tools",
                  "samtools",
                  "htseq/0.6.1"]}

def mkdir(path_to_dir):
    if not os.path.isdir(path_to_dir):
        os.makedirs(path_to_dir)




rule count_tophat_cutadapt:
    input: expand("./{sample}_fastqc.html",sample=SAMPLES),
           expand("cutadapt/{sample}_fastqc.html",sample=SAMPLES),
           "cutadapt/tophat/edger/map_count.txt"

rule qc:
    message: "Quality control, creates "
    input: "{sample}.fastq"
    output: qcout="{sample}_fastqc.html",zipout="{sample}_fastqc.html"
    params:
        modules="module load bioinfo-tools FastQC",
        threads="16",
        noextract="--noextract"
    run:
        abspath = os.path.abspath(output.qcout)
        outdir = os.path.dirname(abspath)
        mkdir(outdir)
        shell("""fastqc --outdir {outdir} {params.noextract} --threads {params.threads} {input}""")

rule cutadapt:
    input: 
        "{sample}.fastq"
    output:
        cutout="cutadapt/{sample}.fastq"
    params:
        cutadatp_exec="/pica/h1/brynjar/miniconda/envs/sci2/bin/cutadapt",
        cut="5",
        qual_cutoff="20",
        adapter="AGATCGGAAGAGC",
        min_len="20",
    run:
        mkdir(os.path.dirname(output.cutout))
        shell("""{params.cutadatp_exec} --adapter {params.adapter} --minimum-length {params.min_len} --cut {params.cut} --quality-cutoff {params.qual_cutoff} --output {output.cutout} {input}""")

rule tophat_index:
    output:
        od="{transcript}.gff"
    params:
        py2=PY2,
        modules="module load bioinfo-tools bowtie2/2.2.3 samtools tophat/2.0.12",
        path="/sw/apps/bioinfo/tophat/2.0.12/milou/bin/tophat",
        ref=WGS_BOWTIE2, 
        ann=ANNOTATION
    run:
        outdir=os.path.dirname(output.od)
        mkdir(outdir)

        shell("""{params.py2} {params.path} -G {params.ann} --transcriptome-index={wildcards.transcript} -o {outdir} {params.ref}""")
        
rule tophat:
    input:
        WGS_BOWTIE2+".1.bt2", 
        fastq = "{dir}/{sample}.fastq", 
        ann   = ANNOTATION, 
        ti    = TRANS_INDEX+".gff"
    output:
        topout="{dir}/tophat/{sample}/accepted_hits.bam"
    params:
        py2=PY2,
        modules="module load bioinfo-tools bowtie2/2.2.3 samtools tophat/2.0.12",
        path="/sw/apps/bioinfo/tophat/2.0.12/milou/bin/tophat",
        p="16"
    run:
        outdir=os.path.dirname(output.topout)
        mkdir(outdir)
        shell("""{params.py2} {params.path} -p {params.p} -o {outdir} --transcriptome-index={TRANS_INDEX} {WGS_BOWTIE2} {input.fastq}""")

rule count:
    input:    
        bam    = "{dir}/{sample}/accepted_hits.bam",
        ti     = TRANS_INDEX+".gff"
    output:
        countout   = "{dir}/edger/map_count_{sample}.txt"
    params:
        py2=PY2,
        modules="module load bioinfo-tools samtools htseq/0.6.1",
        path="/pica/h1/brynjar/miniconda/envs/sci2/bin/htseq-count",
    run:
        outdir=os.path.dirname(output.countout)
        mkdir(outdir)
        shell("""samtools view {input.bam} | {params.path} -m union - {input.ti} | cat > {output.countout}""")

rule merge_count:
    input:    
        bam      = expand("{{dir}}/edger/map_count_{sample}.txt",sample=SAMPLES)
    output:
        countout = "{dir}/edger/map_count.txt"
    run:
        shell("""echo {input.bam} {output.countout}""")

