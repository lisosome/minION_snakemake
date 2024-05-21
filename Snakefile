## Snakemake minION workflow - snakemake file
import os

raw_data = config["paths"]["rawdata"] # path of the folder containing the data transferred from the minION
base_out = config["paths"]["baseout"] # path of the folder containing all the steps of the analysis
sample = config["names"]["sample"]
guppy_config = config["guppy"]["config"]
ref = config["reference_genome"]
#aligner = config["paths"]["ngmlr"]
#aligner = config["paths"]["minimap"]
#cov = config["paths"]["mosdepth"]

# Steps will be divided in:
#   0 <- RAW.DATA
#   1 <- BASECALLING
#   2 <- TRIMMING
#   3 <- ALIGNMENT
#   4 <- VARIANT_CALLING


###################
### START RULES ###
###################
### rule all: define all the outputs of the workflow
rule all:        
    input:
        f"""{base_out}/4.VARIANT_CALLING/{sample}.no_qc.vcf.gz""",
        f"""{base_out}/3.ALIGNMENT/coverage/{sample}.mosdepth.global.dist.txt""",
        f"""{base_out}/3.ALIGNMENT/coverage/{sample}.mosdepth.summary.txt"""

### Step 01: basecalling
# First let's link all the pod5 files to a single directory, to make guppy work
rule linking:
    input:
        raw_data
    output:
        directory(base_out + "/0.RAW.DATA")
    resources:
        partition="THIN",
        runtime=180,
        mem_mb=2000,
        cpus_per_task=2
    message: """ Linking fast5/pod5 files! """
    log:
        config["paths"]["log_dir"] + "/linking_for_baselcalling.log",
        config["paths"]["log_dir"] + "/linking_for_baselcalling.err"
    shell:
        """
        mkdir -p {output}

        if [[ -d {input}/fast5_pass ]];then
            echo "Fast5 files detected! Linking to {output}"
            find {input} -type f -name "*.fast5" -exec ln -s {{}} {output}/ \; 2> {log[1]} 1> {log[0]}
        else
            echo "Pod5 files detected! Linking to {output}"
            find {input} -type f -name "*.pod5" -exec ln -s {{}} {output}/ \; 2> {log[1]} 1> {log[0]}
        fi
        """

rule guppy_basecall:
    input:
        fol = rules.linking.output
    output:
        fst = directory(base_out + "/1.BASECALLING")
    params:
        #guppy_args = "--compress_fastq --recursive --gpu_runners_per_device 48 --chunks_per_caller 2500 --chunks_per_runner 3000 --chunk_size 300",
        guppy_args = "--recursive --num_callers 12 --gpu_runners_per_device 4 --chunks_per_runner 512 --chunk_size 3000",
        config = {guppy_config}
    envmodules:
        "ont-guppy-gpu/6.5.7"
    resources:
        partition="GPU",
        runtime=1080,
        mem_mb=192000,
        cpus_per_task=48
    benchmark:
        base_out + "/guppy_benchmark.tsv"
    message: """ Guppy Basecalling """
    log:
        config["paths"]["log_dir"] + "/guppy_basecalling.log",
        config["paths"]["log_dir"] + "/guppy_basecalling.err"
    shell: 
        """
            guppy_basecaller -i {input.fol} -s {output.fst} -c {params.config} -x 'auto' {params.guppy_args} 2> {log[1]} 1> {log[0]}
        """

rule preprocessing:
    input:
        base_out + "/1.BASECALLING"
    output:
        base_out + "/2.TRIMMING/" + sample + ".trimmed_and_clean.fastq"
    resources:
        partition="THIN",
        runtime=5760,
        mem_mb=150000,
        cpus_per_task=5
    conda:
        "envs/long_read_processing_snakemake.yaml"
    params:
        filt_args = "-q 10 -l 500 --headcrop 50"
    message: """ Preprocessing """
    log:
        config["paths"]["log_dir"] + "/preprocessing.log",
        config["paths"]["log_dir"] + "/preprocessing.err"
    shell:
        """
        cat $(ls -rt {input}/pass/*.fastq) | NanoLyse | NanoFilt {params.filt_args} > {output[0]} 
        find {base_out} -type f -name "Nanolyse.log" -exec rm {{}} \; 2>> {log[1]} 1>> {log[0]}
        """

# rule alignment:
#     input:
#         rules.preprocessing.output
#     output:
#         #sam = temp(base_out + "/3.ALIGNMENT/" + sample + ".ngmlr.sam"),
#         bam = base_out + "/3.ALIGNMENT/" + sample + ".ngmlr.bam",
#         bai = base_out + "/3.ALIGNMENT/" + sample + ".ngmlr.bai"
#     resources:
#         partition="EPYC",
#         runtime=5760,
#         mem_mb=200000,
#         cpus_per_task=50
#     conda: 
#         "envs/ngmlr2.yaml"
#     params:
#         ngmlr = aligner,
#         ngmlr_args=f"""-t 50 -r {ref} -x ont"""
#     message: """ ngmlr Alignment """
#     log:
#         config["paths"]["log_dir"] + "/alignment.log",
#         config["paths"]["log_dir"] + "/alignment.err"
#     benchmark:
#         base_out + "/ngmlr_alignment.benchmark.tsv"
#     shell:
#         """
#         {params.ngmlr} {params.ngmlr_args} -q {input} | samtools view -b /dev/stdin | samtools sort -@ 50 -o {output.bam} 2> {log[1]} 1> {log[0]} &&
#         samtools index -@ 50 -b {output.bai} {output.bam} 2>> {log[1]} 1>> {log[0]}
#         """

rule alignment:
    input:
        rules.preprocessing.output
    output:
        bam = base_out + "/3.ALIGNMENT/" + sample + ".minimap.bam",
        bai = base_out + "/3.ALIGNMENT/" + sample + ".minimap.bai"
    resources:
        partition="EPYC",
        runtime=5760,
        mem_mb=150000,
        cpus_per_task=15
    params:
        minimap = os.path.join(workflow.basedir, "software/minimap2"),
        minimap_args=f"""-y -t 10 -ax map-ont {ref}"""
    message: """ minimap2 Alignment """
    log:
        config["paths"]["log_dir"] + "/alignment.log",
        config["paths"]["log_dir"] + "/alignment.err"
    benchmark:
        base_out + "/minimap2_alignment.benchmark.tsv"
    envmodules:
        "samtools/1.17"
    shell:
        """
        {params.minimap} {params.minimap_args} {input} | samtools view -@ 15 -b | samtools sort -@ 15 --write-index -o {output.bam}\#\#idx\#\#{output.bai} 2> {log[1]} 1> {log[0]}
        """


rule coverage:
    input:
        base_out + "/3.ALIGNMENT/" + sample + ".minimap.bam"
    output:
        base_out + "/3.ALIGNMENT/coverage/" + sample + ".mosdepth.global.dist.txt",
        base_out + "/3.ALIGNMENT/coverage/" + sample + ".mosdepth.summary.txt"
    resources:
        partition="THIN",
        runtime=5760,
        mem_mb=40000,
        cpus_per_task=4
    params:
        #mosdepth = cov,
        mos_opt = "mosdepth -x -n -t 4",
        out_dir = f"""{base_out}/3.ALIGNMENT/coverage/{sample}"""
    conda:
        "envs/mosdepth.yaml"
    message: """ Coverage Calculation! """
    log:
        config["paths"]["log_dir"] + "/coverage_calculation.log",
        config["paths"]["log_dir"] + "/coverage_calculation.err"
    shell:
        """ 
        {params.mos_opt} {params.out_dir} {input} 2> {log[1]} 1> {log[0]}
        """


rule sv_calling:
    input:
        base_out + "/3.ALIGNMENT/" + sample + ".minimap.bam"
    output:
        vcf=base_out + "/4.VARIANT_CALLING/" + sample + ".no_qc.vcf.gz"
    resources:
        partition="THIN",
        runtime=5760,
        mem_mb=40000,
        cpus_per_task=8
    conda:
        "envs/long_read_processing_snakemake.yaml"
    params:
        tr = os.path.join(workflow.basedir, "software/sniffles2/human_GRCh38_no_alt_analysis_set.trf.bed")
    message: """ sniffles2 SV calling """
    log:
        config["paths"]["log_dir"] + "/SV_calling.log",
        config["paths"]["log_dir"] + "/SV_calling.err"
    shell:
        """
        sniffles -t 8 --no-qc --tandem-repeats {params.tr} --input {input} --vcf {output.vcf} 2> {log[1]} 1> {log[0]}
        """
