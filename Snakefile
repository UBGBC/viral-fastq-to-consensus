#!/usr/bin/env python3
import os
import pandas as pd
from snakemake.shell import shell
shell.executable("/bin/bash")

configfile: "config.json"
localrules: all, mkdir

# ---------- sample table ----------
df = pd.read_csv(config["meta_file"], sep="\t", header=0, index_col=0)
sample_ids = list(df.index)
df.index = sample_ids

def get_pair_gz(sample_id):
    d = config["raw_fastq_gz_dir"]
    return tuple(os.path.join(d, df.loc[str(sample_id), x]) for x in ("ForwardFastqGZ","ReverseFastqGZ"))

def get_forward_primer(s): return df.loc[s]["Adapter_1"]
def get_reverse_primer(s): return df.loc[s]["Adapter_2"]
def get_year(s):           return df.loc[s]["Year"]
def get_date(s):           return df.loc[s]["Date"]

# NEW: reference / index values come from samplesheet
def get_reference_fasta(s):
    fa = df.loc[s]["ReferenceFasta"]
    if not isinstance(fa, str) or fa.strip()=="":
        raise ValueError(f"ReferenceFasta missing for sample {s}")
    return fa

def get_index_prefix(s):
    v = df.loc[s].get("IndexPrefix", "")
    return None if (pd.isna(v) or str(v).strip()=="") else str(v)

def chosen_prefix(s):
    p = get_index_prefix(s)
    return p if p else os.path.splitext(get_reference_fasta(s))[0]

# ---------- aligner + index helpers ----------
def aligner():
    # still controlled globally (can be per-sample later if you want)
    return config["params"].get("aligner","bwa").lower()

def bwa_index_files(prefix):
    return [prefix + ext for ext in [".amb",".ann",".bwt",".pac",".sa"]]

def bowtie2_index_files(prefix):
    return [prefix + ext for ext in [".1.bt2",".2.bt2",".3.bt2",".4.bt2",".rev.1.bt2",".rev.2.bt2"]]

def index_targets(prefix):
    return bowtie2_index_files(prefix) if aligner()=="bowtie2" else bwa_index_files(prefix)

# ---------- DAG ----------
rule all:
    input:
        expand("{dir}/{sample_id}.masked_consensus.fasta",
               dir=config["dir_names"]["consensus_dir"], sample_id=sample_ids),
        config["dir_names"]["multiqc_dir"] + "/multiqc_report.html"

rule mkdir:
    output: touch(config["file_names"]["mkdir_done"])
    params: dirs=list(config["dir_names"].values())
    resources: time_min=10, mem_mb=2000, cpus=1
    shell: "mkdir -p {params.dirs}"

# -------------------- QC (pre-trim) --------------------
rule fastqc_pre:
    input:
        rules.mkdir.output,
        r1=lambda wc: get_pair_gz(wc.sample_id)[0],
        r2=lambda wc: get_pair_gz(wc.sample_id)[1]
    output:
        html1 = config["dir_names"]["qc_pre_dir"] + "/{sample_id}_R1_fastqc.html",
        zip1  = config["dir_names"]["qc_pre_dir"] + "/{sample_id}_R1_fastqc.zip",
        html2 = config["dir_names"]["qc_pre_dir"] + "/{sample_id}_R2_fastqc.html",
        zip2  = config["dir_names"]["qc_pre_dir"] + "/{sample_id}_R2_fastqc.zip"
    params: outdir=config["dir_names"]["qc_pre_dir"]
    threads: 2
    shell: "fastqc -t {threads} -o {params.outdir} {input.r1} {input.r2}"

rule fastq_screen_pre:
    input:
        r1=lambda wc: get_pair_gz(wc.sample_id)[0],
        r2=lambda wc: get_pair_gz(wc.sample_id)[1]
    output:
        png1 = config["dir_names"]["qc_pre_dir"] + "/{sample_id}_R1_screen.png",
        png2 = config["dir_names"]["qc_pre_dir"] + "/{sample_id}_R2_screen.png"
    params:
        conf=config["params"]["fastq_screen"]["conf"],
        outdir=config["dir_names"]["qc_pre_dir"]
    threads: 2
    shell: r"""
        fastq_screen --aligner bowtie2 --threads {threads} --conf {params.conf} --outdir {params.outdir} {input.r1} {input.r2}
        [[ -f {params.outdir}/{wildcards.sample_id}_R1_screen.png ]] && ln -sf {params.outdir}/{wildcards.sample_id}_R1_screen.png {output.png1} || touch {output.png1}
        [[ -f {params.outdir}/{wildcards.sample_id}_R2_screen.png ]] && ln -sf {params.outdir}/{wildcards.sample_id}_R2_screen.png {output.png2} || touch {output.png2}
    """

# -------------------- Trimming --------------------
rule trim:
    input:
        rules.fastqc_pre.output,
        rules.fastq_screen_pre.output,
        all_read1=lambda wc: get_pair_gz(wc.sample_id)[0],
        all_read2=lambda wc: get_pair_gz(wc.sample_id)[1]
    output:
        trimmed_read1=config["dir_names"]["trimmed_dir"] + "/{sample_id}.trimmed.R1.fastq.gz",
        trimmed_read2=config["dir_names"]["trimmed_dir"] + "/{sample_id}.trimmed.R2.fastq.gz",
        trimmed_stats=config["dir_names"]["trimmed_dir"] + "/{sample_id}.trimmed.stats"
    resources: time_min=360, mem_mb=2000, cpus=4
    version: config["tool_version"]["cutadapt"]
    params:
        adapter1=lambda wc: get_forward_primer(wc.sample_id),
        adapter2=lambda wc: get_reverse_primer(wc.sample_id)
    shell: r"""
        cutadapt -j 4 -m 15 -a {params.adapter1} -A {params.adapter2} -n 2 \
          -o {output.trimmed_read1} -p {output.trimmed_read2} {input.all_read1} {input.all_read2} > {output.trimmed_stats}
    """

# -------------------- QC (post-trim) --------------------
rule fastqc_post:
    input:
        r1=rules.trim.output.trimmed_read1,
        r2=rules.trim.output.trimmed_read2
    output:
        html1=config["dir_names"]["qc_post_dir"] + "/{sample_id}.trimmed_R1_fastqc.html",
        zip1 =config["dir_names"]["qc_post_dir"] + "/{sample_id}.trimmed_R1_fastqc.zip",
        html2=config["dir_names"]["qc_post_dir"] + "/{sample_id}.trimmed_R2_fastqc.html",
        zip2 =config["dir_names"]["qc_post_dir"] + "/{sample_id}.trimmed_R2_fastqc.zip"
    params: outdir=config["dir_names"]["qc_post_dir"]
    threads: 2
    shell: "fastqc -t {threads} -o {params.outdir} {input.r1} {input.r2}"

rule fastq_screen_post:
    input:
        r1=rules.trim.output.trimmed_read1,
        r2=rules.trim.output.trimmed_read2
    output:
        png1=config["dir_names"]["qc_post_dir"] + "/{sample_id}.trimmed_R1_screen.png",
        png2=config["dir_names"]["qc_post_dir"] + "/{sample_id}.trimmed_R2_screen.png"
    params:
        conf=config["params"]["fastq_screen"]["conf"],
        outdir=config["dir_names"]["qc_post_dir"]
    threads: 2
    shell: r"""
        fastq_screen --aligner bowtie2 --threads {threads} --conf {params.conf} --outdir {params.outdir} {input.r1} {input.r2}
        [[ -f {params.outdir}/{wildcards.sample_id}.trimmed_R1_screen.png ]] && ln -sf {params.outdir}/{wildcards.sample_id}.trimmed_R1_screen.png {output.png1} || touch {output.png1}
        [[ -f {params.outdir}/{wildcards.sample_id}.trimmed_R2_screen.png ]] && ln -sf {params.outdir}/{wildcards.sample_id}.trimmed_R2_screen.png {output.png2} || touch {output.png2}
    """

# -------------------- Host read removal (Kraken2) --------------------
rule kraken2_dehost:
    input:
        r1=rules.trim.output.trimmed_read1,
        r2=rules.trim.output.trimmed_read2
    output:
        nh1    = config["dir_names"]["dehost_dir"] + "/{sample_id}.dehost.R1.fastq.gz",
        nh2    = config["dir_names"]["dehost_dir"] + "/{sample_id}.dehost.R2.fastq.gz",
        report = config["dir_names"]["dehost_dir"] + "/{sample_id}.kraken2.report.txt"
    params:
        db=config["params"]["kraken2"]["db"]
    threads: lambda wc: int(config["params"]["kraken2"].get("threads", 8))
    resources: time_min=360, mem_mb=4000, cpus=lambda wc: int(config["params"]["kraken2"].get("threads", 8))
    shell: r"""
        tmpdir=$(mktemp -d)
        kraken2 --db {params.db} --threads {threads} --paired \
            --report {output.report} \
            --unclassified-out $tmpdir/{wildcards.sample_id}_unclass#.fastq \
            {input.r1} {input.r2}
        pigz -c $tmpdir/{wildcards.sample_id}_unclass_1.fastq > {output.nh1}
        pigz -c $tmpdir/{wildcards.sample_id}_unclass_2.fastq > {output.nh2}
        rm -rf "$tmpdir"
    """

# -------------------- Index building (per-sample) --------------------
rule build_bwa_index:
    input:
        fa=lambda wc: get_reference_fasta(wc.sample_id)
    output:
        lambda wc: bwa_index_files(chosen_prefix(wc.sample_id)) if aligner()=="bwa" else [temp("ignore_bwa_index")]
    params:
        prefix=lambda wc: chosen_prefix(wc.sample_id)
    threads: 2
    run:
        if aligner()!="bwa":
            shell("touch {output[0]}")
        else:
            shell("bwa index -p {params.prefix} {input.fa}")

rule build_bowtie2_index:
    input:
        fa=lambda wc: get_reference_fasta(wc.sample_id)
    output:
        lambda wc: bowtie2_index_files(chosen_prefix(wc.sample_id)) if aligner()=="bowtie2" else [temp("ignore_bt2_index")]
    params:
        prefix=lambda wc: chosen_prefix(wc.sample_id)
    threads: 2
    run:
        if aligner()!="bowtie2":
            shell("touch {output[0]}")
        else:
            shell("bowtie2-build {input.fa} {params.prefix}")

rule ensure_index:
    input:
        lambda wc: index_targets(chosen_prefix(wc.sample_id))
    output:
        touch(config["file_names"]["index_ready"])
    shell: "touch {output}"

# -------------------- Mapping (bwa or bowtie2) --------------------
rule map:
    input:
        p1 = rules.kraken2_dehost.output.nh1,
        p2 = rules.kraken2_dehost.output.nh2,
        idx_ready = rules.ensure_index.output
    output:
        mapped_bam_file = config["dir_names"]["mapped_dir"] + "/{sample_id}.bam"
    resources: time_min=360, mem_mb=20000, cpus=lambda wc: int(config["params"].get(aligner(), {}).get("threads", 6))
    params:
        threads=lambda wc: int(config["params"].get(aligner(), {}).get("threads", 6)),
        prefix =lambda wc: chosen_prefix(wc.sample_id)
    run:
        if aligner()=="bowtie2":
            shell("""
                bowtie2 -x {params.prefix} -1 {input.p1} -2 {input.p2} -p {params.threads} \
                | samtools view -F 12 -Sb - \
                | samtools sort -T {output.mapped_bam_file}.align -o {output.mapped_bam_file}
            """)
        else:
            shell("""
                bwa mem -t {params.threads} {params.prefix} {input.p1} {input.p2} \
                | samtools view -F 12 -Sb - \
                | samtools sort -T {output.mapped_bam_file}.align -o {output.mapped_bam_file}
            """)

# -------------------- Downstream steps --------------------
rule ivar_filter:
    input:
        sorted_bam = rules.map.output.mapped_bam_file,
        primer_bed = config["params"]["ivar"]["primer_bed"]
    output:
        filtered_bam_file = config["dir_names"]["filtered_dir"] + "/{sample_id}.filtered.bam"
    resources: time_min=360, mem_mb=20000, cpus=6
    shell: r"""
        #ivar trim -i {input.sorted_bam} -b {input.primer_bed} -p {output.filtered_bam_file}
        ivar trim -e -k -i {input.sorted_bam} -b {input.primer_bed} -p {output.filtered_bam_file}
    """

rule second_sort:
    input:  filtered_bam = rules.ivar_filter.output.filtered_bam_file
    output: second_sorted_bam_file = config["dir_names"]["filtered_dir"] + "/{sample_id}.filtered.sorted.bam"
    resources: time_min=360, mem_mb=20000, cpus=1
    shell:  "samtools view -hb {input.filtered_bam} | samtools sort -T {input.filtered_bam}.tmp -o {output.second_sorted_bam_file}"

rule pileup:
    input:
        second_sorted_bam = rules.second_sort.output.second_sorted_bam_file,
    output:
        pileup = config["dir_names"]["mpileup_dir"] + "/{sample_id}.mpileup"
    params:
        depth = config["params"]["mpileup"]["depth"],
        min_base_qual = config["params"]["varscan"]["snp_qual_threshold"],
        reference=lambda wc: get_reference_fasta(wc.sample_id)
    resources: time_min=360, mem_mb=10000, cpus=1
    shell: r"""
        bcftools mpileup -A -Q {params.min_base_qual} --max-depth 5000000 -L 5000000 \
            -f {params.reference} {input.second_sorted_bam} > {output.pileup}
    """

rule call_snps:
    input:
        pileup = rules.pileup.output.pileup,
        second_sorted_bam = rules.second_sort.output.second_sorted_bam_file
    output:
        vcf = config["dir_names"]["varscan_dir"] + "/{sample_id}.vcf.gz"
    params:
        depth = config["params"]["mpileup"]["depth"],
        min_base_qual = config["params"]["varscan"]["snp_qual_threshold"],
        min_base_cov = config["params"]["varscan"]["min_cov"],
        reference=lambda wc: get_reference_fasta(wc.sample_id)
    resources: time_min=360, mem_mb=10000, cpus=1
    shell: r"""
        bcftools mpileup -A -a "INFO/AD,INFO/ADF,INFO/ADR,FORMAT/ADF,FORMAT/ADR,FORMAT/SP" \
            -Q {params.min_base_qual} -f {params.reference} -L 5000000 --max-depth 5000000 \
            -Ou {input.second_sorted_bam} \
        | bcftools call -Ou -mv \
        | bcftools norm -f {params.reference} -Ou \
        | bcftools filter --include '(TYPE="INDEL" && IMF > .3 && IDV > 30) || (TYPE="SNP" && DP > {params.min_base_cov})' -Oz -o {output.vcf}
    """

rule call_consensus_snps:
    input:  vcf = rules.call_snps.output.vcf
    output: consensus_vcf = config["dir_names"]["varscan_dir"] + "/{sample_id}.consensus.vcf.gz"
    resources: time_min=360, mem_mb=10000, cpus=1
    shell:  "bcftools filter --exclude '(AD[0])/ (AD[0] + AD[1]) >= 0.5' {input.vcf} -Oz -o {output.consensus_vcf}"

rule index_bcf:
    input:
        vcf = rules.call_snps.output.vcf,
        consensus_vcf = rules.call_consensus_snps.output.consensus_vcf
    output:
        index = config["dir_names"]["varscan_dir"]+"/{sample_id}.vcf.gz.tbi",
        consensus_index = config["dir_names"]["varscan_dir"]+"/{sample_id}.consensus.vcf.gz.tbi"
    resources: time_min=360, mem_mb=10000, cpus=1
    shell: "tabix -f {input.vcf} > {output.index}; tabix -f {input.consensus_vcf} > {output.consensus_index};"

rule bedtools_mask:
    input:
        second_sorted_bam = rules.second_sort.output.second_sorted_bam_file,
        vcf = rules.call_snps.output.vcf
    output:
        bed_file = config["dir_names"]["varscan_dir"]+"/{sample_id}.bed"
    resources: time_min=360, mem_mb=10000, cpus=1
    shell: r"""
        bedtools genomecov -ibam {input.second_sorted_bam} -bga \
        | awk '{{if($4 < 50) print $0}}' \
        | bedtools intersect -v -a - -b {input.vcf} > {output.bed_file}
    """

rule mask_consensus:
    input:
        vcf = rules.call_consensus_snps.output.consensus_vcf,
        second_sorted_bam = rules.second_sort.output.second_sorted_bam_file,
        index = rules.index_bcf.output.consensus_index,
        bed_file = rules.bedtools_mask.output.bed_file
    output:
        consensus_genome = config["dir_names"]["consensus_dir"] + "/{sample_id}.masked_consensus.fasta"
    resources: time_min=360, mem_mb=10000, cpus=1
    params:
        reference=lambda wc: get_reference_fasta(wc.sample_id),
        year=lambda wc: get_year(wc.sample_id),
        sample="{sample_id}"
    shell: r"""
        bcftools consensus -f {params.reference} -m {input.bed_file} {input.vcf} > {output.consensus_genome}
    """

# -------------------- MultiQC --------------------
rule multiqc:
    input:
        expand(config["dir_names"]["consensus_dir"]+"/{sample_id}.masked_consensus.fasta", sample_id=sample_ids)
    output:
        report = config["dir_names"]["multiqc_dir"] + "/multiqc_report.html"
    params: 
        outdir=config["dir_names"]["multiqc_dir"],
        multiqc_run_dir=config["dir_names"]["results_dir"]
    threads: 2
    shell: """
    # move to main results directory to gather all reports
    cd {params.multiqc_run_dir}
    multiqc -o {params.outdir} . -T "ViralRecon MultiQC Report {sample_id}"
    """
