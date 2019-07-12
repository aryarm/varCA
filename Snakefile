import os
import re
import hashlib
import warnings
from snakemake.utils import min_version

##### set minimum snakemake version #####
min_version("5.5.0")

configfile: "config.yaml"


def read_samples():
    """Function to get names and dna fastq paths from a sample file
    specified in the configuration. Input file is expected to have 3
    columns: <unique_sample_id> <fastq1_path> <fastq2_path>. Modify
    this function as needed to provide a dictionary of sample_id keys and
    (fastq1, fastq2) values"""
    f = open(config['sample_file'], "r")
    samp_dict = {}
    for line in f:
        words = line.strip().split("\t")
        samp_dict[words[0]] = (words[1], words[2])
    return samp_dict
SAMP = read_samples()

# the user can change config['SAMP_NAMES'] here (or define it in the config
# file) to contain whichever sample names they'd like to run the pipeline on
if 'SAMP_NAMES' not in config:
    config['SAMP_NAMES'] = list(SAMP.keys())
else:
    # double check that the user isn't asking for samples they haven't provided
    user_samps = set(config['SAMP_NAMES'])
    config['SAMP_NAMES'] = list(set(SAMP.keys()).intersection(user_samps))
    if len(config['SAMP_NAMES']) != len(user_samps):
        warnings.warn("Not all of the samples requested have provided input. Proceeding with as many samples as is possible...")


def hash_str(to_hash):
    """ hash a str to get a unique value """
    return hashlib.md5(to_hash.encode('utf-8')).hexdigest()[:8]


rule all:
    # if you'd like to run the pipeline on only a subset of the samples,
    # you should specify them in the config['SAMP_NAMES'] variable above
    input:
        expand(
            config['output_dir'] + "/merged_{type}/{sample}.tsv.gz",
            sample=config['SAMP_NAMES'],
            type=[
                i for i in ["snp", "indel"]
                if i+"_callers" in config and config[i+"_callers"]
            ]
        )

rule align:
    """Align reads using BWA-MEM. Note that we use -R to specify read group
    info for haplotype caller"""
    input:
        ref = config['genome'],
        fastq1 = lambda wildcards: SAMP[wildcards.sample][0],
        fastq2 = lambda wildcards: SAMP[wildcards.sample][1]
    output:
        config['output_dir'] + "/align/{sample}/aln.bam"
    threads: config['num_threads']
    conda: "envs/default.yml"
    shell:
        "bwa mem -t {threads} {input.ref} {input.fastq1} {input.fastq2} "
        "-R '@RG\\tID:{wildcards.sample}\\tLB:lib1\\tPL:Illumina\\tPU:unit1\\tSM:{wildcards.sample}' | "
        "samtools view -S -b -h -F 4 -q 20 -> {output}"

rule add_mate_info:
    """Use fixmate to fill in mate coordinates and mate related flags, since
    our data is pair-ended. We need the MC tags (included because we used the
    -m flag) that it creates for markdup"""
    input:
        rules.align.output
    output:
        config['output_dir'] + "/align/{sample}/sorted.mated.bam"
    threads: config['num_threads']
    conda: "envs/default.yml"
    shell:
        "samtools sort -n -@ {threads} {input} | "
        "samtools fixmate -m -@ {threads} -O bam - - | "
        "samtools sort -@ {threads} -o {output} -"

rule rm_dups:
    """Remove duplicates that may have occurred from PCR and index the
    resulting file."""
    input:
        rules.add_mate_info.output
    output:
        final_bam = config['output_dir'] + "/align/{sample}/rmdup.bam",
        final_bam_index = config['output_dir'] + "/align/{sample}/rmdup.bam.bai"
    threads: config['num_threads']
    conda: "envs/default.yml"
    shell:
        "samtools markdup -@ {threads} {input} {output.final_bam} && "
        "samtools index -b -@ {threads} {output.final_bam}"

rule call_peaks:
    """Call peaks in the bam files using macs2"""
    input:
        rules.rm_dups.output.final_bam
    params:
        output_dir = lambda wildcards, output: os.path.dirname(output[0])
    output:
        config['output_dir'] + "/peaks/{sample}/{sample}_peaks.narrowPeak"
    conda: "envs/default.yml"
    shell:
        "macs2 callpeak --nomodel --extsize 200 --slocal 1000 --qvalue 0.05 "
        "-g hs -f BAMPE -t {input} -n {wildcards.sample} --outdir {params.output_dir}"

rule bed_peaks:
    """Convert the BAMPE file to a sorted BED file (with the ref allele)"""
    input:
        ref = config['genome'],
        peaks = rules.call_peaks.output
    output:
        config['output_dir'] + "/peaks/{sample}/peaks.bed"
    conda: "envs/default.yml"
    shell:
        # to convert to BED, we must extract the first three columns (chr, start, stop)
        "cut -f -3 \"{input.peaks}\" | "
        "bedtools getfasta -fi {input.ref} -bedOut -bed stdin | "
        "sort -t $'\t' -k1,1V -k2,2n > \"{output}\""


def get_special_script_path(wildcards):
    """ retrieve the path to the output of a special script if needed """
    if '-' in wildcards.caller:
        special_caller = "-".join(wildcards.caller.split('-')[:-1])
        if os.path.exists("callers/"+special_caller):
            return rules.prepare_caller.output[0].format(
                sample=wildcards.sample,
                caller=special_caller
            )
    return []


rule prepare_caller:
    """Run any scripts that must be run before the caller scripts"""
    input:
        bam = rules.rm_dups.output.final_bam,
        peaks = rules.bed_peaks.output,
        genome = config['genome'],
        shared = get_special_script_path,
        caller_script = "callers/{caller}"
    params:
        caller_params = lambda wildcards: config[wildcards.caller]['params'] if wildcards.caller in config and 'params' in config[wildcards.caller] else ""
    output:
        directory(config['output_dir'] + "/callers/{sample}/{caller}")
    wildcard_constraints:
        sample = "[^\/]*",
        caller = "[^\/]*"
    threads: config['num_threads']
    conda: "envs/default.yml"
    shell:
        "mkdir -p \"{output}\" && "
        "{input.caller_script} {input.bam} {input.peaks} "
        "{input.genome} {output} {wildcards.sample} "
        "{threads} {input.shared} {params.caller_params}"

rule run_caller:
    """Run any callers that are needed"""
    input:
        bam = rules.rm_dups.output.final_bam,
        peaks = rules.bed_peaks.output,
        genome = config['genome'],
        shared = get_special_script_path,
        caller_script = "callers/{caller}"
    params:
        caller_params = lambda wildcards: config[wildcards.caller]['params'] if wildcards.caller in config and 'params' in config[wildcards.caller] else "",
        out_dir = config['output_dir'] + "/callers/{sample}/{caller}"
    output:
        vcf = config['output_dir'] + "/callers/{sample}/{caller}/{caller}.{ext}"
    wildcard_constraints:
        ext = "(tsv|vcf)"
    threads: config['num_threads']
    conda: "envs/default.yml"
    shell:
        "mkdir -p \"{params.out_dir}\" && "
        "{input.caller_script} {input.bam} {input.peaks} "
        "{input.genome} {params.out_dir} {wildcards.sample} "
        "{threads} {input.shared} {params.caller_params}"


def sorted_cols(cols):
    """ sort and flatten the cols """
    return [
        c
        for col in ['other', 'info', 'format'] if col in cols
        for c in cols[col]
    ]


def caller_out(ext):
    # I'm embarassed by this code but it's the simplest way to achieve this
    if not ext:
        ext = 'vcf'
    return re.sub("{ext}$", ext, rules.run_caller.output.vcf)


rule filter_vcf:
    """ filter the vcf if needed """
    input:
        vcf = caller_out('vcf')
    output:
        vcf = pipe(caller_out('vcf')+".filter")
    conda: "envs/bcftools.yml"
    shell:
        "bcftools view {config[bcftools_params]} {input.vcf} > {output.vcf}"

rule normalize_vcf:
    """ normalize the vcf, split multiallelic sites, and remove duplicate sites if needed """
    input:
        vcf = rules.filter_vcf.output.vcf if 'bcftools_params' in config \
            and config['bcftools_params'] else caller_out('vcf'),
        ref = config['genome']
    output:
        vcf = pipe(caller_out('vcf')+".norm")
    conda: "envs/bcftools.yml"
    shell:
        "bcftools norm -m -any {input.vcf} | bcftools norm --check-ref xw -d "
        "all -f {input.ref} > {output.vcf}"


rule prepare_vcf:
    """ bgzip and index the vcf """
    input:
        vcf = rules.normalize_vcf.output.vcf if 'normalize' in config \
            and config['normalize'] else rules.filter_vcf.output.vcf \
            if 'bcftools_params' in config and config['bcftools_params'] else \
            caller_out('vcf')
    output:
        gzvcf = temp(caller_out('vcf')+".gz"),
        index = temp(caller_out('vcf')+".gz.tbi")
    conda: "envs/default.yml"
    shell:
        "bgzip <{input.vcf} >{output.gzvcf} && tabix -p vcf -f {output.gzvcf}"


def bcftools_query_str(wildcards):
    """ return the bcftools query string for this caller's columns """
    if wildcards.caller in config and config[wildcards.caller] and 'cols' in config[wildcards.caller]:
        cols = config[wildcards.caller]['cols'].copy()
    else:
        cols = {}
    if 'info' in cols:
        cols['info'] = ["INFO/"+c for c in cols['info']]
    if 'other' not in cols:
        cols['other'] = []
    cols['other'] = ['CHROM', 'POS', 'REF', 'ALT'] + cols['other']
    for col in cols:
        cols[col] = ["%"+c for c in cols[col]]
    if 'format' in cols:
        cols['format'] = ["["+c+"]" for c in cols['format']]
    return "\\t".join(sorted_cols(cols))+"\\n"


def cols_str(wildcards):
    col_str = "\\t".join(
        ['CHROM', 'POS', 'REF', 'ALT'] +
        sorted_cols(
            config[wildcards.caller]['cols']
            if wildcards.caller in config and config[wildcards.caller]
            and 'cols' in config[wildcards.caller]
            else {}
        )
    )
    if not hasattr(wildcards, 'hash'):
        return hash_str(col_str)
    return col_str


rule vcf2tsv:
    """Convert from vcf to tsv format, extracting relevant columns"""
    input:
        gzvcf = rules.prepare_vcf.output.gzvcf,
        index = rules.prepare_vcf.output.index
    params:
        cols = cols_str,
        qstr = bcftools_query_str
    output:
        tsv = config['output_dir'] + "/callers/{sample}/{caller}/{caller}.{hash}.tsv"
    conda: "envs/bcftools.yml"
    shell:
        "echo -e '{params.cols}' > {output.tsv} && "
        "bcftools query -f '{params.qstr}' {input.gzvcf} >> {output.tsv}"


def caller_tsv(wildcards):
    if wildcards.caller in config and 'ext' in config[wildcards.caller] and config[wildcards.caller]['ext'] == 'tsv':
        return expand(
            caller_out('tsv'),
            sample=wildcards.sample, caller=wildcards.caller
        )[0]
    return expand(
        rules.vcf2tsv.output.tsv, hash=cols_str(wildcards),
        sample=wildcards.sample, caller=wildcards.caller
    )[0]


rule prepare_merge:
    """
        1) add the caller as a prefix of every column name
        2) sort the file by CHROM and POS
        3) separate chrom and pos cols by comma instead of tab
        4) replace NA with .
        5) remove the header
        (not necessarily in that order)
    """
    input:
        tsv = caller_tsv
    output:
        pipe(config['output_dir'] + "/callers/{sample}/{caller}/prepared.tsv")
    conda: "envs/default.yml"
    shell:
        "tail -n+2 {input} | awk -F '\\t' -v 'OFS=\\t' '{{for (i=1; i<=NF; i++) if ($i==\"NA\") $i=\".\"}}1' | "
        "sed 's/\\t\+/,/' | LC_ALL=C sort -t $'\\t' -k1,1 > {output}"

rule get_all_sites:
    """retrieve all sites for output in the merged table"""
    input:
        rules.bed_peaks.output
    output:
        temp(config['output_dir'] + "/peaks/{sample}/all_sites.csv")
    shell:
        "awk '{{printf(\"%s\\t%d\\t%d\\t%s\\n\",$1,int($2)+1,int($3),$4);}}' {input} | "
        "awk -F '\\t' -v 'OFS=\\t' '{{for (i=$2;i<$3;i++) print $1\",\"i,substr($4,i-$2+1,1)}}' | "
        "LC_ALL=C sort -t $'\\t' -k1,1 > {output}"

rule join_all_sites:
    """
        1) add all sites to the prepared caller output using an outer join
        2) rename the column headers so we know which caller they came from
        3) get rid of the CHROM and POS cols so we can merge later
        4) rename the ref and alt columns as REF and ALT (to standardize them)
        (not necessarily in that order)
    """
    input:
        sites = rules.get_all_sites.output,
        tsv = caller_tsv,
        prepared_tsv = rules.prepare_merge.output
    output:
        pipe(config['output_dir'] + "/merged_{type}/{sample}.{caller}.tsv")
    conda: "envs/default.yml"
    shell:
        "LC_ALL=C join -t $'\\t' -e. -a1 -j1 -o auto --nocheck-order "
        "<(cut -f 1 {input.sites}) {input.prepared_tsv} | cut -f 2- | cat "
        "<(head -n1 {input.tsv} | cut -f 5- | tr '\\t' '\\n' | "
        "sed 's/^/{wildcards.caller}~/' | cat "
        "<(echo -e \'{wildcards.caller}~REF\\n{wildcards.caller}~ALT\') - | "
        "paste -s) - > {output}"

rule merge_callers:
    """merge the columns of each snp caller into a single file"""
    input:
        all_sites = rules.get_all_sites.output,
        caller_output = lambda wildcards: expand(
            rules.join_all_sites.output,
            caller=config[wildcards.type+'_callers'],
            sample=wildcards.sample,
            type=wildcards.type
        )
    output:
        config['output_dir'] + "/merged_{type}/{sample}.tsv.gz"
    conda: "envs/default.yml"
    shell:
        "paste <(echo -e 'CHROM\\tPOS\\tREF'; sed 's/,/\\t/' "
        "{input.all_sites}) {input.caller_output} | "
        "(read -r head; echo \"$head\"; sort -t $'\\t' -k1,1V -k2,2n) | gzip > {output}"
