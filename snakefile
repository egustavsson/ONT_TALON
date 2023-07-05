import os
import Bio
import shutil
from os import path
from re import search
from pathlib import Path

from snakemake.utils import validate
from snakemake.utils import min_version
min_version("7.18")

# ----------------------------------------------------------------

configfile: "config.yml"
# validate(config=config, schema="schema/config_schema.yaml")
workdir: config["workdir"]

WORKDIR = config["workdir"]
SNAKEDIR = path.dirname(workflow.snakefile)

shutil.copy2(SNAKEDIR + "/config.yml", WORKDIR)

sample = config["sample_name"]


# ----------------------------------------------------------------

target_list = [
    "Nanostat/stat_out.txt",
    path.join("Talon", f"{sample}_talon_abundance_filtered.tsv"),
    path.join("Talon", f"{sample}_talon.gtf")
]

if config.get("reads_fastq") == "":
    target_list.remove("Nanostat/stat_out.txt")

rule all:
    input:
        target_list

# ----------------------------------------------------------------

rule nanostat:
    input:
        fq=config["reads_fastq"]

    output:
        ns="Nanostat/stat_out.txt"

    threads: config["threads"]

    shell:
        """
        NanoStat -n {output.ns} -t {threads} --tsv --fastq {input.fq}
        """


# ----------------------------------------------------------------

rule pychopper:
    input:
        fq=config["reads_fastq"]

    output:
        pyfq=path.join("Pychopper", f"{sample}_full_length_reads.fastq")

    params:
        outpath="Pychopper",
        prefix=f"{sample}_full_length_reads.fastq",
        pc="True" if config["run_pychopper"] else "False",
        pc_opts=config["pychopper_opts"]

    log: path.join(WORKDIR, "Pychopper", f"{sample}_pychopped.log")

    threads: config["threads"]

    run:
        if params.pc == "True":
            shell(
                """
                cd {params.outpath};
                pychopper {params.pc_opts} -t {threads} -r report.pdf -u unclassified.fq -w rescued.fq {input.fq} {params.prefix} 2> {log}
                """
            )
        else:
            shell(
                """
                ln -s `realpath {input.fq}` {output.pyfq}
                """
            )

# ----------------------------------------------------------------

rule minimap_mapping:
    input:
        genome=config["genome"],
        fq=rules.pychopper.output.pyfq

    output:
        sam=path.join("Mapping", f"{sample}_minimap.sam")

    params:
        opts=config["minimap2_opts"]

    log: path.join("Mapping", f"{sample}_minimap.log")

    threads: config["threads"]

    shell:
        """
        minimap2 -ax splice {params.opts} --secondary=no -t {threads} {input.genome} {input.fq} \
        | samtools sort -@ {threads} -o {output.sam}
        """
# ----------------------------------------------------------------

rule aln_stats:
    input:
        sam=rules.minimap_mapping.output.sam

    output:
        tsv=path.join("Mapping", f"{sample}_alignment_stats.tsv")

    threads: config["threads"]

    shell:
        """
        samtools stats -@ {threads} {input.sam} > {output.tsv}
        """
# ----------------------------------------------------------------

rule talon_initialize_database:
    input:
        gtf=config["gtf"]

    output:
        db=path.join("Talon", f"talon_{Path(config['gtf']).stem}.db")

    params:
        gtf_name=Path(config["gtf"]).stem,
        gbuild="hg38",
        prefix=path.join("Talon", f"talon_{Path(config['gtf']).stem}")

    log: path.join("Talon", f"{sample}_talon_db.log")

    shell:
        """
        (talon_initialize_database --f {input.gtf} --a {params.gtf_name} --g {params.gbuild} --o {params.prefix}) &> {log}
        """

# ----------------------------------------------------------------

rule talon_label_reads:
    input:
        sam=rules.minimap_mapping.output.sam,
        fa=config["genome"]

    output:
        sam=path.join("Talon", f"{sample}_labeled.sam")

    params:
        prefix=path.join("Talon", sample)

    log: path.join("Talon", f"{sample}_talon_labelReads.log")

    threads: config["threads"]

    shell:
        """
        (talon_label_reads --f {input.sam} --g {input.fa} --t {threads} --deleteTmp --o {params.prefix}) &> {log}
        """
# ----------------------------------------------------------------

rule talon_annotate:
    input:
        sam=rules.talon_label_reads.output.sam,
        talon_db=rules.talon_initialize_database.output.db

    output:
        qc=path.join("Talon", f"{sample}_QC.log"),
        tsv=path.join("Talon", f"{sample}_talon_read_annot.tsv")

    params:
        prefix1=sample,
        prefix2=path.join("Talon", sample),
        description=f"{sample}_sam",
        gbuild="hg38"

    log: path.join("Talon", f"{sample}_talon_annotate.log")

    threads: config["threads"]

    shell:
        """
        echo '{params.prefix1},{params.description},ONT,{input.sam}' > Talon/config_file.txt &&
        (talon --f Talon/config_file.txt --db {input.talon_db} --build {params.gbuild} --threads {threads} --o {params.prefix2}) &> {log}
        """
# ----------------------------------------------------------------

rule talon_filter_transcripts:
    input:
        talon_db=rules.talon_initialize_database.output.db,
        job_hold=rules.talon_annotate.output.tsv

    output:
        csv=path.join("Talon", f"{sample}_whitelist.csv")

    params:
        gtf_name=Path(config["gtf"]).stem,
        min_count=config["MIN_COUNT"]

    log: path.join("Talon", f"{sample}_talon_filter.log")

    shell:
        """
        (talon_filter_transcripts --db {input.talon_db} -a {params.gtf_name} --maxFracA 0.5 --minCount {params.min_count} --o {output.csv}) &> {log}
        """
# ----------------------------------------------------------------

rule talon_abundance:
    input:
        talon_db=rules.talon_initialize_database.output.db,
        whitelist=rules.talon_filter_transcripts.output.csv

    output:
        tsv=path.join("Talon", f"{sample}_talon_abundance_filtered.tsv")

    params:
        gtf_name=Path(config["gtf"]).stem,
        gbuild="hg38",
        prefix=path.join("Talon", sample)

    log: path.join("Talon", f"{sample}_talon_abundance.log")

    shell:
        """
        (talon_abundance --db {input.talon_db} -a {params.gtf_name} --build {params.gbuild} --whitelist {input.whitelist} --o {params.prefix}) &> {log}
        """
# ----------------------------------------------------------------

rule talon_create_GTF:
    input:
        talon_db=rules.talon_initialize_database.output.db,
        whitelist=rules.talon_filter_transcripts.output.csv

    output:
        talon_gtf=path.join("Talon", f"{sample}_talon.gtf")

    params:
        gtf_name=Path(config["gtf"]).stem,
        gbuild="hg38",
        prefix=path.join("Talon", sample)

    log: path.join("Talon", f"{sample}_talon_gtf.log")

    shell:
        """
        (talon_create_GTF --db {input.talon_db} -a {params.gtf_name} --build {params.gbuild} --whitelist {input.whitelist} --o {params.prefix}) &> {log}
        """