import os
import yaml

# ----------------------------------------------------------------------------------- #
# Load VCF file paths from input file (vcfs.txt)
vcf_files = [line.strip() for line in open('vcfs.txt')]

# Load config for other files like reference, filters, bed files, etc.
config = yaml.safe_load(open('config.yaml'))

# Create a dictionary to store the jobs
jobs = {}
for vcf in vcf_files:
    vcf_basename = os.path.basename(vcf).replace(".vcf.gz", "")
    output_vcf = os.path.join(config['output_dir'], f"{vcf_basename}.hardnormly.vcf.gz")
    stats_dir = os.path.join(config['output_dir'], vcf_basename, "stats")
    jobs[vcf_basename] = {
        "vcf": vcf,
        "vcf_basename": vcf_basename,
        "output_vcf": output_vcf,
        "stats_dir": stats_dir
    }

# ----------------------------------------------------------------------------------- #
# Define the rules
rule all:
    input:
        expand("logs/{vcf_basename}.log", vcf_basename=[jobs[key]['vcf_basename'] for key in jobs.keys()])

rule run_hardnormly_pipeline:
    input:
        vcf=lambda wildcards: jobs[wildcards.vcf_basename]['vcf'],
    output:
        log="logs/{vcf_basename}.log",
    params:
        hardnormly_script=config['hardnormly_script'],  # Path to the hardnormly script
        fasta=config['reference_fasta'],
        include_beds=" ".join(f"-b {bed}" for bed in config['include_beds']),
        exclude_beds=" ".join(f"-e {bed}" for bed in config['exclude_beds']),
        filters_file=config['filters_file'],
        genome_file=config['genome_file'],
        stats_dir=lambda wildcards: jobs[wildcards.vcf_basename]['stats_dir'],
        output_vcf=lambda wildcards: jobs[wildcards.vcf_basename]['output_vcf']
    threads: 8
    resources:
        mem_mb=8000,
        time="72:00:00"
    conda:
        "hardnormly"
    shell:
        """
        {params.hardnormly_script} -v {input.vcf} -f {params.fasta} {params.include_beds} \
        {params.exclude_beds} --filters-file {params.filters_file} \
        --generate-stats -g {params.genome_file} -o {params.output_vcf} &> {output.log}
        """
