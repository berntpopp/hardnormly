import os
import yaml

# ----------------------------------------------------------------------------------- #
# Load config for other files like reference, filters, bed files, etc.
config = yaml.safe_load(open('config.yaml'))

# ----------------------------------------------------------------------------------- #
# Ensure the main output directory exists
os.makedirs(config['output_dir'], exist_ok=True)

# Create a logs subfolder inside the output directory
log_dir = os.path.join(config['output_dir'], "logs")
os.makedirs(log_dir, exist_ok=True)

# ----------------------------------------------------------------------------------- #
# Load VCF file paths from the file specified in config['vcf_files']
vcf_files = [line.strip() for line in open(config['vcf_files'])]

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
        expand(f"{log_dir}/{{vcf_basename}}.log", vcf_basename=[jobs[key]['vcf_basename'] for key in jobs.keys()])

rule run_hardnormly_pipeline:
    input:
        vcf=lambda wildcards: jobs[wildcards.vcf_basename]['vcf'],
    output:
        log = f"{log_dir}/{{vcf_basename}}.log",
    params:
        hardnormly_script=config['hardnormly_script'],  # Path to the hardnormly script
        fasta=config['reference_fasta'],
        include_beds=" ".join(f"-b {bed}" for bed in config['include_beds']),
        exclude_beds=" ".join(f"-e {bed}" for bed in config['exclude_beds']),
        filters_file=config['filters_file'],
        genome_file=config['genome_file'],
        slop=config['slop'],
        stats_dir=lambda wildcards: jobs[wildcards.vcf_basename]['stats_dir'],
        output_vcf=lambda wildcards: jobs[wildcards.vcf_basename]['output_vcf']
    threads: 2
    resources:
        mem_mb=8000,
        time="8:00:00"
    conda:
        "hardnormly"
    shell:
        """
        {params.hardnormly_script} -v {input.vcf} -f {params.fasta} {params.include_beds} \
        {params.exclude_beds} --filters-file {params.filters_file} \
        --generate-stats -g {params.genome_file} --slop {params.slop} \
        -o {params.output_vcf} &> {output.log}
        """
