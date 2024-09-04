import os
import yaml  # or json depending on config format

# ----------------------------------------------------------------------------------- #
# Load VCF file paths from input file (vcfs.txt)
vcf_files = [line.strip() for line in open('vcfs.txt')]

# Load config for other files like reference, filters, bed files, etc.
config = yaml.safe_load(open('config.yaml'))  # Alternatively, use json.load() for JSON

# Create a dictionary to store the jobs
jobs = {}
for vcf in vcf_files:
    vcf_basename = os.path.basename(vcf).replace(".vcf.gz", "")
    output_vcf = f"results/{vcf_basename}.vcf.gz"
    stats_dir = f"stats/{vcf_basename}/"
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
        expand("results/{vcf_basename}.vcf.gz", vcf_basename=[jobs[key]['vcf_basename'] for key in jobs.keys()])

rule run_hardnormly_pipeline:
    input:
        vcf=lambda wildcards: jobs[wildcards.vcf_basename]['vcf'],
    output:
        vcf="results/{vcf_basename}.vcf.gz",
    params:
        fasta=config['reference_fasta'],
        include_beds=" ".join(f"-b {bed}" for bed in config['include_beds']),
        exclude_beds=" ".join(f"-e {bed}" for bed in config['exclude_beds']),
        filters_file=config['filters_file'],
        genome_file=config['genome_file'],
        stats_dir=lambda wildcards: jobs[wildcards.vcf_basename]['stats_dir']
    threads: 8
    resources:
        mem_mb=16000,
        time="72:00:00"
    shell:
        """
        ./hardnormly.sh -v {input.vcf} -f {params.fasta} {params.include_beds} \
        {params.exclude_beds} --filters-file {params.filters_file} \
        --generate-stats -g {params.genome_file} -o {output.vcf} \
        --only-pass --plot-stats --plot-output-dir {params.stats_dir}
        """
