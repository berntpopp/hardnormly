import os


# ---------------------------------------------------------------------------
# Config shortcuts
# ---------------------------------------------------------------------------

REF = config["ref"]["fasta"]
GENOME_FILE = config["ref"].get("genome_file", "")

INCLUDE_BEDS = config.get("regions", {}).get("include_beds", [])
EXCLUDE_BEDS = config.get("regions", {}).get("exclude_beds", [])
SLOP = config.get("regions", {}).get("slop", 100)

CALLER = config.get("filtering", {}).get("caller", "")
FILTERS_FILE = config.get("filtering", {}).get("filters_file", "")
STRIP_ANNOTATIONS = config.get("filtering", {}).get("strip_annotations", "")
ONLY_PASS = config.get("filtering", {}).get("only_pass", False)

OUTPUT_DIR = config["paths"]["output_folder"]
LOG_DIR = os.path.join(OUTPUT_DIR, config["paths"].get("log_subdir", "logs"))

GENERATE_STATS = config.get("processing", {}).get("generate_stats", True)
AUTO_INDEX = config.get("processing", {}).get("auto_index", True)
PLOT_STATS = config.get("processing", {}).get("plot_stats", False)
PLOT_OUTPUT_DIR = config.get("processing", {}).get(
    "plot_output_dir", os.path.join(OUTPUT_DIR, "plots")
)


# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------

def get_vcf_jobs():
    """Read VCF list file and return dict of {basename: vcf_path}."""
    vcf_list_file = config["paths"]["vcf_list"]
    jobs = {}
    with open(vcf_list_file) as fh:
        for line in fh:
            vcf = line.strip()
            if not vcf or vcf.startswith("#"):
                continue
            vcf_basename = os.path.basename(vcf).replace(".vcf.gz", "")
            jobs[vcf_basename] = vcf
    return jobs


def get_hardnormly_script():
    """Resolve path to hardnormly.sh relative to the repository root."""
    return os.path.join(workflow.basedir, "..", "hardnormly.sh")
