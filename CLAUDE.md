# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

**hardnormly** is a bioinformatics toolkit (v0.6.0) for VCF normalization and hard filtering. It uses `bcftools` for variant filtering/normalization and `bedtools` for region-based operations on genomic intervals. Designed for whole-exome sequencing (WES) variant processing pipelines.

## Architecture

The project has two execution modes:

1. **Standalone script** (`hardnormly.sh`) — the core bash tool that performs the full pipeline: BED normalization → genome file creation → region annotation → VCF normalization → hard filtering → optional stats/plots.

2. **Snakemake workflow** (`workflow/Snakefile`) — batch-processes multiple VCF files in parallel via SLURM, calling `hardnormly.sh` for each sample. Uses Snakemake 8+ conventions with config validation, profile-based resources, and per-rule conda environments.

### Pipeline Steps (hardnormly.sh)

1. Create genome file (from UCSC MySQL if not provided)
2. Normalize and merge include/exclude BED files (bedtools intersect, slop, multiinter)
3. Compress and index merged BED files (bgzip + tabix)
4. Annotate VCF with INCLUDE_REGION/EXCLUDE_REGION INFO fields (bcftools annotate)
5. Normalize VCF (bcftools norm with multiallelic splitting)
6. Apply filters: region-based → inline/file-based → optional PASS-only (bcftools filter pipeline built via eval)
7. Optional: generate stats (bcftools stats) and plots (plot-vcfstats)

### Filter File Format

Filter files (e.g., `defaults/gatk_filters.txt`) use three space-separated columns per line:
```
<filter_name> <action: e=exclude/i=include> <bcftools_expression>
```
Example: `DPu10het e FORMAT/DP<10 && GT!="hom"`

Two default filter sets are provided: `defaults/gatk_filters.txt` (GATK HaplotypeCaller) and `defaults/freebayes_filters.txt`.

## Directory Structure

```
hardnormly/
├── hardnormly.sh                    # Main pipeline script (standalone)
├── defaults/                        # Default filter files and hg19 genome file
├── ref/                             # Reference FASTA, target/exclusion BED files
├── conda/                           # Full conda environment (standalone usage)
│   └── hardnormly_environment.yml
├── workflow/                        # Snakemake workflow (standard layout)
│   ├── Snakefile                    # Entry point (min_version 8.0, config validation)
│   ├── rules/
│   │   ├── common.smk              # Config shortcuts and helper functions
│   │   └── hardnormly.smk          # Pipeline rule (calls hardnormly.sh)
│   ├── envs/
│   │   └── hardnormly.yaml         # Lightweight conda env (bcftools, bedtools, htslib)
│   └── schemas/
│       └── config.schema.yaml      # JSON Schema for config validation
├── config/
│   └── config.yaml                 # Workflow configuration (hierarchical)
├── profiles/
│   ├── default/
│   │   └── config.yaml             # Default resources and execution settings
│   └── charite/
│       └── config.yaml             # Charite SLURM cluster settings
└── scripts/
    └── run_snakemake.sh            # Launcher with cluster auto-detection
```

## Key Files

- `hardnormly.sh` — main script, all pipeline logic
- `workflow/Snakefile` — Snakemake workflow entry point
- `workflow/rules/common.smk` — config shortcuts, VCF list loading helpers
- `workflow/rules/hardnormly.smk` — pipeline rule definition
- `workflow/schemas/config.schema.yaml` — config validation schema
- `config/config.yaml` — workflow configuration (reference paths, BED files, filters, output)
- `profiles/default/config.yaml` — resource allocations and execution settings
- `conda/hardnormly_environment.yml` — full conda environment specification (standalone)
- `workflow/envs/hardnormly.yaml` — lightweight conda env (Snakemake per-rule)

## Environment Setup

```bash
# For standalone hardnormly.sh usage (full environment with plotting tools)
conda env create -f conda/hardnormly_environment.yml
conda activate hardnormly

# For Snakemake workflow (manages per-rule envs automatically)
conda activate snakemake  # needs snakemake 8+
```

Key dependencies: bcftools >=1.21 (1.21 required for `--write-index=tbi`), bedtools 2.31, htslib (bgzip/tabix), mysql client, matplotlib, tectonic.

## Running

### Standalone
```bash
./hardnormly.sh -v input.vcf.gz -f ref.fasta -b include.bed -e exclude.bed \
  --filters-file defaults/gatk_filters.txt -g defaults/hg19.genome -o output.vcf.gz
```

### Snakemake (local)
```bash
# From repository root
snakemake --snakefile workflow/Snakefile --configfile config/config.yaml \
  --workflow-profile profiles/default -n  # dry run
```

### Snakemake (SLURM)
```bash
# Submit via launcher script (auto-detects cluster)
sbatch scripts/run_snakemake.sh

# Or run directly with cluster profile
snakemake --snakefile workflow/Snakefile --configfile config/config.yaml \
  --workflow-profile profiles/default --profile profiles/charite
```

The workflow reads VCF paths from the file specified in `config/config.yaml` → `paths.vcf_list` (one path per line).

### Config Structure

The config (`config/config.yaml`) uses hierarchical sections:
- `ref` — reference genome paths and build
- `paths` — input VCF list, output directory, log subdirectory
- `regions` — include/exclude BED files, slop value
- `filtering` — filter definitions file, PASS-only flag
- `processing` — stats generation, auto-indexing, plot generation

Config is validated against `workflow/schemas/config.schema.yaml` at workflow start.

## Development Notes

- The default genome build is hg19 (GRCh37); BED files under `ref/target_files/` include both hg19 and hg38 versions, with `_nochr` variants for chromosome naming without "chr" prefix.
- The script builds a filter pipeline as a string and executes it via `eval` — be careful with quoting in filter expressions.
- `--debug` enables `set -x` tracing and verbose logging throughout.
- `--no-cleanup` preserves the temp directory for inspecting intermediate files.
- The `stats/plot.py` file is auto-generated by `plot-vcfstats`, not hand-written code.
- Resources (threads, memory, runtime) are managed via profiles, not hardcoded in rules.
- The Snakemake workflow tracks actual output VCFs (not just logs) for proper dependency tracking.
