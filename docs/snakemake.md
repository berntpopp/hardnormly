# Snakemake Workflow

For batch processing multiple VCF files in parallel, hardnormly includes a Snakemake 8+ workflow with SLURM support.

## Overview

The workflow reads a list of VCF file paths and runs `hardnormly.sh` on each one independently. Configuration (reference paths, BED files, filters) is shared across all samples via `config.yaml`.

## Setup

1. **Edit configuration**: Adjust `config/config.yaml` with your reference paths, BED files, and filter settings.

2. **Create VCF list**: Create a text file with one VCF path per line. Set `paths.vcf_list` in the config to point to this file.

3. **Choose profile**: Use `profiles/default/` for local execution or create a cluster profile (see `profiles/charite/` for a SLURM example).

## Running

```bash
# Dry run — preview what will be executed
snakemake --snakefile workflow/Snakefile --configfile config/config.yaml \
  --workflow-profile profiles/default -n

# Local execution
snakemake --snakefile workflow/Snakefile --configfile config/config.yaml \
  --workflow-profile profiles/default

# SLURM submission (auto-detects cluster environment)
sbatch scripts/run_snakemake.sh
```

## Workflow Structure

```
workflow/
├── Snakefile                # Entry point (config validation, rule imports)
├── rules/
│   ├── common.smk           # Config shortcuts and helper functions
│   └── hardnormly.smk       # Pipeline rule (calls hardnormly.sh per sample)
├── envs/
│   └── hardnormly.yaml      # Lightweight conda env for the pipeline rule
└── schemas/
    └── config.schema.yaml   # JSON Schema for config validation

config/
└── config.yaml              # Workflow configuration

profiles/
├── default/
│   └── config.yaml          # Default resources and execution settings
└── charite/
    └── config.yaml          # Charite SLURM cluster settings
```

## Configuration

Key sections in `config/config.yaml`:

| Section | What it controls |
|---------|-----------------|
| `ref` | Reference genome FASTA and genome build |
| `paths` | Input VCF list, output directory, log directory |
| `regions` | Include/exclude BED files, slop value |
| `filtering` | Filter file path, PASS-only flag |
| `processing` | Stats generation, auto-indexing, plot generation |

Config is validated against `workflow/schemas/config.schema.yaml` at workflow start. Invalid configuration fails fast with a clear error message.

## Profiles

Profiles control resource allocation (threads, memory, runtime) without modifying rules:

- **`profiles/default/`** — Sensible defaults for local execution
- **`profiles/charite/`** — SLURM cluster settings (partition, account, walltime)

Create your own profile by copying `profiles/default/` and adjusting values for your cluster.

## Environment

The workflow uses Snakemake's `--use-conda` feature. Only a `snakemake` (8+) base environment is needed — per-rule conda environments are created automatically from `workflow/envs/hardnormly.yaml`.
