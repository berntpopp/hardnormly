# hardnormly

## What This Is

A bioinformatics toolkit for VCF normalization and hard filtering in whole-exome sequencing (WES) variant processing pipelines. Uses bcftools for variant filtering/normalization and bedtools for region-based operations. Runs as a standalone bash script or via Snakemake workflow for batch processing.

## Core Value

Reliably normalize and filter VCF files so downstream clinical/research analysis starts from clean, consistent variant calls.

## Current Milestone: v0.7.0 — Code Quality, Testing & Features

**Goal:** Transform hardnormly from a monolithic script with no tests or CI into a well-structured, tested, linted, and documented tool with new capabilities.

**Target features:**
- Linting infrastructure (ShellCheck, shfmt, CI, pre-commit)
- Bash strict mode (`set -euo pipefail`)
- Synthetic + public test data suite
- BATS automated test framework
- Modularization into `lib/` function libraries
- Safer pipeline construction (replace `eval`)
- Plot-vcfstats error handling (non-fatal)
- Subcommands (generate-inclusion-bed, generate-exclusion-bed, run-pipeline)
- Combined exclusion BED generation
- `--caller` flag for auto-selecting filter files
- Enhanced help text and documentation

## Requirements

### Validated

- VCF normalization via bcftools norm (multiallelic splitting, left-alignment)
- Region-based annotation (INCLUDE_REGION/EXCLUDE_REGION INFO fields via bcftools annotate)
- Hard filtering with caller-specific filter files (GATK HaplotypeCaller, Freebayes)
- BED file normalization and merging (bedtools intersect, slop, multiinter)
- Genome file creation from UCSC MySQL
- Optional stats generation (bcftools stats) and plotting (plot-vcfstats)
- Compressed/indexed output support (bgzip + tabix)
- Snakemake workflow for batch processing via SLURM
- Config validation via JSON Schema
- Profile-based resource management (default, charite)

### Active

- [ ] ShellCheck + shfmt linting with CI enforcement
- [ ] Bash strict mode (`set -euo pipefail`)
- [ ] GitHub Actions CI (lint + test workflows)
- [ ] Pre-commit hooks (shellcheck, shfmt, trailing whitespace)
- [ ] Synthetic test VCFs (GATK + Freebayes styles, edge cases)
- [ ] Real data test subsets (GIAB, 1000 Genomes, Freebayes tiny)
- [ ] BATS test framework (smoke, unit, integration, regression)
- [ ] Modularization into lib/ (logging, cli, genome, bed, annotate, normalize, filter, stats)
- [ ] Strict error handling throughout
- [ ] Replace eval pipeline construction
- [ ] Plot-vcfstats non-fatal error handling
- [ ] Subcommands (generate-inclusion-bed, generate-exclusion-bed, run-pipeline)
- [ ] Combined exclusion BED generation script
- [ ] `--caller` flag for auto-selecting filter file
- [ ] Enhanced show_help() with examples
- [ ] Exclusion BED source documentation

### Out of Scope

- Snakemake workflow refactoring — works but needs cluster testing separately
- hg38 default support — hg19/GRCh37 remains default
- GUI or web interface — CLI tool only
- VCF annotation beyond region-based — not a variant annotator
- Real clinical test data in repo — replaced with synthetic + public

## Context

- Current script is v0.6.0, 463 lines, remarkably clean for ShellCheck (5 findings only)
- Three existing functions: `cleanup_tmp_dir`, `show_help`, `normalize_bed`; rest is procedural
- Filter pipeline built as string and executed via `eval` (L401) — works but fragile
- Snakemake workflow recently restructured (new config/, profiles/, workflow/ layout) — not in scope
- Existing `test/` directory has real clinical VCFs — will be replaced with synthetic + public data
- Two default filter sets: `defaults/gatk_filters.txt` and `defaults/freebayes_filters.txt`
- Issue #13 (genome file option) already implemented — verify with test and close
- Issues #8/#18, #5/#14, #15/#6 are duplicate pairs

## Constraints

- **Tech stack**: Bash (bcftools, bedtools, htslib) — no language change
- **Testing**: BATS framework (bash-native testing)
- **Linting**: ShellCheck + shfmt (industry standard for bash)
- **CI**: GitHub Actions
- **Test data size**: <5 MB total, no Git LFS
- **Backward compatibility**: Existing CLI interface must keep working (subcommands are additive)
- **Focus**: Main script only — Snakemake workflow changes out of scope

## Key Decisions

| Decision | Rationale | Outcome |
|----------|-----------|---------|
| Focus on main script, not Snakemake | Snakemake needs cluster testing; script is the core | -- Pending |
| Replace real test data with synthetic + public | No clinical data in repo; reproducible tests | -- Pending |
| BATS for testing framework | Bash-native, well-supported, standard for bash projects | -- Pending |
| lib/ modularization pattern | Source-able functions, each module independently testable | -- Pending |
| Subcommands: generate-inclusion-bed, generate-exclusion-bed, run-pipeline | Matches existing BED generation + pipeline steps | -- Pending |

---
*Last updated: 2026-02-18 after milestone v0.7.0 initialization*
