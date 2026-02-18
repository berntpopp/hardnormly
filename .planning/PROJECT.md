# hardnormly

## What This Is

A bioinformatics toolkit for VCF normalization and hard filtering in whole-exome sequencing (WES) variant processing pipelines. Uses bcftools for variant filtering/normalization and bedtools for region-based operations. Runs as a standalone bash script or via Snakemake workflow for batch processing. The script is modularized into 8 lib/ function libraries with a comprehensive BATS test suite (130 tests) and full CI/linting infrastructure.

## Core Value

Reliably normalize and filter VCF files so downstream clinical/research analysis starts from clean, consistent variant calls.

## Requirements

### Validated

- ✓ VCF normalization via bcftools norm (multiallelic splitting, left-alignment) — v0.6.0
- ✓ Region-based annotation (INCLUDE_REGION/EXCLUDE_REGION INFO fields via bcftools annotate) — v0.6.0
- ✓ Hard filtering with caller-specific filter files (GATK HaplotypeCaller, Freebayes) — v0.6.0
- ✓ BED file normalization and merging (bedtools intersect, slop, multiinter) — v0.6.0
- ✓ Genome file creation from UCSC MySQL — v0.6.0
- ✓ Optional stats generation (bcftools stats) and plotting (plot-vcfstats) — v0.6.0
- ✓ Compressed/indexed output support (bgzip + tabix) — v0.6.0
- ✓ Snakemake workflow for batch processing via SLURM — v0.6.0
- ✓ Config validation via JSON Schema — v0.6.0
- ✓ Profile-based resource management (default, charite) — v0.6.0
- ✓ ShellCheck + shfmt linting with CI enforcement — v0.7.0
- ✓ Bash strict mode (`set -euo pipefail`) — v0.7.0
- ✓ GitHub Actions CI (lint + test workflows) — v0.7.0
- ✓ Pre-commit hooks (shellcheck, shfmt, trailing whitespace) — v0.7.0
- ✓ Synthetic test VCFs (GATK + Freebayes styles, edge cases) — v0.7.0
- ✓ Real data test subsets (GIAB, 1000 Genomes, Freebayes tiny) — v0.7.0
- ✓ BATS test framework (smoke, unit, integration, regression) — v0.7.0
- ✓ Modularization into lib/ (logging, cli, genome, bed, annotate, normalize, filter, stats) — v0.7.0
- ✓ Strict error handling throughout — v0.7.0
- ✓ Replace eval pipeline construction — v0.7.0
- ✓ Plot-vcfstats non-fatal error handling — v0.7.0
- ✓ Subcommands (generate-inclusion-bed, generate-exclusion-bed, run-pipeline) — v0.7.0
- ✓ Combined exclusion BED generation script — v0.7.0
- ✓ `--caller` flag for auto-selecting filter file — v0.7.0
- ✓ Enhanced show_help() with examples and filter format — v0.7.0
- ✓ Exclusion BED source documentation — v0.7.0

### Active

(None — define with `/gsd:new-milestone`)

### Out of Scope

- Snakemake workflow refactoring — works but needs cluster testing on nserver separately
- hg38 as default genome build — hg19/GRCh37 remains default; hg38 support exists for BED files
- GUI or web interface — CLI tool only
- VCF annotation beyond regions — not a variant annotator (use VEP/SnpEff for that)
- Real clinical test data in repo — replaced with synthetic + public data for reproducibility
- Python replacement — staying in bash; the ecosystem (bcftools, bedtools) is shell-native

## Context

- Shipped v0.7.0 with ~14,043 LOC (bash + bats), 169 files
- Tech stack: Bash, bcftools, bedtools, htslib, BATS, GitHub Actions, GNU Make
- 130 tests pass (8 BATS files): smoke, GATK filter units, Freebayes filter units, integration, regression, lib/ unit tests
- Main script: hardnormly.sh (~219 lines orchestrator) + 8 lib/ modules
- Issue #13 (genome flag) verified closed by TFWK-06 tests
- Two default filter sets: defaults/gatk_filters.txt and defaults/freebayes_filters.txt
- One deferred tech debt item: hardcoded PATH in tests/data/generate_expected.sh (dev-only, non-blocking)

## Constraints

- **Tech stack**: Bash (bcftools, bedtools, htslib) — no language change
- **Testing**: BATS framework (bash-native testing)
- **Linting**: ShellCheck + shfmt (industry standard for bash)
- **CI**: GitHub Actions
- **Test data size**: <5 MB total, no Git LFS
- **Backward compatibility**: Existing CLI interface keeps working (subcommands are additive)
- **Focus**: Main script only — Snakemake workflow changes out of scope

## Key Decisions

| Decision | Rationale | Outcome |
|----------|-----------|---------|
| Focus on main script, not Snakemake | Snakemake needs cluster testing; script is the core | ✓ Good — delivered full quality + features milestone without cluster dependency |
| Replace real test data with synthetic + public | No clinical data in repo; reproducible tests | ✓ Good — 1000G + GIAB + synthetic covers all filter triggers |
| BATS for testing framework | Bash-native, well-supported, standard for bash projects | ✓ Good — 130 tests, clear filter-level assertions |
| lib/ modularization pattern | Source-able functions, each module independently testable | ✓ Good — 8 modules, all with unit tests, main script ~219 lines |
| Array-based filter pipeline (no eval) | Safer, shellcheck-clean, eliminates string injection risk | ✓ Good — `|` delimiter works; filter_stages array pattern established |
| Subcommands as additive (backward compat) | Existing users unaffected; new capabilities discoverable | ✓ Good — legacy invocation still works; `-*` falls through to parse_args |
| `.githooks/` over `.git/hooks/` | Tracked in git, shareable; `setup-hooks.sh` for activation | ✓ Good — hooks shareable across developers |
| ERR trap removed | cleanup_handler trap is sufficient; ERR caused false-positive exits | ✓ Good — EXIT trap-only pattern cleaner |
| shfmt flags `-i 0 -bn -ci` | Tabs, binary ops at line start, case indent — project standard | ✓ Good — consistent throughout codebase |

---
*Last updated: 2026-02-18 after v0.7.0 milestone*
