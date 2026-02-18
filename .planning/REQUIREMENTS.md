# Requirements: hardnormly

**Defined:** 2026-02-18
**Core Value:** Reliably normalize and filter VCF files for clean variant calls

## v1 Requirements

Requirements for milestone v0.7.0. Each maps to roadmap phases.

### Infrastructure

- [ ] **INFR-01**: Project has ShellCheck linting with `.shellcheckrc` config and all existing findings fixed
- [ ] **INFR-02**: Project has shfmt formatting applied with `.editorconfig` for consistent style
- [ ] **INFR-03**: Script runs with `set -euo pipefail` and targeted exception handling for expected failures
- [ ] **INFR-04**: GitHub Actions workflow runs ShellCheck + shfmt checks on every push/PR
- [ ] **INFR-05**: GitHub Actions workflow runs BATS tests on every push/PR
- [ ] **INFR-06**: Pre-commit hooks enforce ShellCheck, shfmt, and trailing-whitespace checks locally

### Testing — Data

- [ ] **TDAT-01**: Synthetic GATK-style VCFs exist (samples A and B) with variants covering all GATK filter triggers
- [ ] **TDAT-02**: Synthetic Freebayes-style VCFs exist (samples A and B) with variants covering all Freebayes filter triggers
- [ ] **TDAT-03**: Synthetic utility VCFs exist (multiallelic, minimal, empty_variants) for edge case testing
- [ ] **TDAT-04**: Synthetic reference FASTA and BED files exist for test region (chr22:16M-16.1M)
- [ ] **TDAT-05**: Real data subsets exist (GIAB NA12878, 1000 Genomes 3 samples, Freebayes tiny) with indexes
- [ ] **TDAT-06**: Test data generation script (`generate_test_data.sh`) can reproducibly create all real data subsets
- [ ] **TDAT-07**: Expected output files exist for regression testing (pre-validated pipeline output)
- [ ] **TDAT-08**: VCF list files and test configs exist for Snakemake batch testing

### Testing — Framework

- [ ] **TFWK-01**: BATS smoke tests verify --help exits 0, --version prints version, missing args exits non-zero
- [ ] **TFWK-02**: BATS filter tests verify each GATK filter correctly tags matching variants
- [ ] **TFWK-03**: BATS filter tests verify each Freebayes filter correctly tags matching variants
- [ ] **TFWK-04**: BATS integration test runs full pipeline on real data and produces valid output
- [ ] **TFWK-05**: BATS regression tests compare pipeline output against expected output files
- [ ] **TFWK-06**: BATS test verifies -g/--genome flag works (issue #13 — close after verification)
- [ ] **TFWK-07**: BATS edge case tests verify empty VCF, no filters, no BED files don't crash

### Refactoring

- [ ] **REFR-01**: Logging module (`lib/logging.sh`) provides log_msg, debug_msg, error handling functions
- [ ] **REFR-02**: CLI module (`lib/cli.sh`) handles argument parsing, validation, show_help
- [ ] **REFR-03**: Genome module (`lib/genome.sh`) handles genome file creation and validation
- [ ] **REFR-04**: BED module (`lib/bed.sh`) handles normalize_bed, merge, intersect, compress, index
- [ ] **REFR-05**: Annotate module (`lib/annotate.sh`) handles VCF annotation with BED regions
- [ ] **REFR-06**: Normalize module (`lib/normalize.sh`) wraps bcftools norm operations
- [ ] **REFR-07**: Filter module (`lib/filter.sh`) handles filter pipeline construction and execution
- [ ] **REFR-08**: Stats module (`lib/stats.sh`) handles stats generation and plotting
- [ ] **REFR-09**: Main script (`hardnormly.sh`) sources lib/ modules and orchestrates pipeline steps
- [ ] **REFR-10**: Filter pipeline uses safer construction (no `eval` on dynamically built strings)
- [ ] **REFR-11**: All external commands use consistent error handling pattern (`run_cmd` or equivalent)
- [ ] **REFR-12**: Duplicate cleanup logic eliminated (single trap-based pattern)
- [ ] **REFR-13**: Unit tests exist for each lib/ module's key functions

### Features

- [ ] **FEAT-01**: Plot-vcfstats errors are non-fatal (logged but don't exit the pipeline)
- [ ] **FEAT-02**: Subcommand dispatcher routes to generate-inclusion-bed, generate-exclusion-bed, or run-pipeline
- [ ] **FEAT-03**: `generate-inclusion-bed` subcommand produces a merged include BED from inputs
- [ ] **FEAT-04**: `generate-exclusion-bed` subcommand produces a merged exclusion BED from inputs
- [ ] **FEAT-05**: Helper script generates combined exclusion BED files for hg19 and hg38 from public sources
- [ ] **FEAT-06**: `--caller` flag auto-selects the matching filter file (gatk → gatk_filters.txt, freebayes → freebayes_filters.txt)
- [ ] **FEAT-07**: Running without a subcommand defaults to `run-pipeline` (backward compatible)

### Documentation

- [ ] **DOCS-01**: `--help` output includes usage examples section
- [ ] **DOCS-02**: `--help` output describes filter file format
- [ ] **DOCS-03**: Exclusion BED file sources documented (origin URL, genome build, regions covered)

## v2 Requirements

Deferred to future milestones. Tracked but not in current roadmap.

### Snakemake

- **SMKF-01**: Snakemake workflow uses context managers for file handles
- **SMKF-02**: Snakemake lint + snakefmt checks added to CI
- **SMKF-03**: Snakemake batch tests validate parallel processing without collisions

### Features

- **FEAT-08**: `--plot-backend` option to choose between tectonic and pdflatex
- **FEAT-09**: `--help-filters` subcommand prints filter file format in detail

## Out of Scope

| Feature | Reason |
|---------|--------|
| Snakemake workflow refactoring | Works but needs cluster testing on nserver separately |
| hg38 as default genome build | hg19/GRCh37 remains default; hg38 support exists for BED files |
| GUI or web interface | CLI tool only |
| VCF annotation beyond regions | Not a variant annotator (use VEP/SnpEff for that) |
| Real clinical test data in repo | Replaced with synthetic + public data for reproducibility |
| Python replacement | Staying in bash — the ecosystem (bcftools, bedtools) is shell-native |

## Traceability

| Requirement | Phase | Status |
|-------------|-------|--------|
| INFR-01 | TBD | Pending |
| INFR-02 | TBD | Pending |
| INFR-03 | TBD | Pending |
| INFR-04 | TBD | Pending |
| INFR-05 | TBD | Pending |
| INFR-06 | TBD | Pending |
| TDAT-01 | TBD | Pending |
| TDAT-02 | TBD | Pending |
| TDAT-03 | TBD | Pending |
| TDAT-04 | TBD | Pending |
| TDAT-05 | TBD | Pending |
| TDAT-06 | TBD | Pending |
| TDAT-07 | TBD | Pending |
| TDAT-08 | TBD | Pending |
| TFWK-01 | TBD | Pending |
| TFWK-02 | TBD | Pending |
| TFWK-03 | TBD | Pending |
| TFWK-04 | TBD | Pending |
| TFWK-05 | TBD | Pending |
| TFWK-06 | TBD | Pending |
| TFWK-07 | TBD | Pending |
| REFR-01 | TBD | Pending |
| REFR-02 | TBD | Pending |
| REFR-03 | TBD | Pending |
| REFR-04 | TBD | Pending |
| REFR-05 | TBD | Pending |
| REFR-06 | TBD | Pending |
| REFR-07 | TBD | Pending |
| REFR-08 | TBD | Pending |
| REFR-09 | TBD | Pending |
| REFR-10 | TBD | Pending |
| REFR-11 | TBD | Pending |
| REFR-12 | TBD | Pending |
| REFR-13 | TBD | Pending |
| FEAT-01 | TBD | Pending |
| FEAT-02 | TBD | Pending |
| FEAT-03 | TBD | Pending |
| FEAT-04 | TBD | Pending |
| FEAT-05 | TBD | Pending |
| FEAT-06 | TBD | Pending |
| FEAT-07 | TBD | Pending |
| DOCS-01 | TBD | Pending |
| DOCS-02 | TBD | Pending |
| DOCS-03 | TBD | Pending |

**Coverage:**
- v1 requirements: 43 total
- Mapped to phases: 0
- Unmapped: 43

---
*Requirements defined: 2026-02-18*
*Last updated: 2026-02-18 after initial definition*
