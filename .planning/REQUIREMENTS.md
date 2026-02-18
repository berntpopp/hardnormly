# Requirements: hardnormly

**Defined:** 2026-02-18
**Core Value:** Reliably normalize and filter VCF files for clean variant calls

## v1 Requirements

Requirements for milestone v0.7.0. Each maps to roadmap phases.

### Infrastructure

- [x] **INFR-01**: Project has ShellCheck linting with `.shellcheckrc` config and all existing findings fixed
- [x] **INFR-02**: Project has shfmt formatting applied with `.editorconfig` for consistent style
- [x] **INFR-03**: Script runs with `set -euo pipefail` and targeted exception handling for expected failures
- [x] **INFR-04**: GitHub Actions workflow runs ShellCheck + shfmt checks on every push/PR
- [x] **INFR-05**: GitHub Actions workflow runs BATS tests on every push/PR
- [x] **INFR-06**: Pre-commit hooks enforce ShellCheck, shfmt, and trailing-whitespace checks locally

### Testing — Data

- [x] **TDAT-01**: Synthetic GATK-style VCFs exist (samples A and B) with variants covering all GATK filter triggers
- [x] **TDAT-02**: Synthetic Freebayes-style VCFs exist (samples A and B) with variants covering all Freebayes filter triggers
- [x] **TDAT-03**: Synthetic utility VCFs exist (multiallelic, minimal, empty_variants) for edge case testing
- [x] **TDAT-04**: Synthetic reference FASTA and BED files exist for test region (chr22:16M-16.1M)
- [x] **TDAT-05**: Real data subsets exist (GIAB NA12878, 1000 Genomes 3 samples, Freebayes tiny) with indexes
- [x] **TDAT-06**: Test data generation script (`generate_test_data.sh`) can reproducibly create all real data subsets
- [x] **TDAT-07**: Expected output files exist for regression testing (pre-validated pipeline output)
- [x] **TDAT-08**: VCF list files and test configs exist for Snakemake batch testing

### Testing — Framework

- [x] **TFWK-01**: BATS smoke tests verify --help exits 0, --version prints version, missing args exits non-zero
- [x] **TFWK-02**: BATS filter tests verify each GATK filter correctly tags matching variants
- [x] **TFWK-03**: BATS filter tests verify each Freebayes filter correctly tags matching variants
- [x] **TFWK-04**: BATS integration test runs full pipeline on real data and produces valid output
- [x] **TFWK-05**: BATS regression tests compare pipeline output against expected output files
- [x] **TFWK-06**: BATS test verifies -g/--genome flag works (issue #13 — close after verification)
- [x] **TFWK-07**: BATS edge case tests verify empty VCF, no filters, no BED files don't crash

### Refactoring

- [x] **REFR-01**: Logging module (`lib/logging.sh`) provides log_msg, debug_msg, error handling functions
- [x] **REFR-02**: CLI module (`lib/cli.sh`) handles argument parsing, validation, show_help
- [x] **REFR-03**: Genome module (`lib/genome.sh`) handles genome file creation and validation
- [x] **REFR-04**: BED module (`lib/bed.sh`) handles normalize_bed, merge, intersect, compress, index
- [x] **REFR-05**: Annotate module (`lib/annotate.sh`) handles VCF annotation with BED regions
- [x] **REFR-06**: Normalize module (`lib/normalize.sh`) wraps bcftools norm operations
- [x] **REFR-07**: Filter module (`lib/filter.sh`) handles filter pipeline construction and execution
- [x] **REFR-08**: Stats module (`lib/stats.sh`) handles stats generation and plotting
- [x] **REFR-09**: Main script (`hardnormly.sh`) sources lib/ modules and orchestrates pipeline steps
- [x] **REFR-10**: Filter pipeline uses safer construction (no `eval` on dynamically built strings)
- [x] **REFR-11**: All external commands use consistent error handling pattern (`run_cmd` or equivalent)
- [x] **REFR-12**: Duplicate cleanup logic eliminated (single trap-based pattern)
- [x] **REFR-13**: Unit tests exist for each lib/ module's key functions

### Features

- [ ] **FEAT-01**: Plot-vcfstats errors are non-fatal (logged but don't exit the pipeline)
- [ ] **FEAT-02**: Subcommand dispatcher routes to generate-inclusion-bed, generate-exclusion-bed, or run-pipeline
- [ ] **FEAT-03**: `generate-inclusion-bed` subcommand produces a merged include BED from inputs
- [ ] **FEAT-04**: `generate-exclusion-bed` subcommand produces a merged exclusion BED from inputs
- [ ] **FEAT-05**: Helper script generates combined exclusion BED files for hg19 and hg38 from public sources
- [ ] **FEAT-06**: `--caller` flag auto-selects the matching filter file (gatk → gatk_filters.txt, freebayes → freebayes_filters.txt)
- [ ] **FEAT-07**: Running without a subcommand defaults to `run-pipeline` (backward compatible)
- [ ] **FEAT-08**: `--strip-annotations` flag removes specified INFO fields (e.g., INFO/CSQ) via `bcftools annotate -x` before filtering

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

- **FEAT-10**: `--plot-backend` option to choose between tectonic and pdflatex
- **FEAT-11**: `--help-filters` subcommand prints filter file format in detail

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
| INFR-01 | Phase 1 | Complete |
| INFR-02 | Phase 1 | Complete |
| INFR-03 | Phase 1 | Complete |
| INFR-04 | Phase 1 | Complete |
| INFR-05 | Phase 1 | Complete |
| INFR-06 | Phase 1 | Complete |
| TDAT-01 | Phase 2 | Complete |
| TDAT-02 | Phase 2 | Complete |
| TDAT-03 | Phase 2 | Complete |
| TDAT-04 | Phase 2 | Complete |
| TDAT-05 | Phase 2 | Complete |
| TDAT-06 | Phase 2 | Complete |
| TDAT-07 | Phase 2 | Complete |
| TDAT-08 | Phase 2 | Complete |
| TFWK-01 | Phase 3 | Complete |
| TFWK-02 | Phase 3 | Complete |
| TFWK-03 | Phase 3 | Complete |
| TFWK-04 | Phase 3 | Complete |
| TFWK-05 | Phase 3 | Complete |
| TFWK-06 | Phase 3 | Complete |
| TFWK-07 | Phase 3 | Complete |
| REFR-01 | Phase 4 | Complete |
| REFR-02 | Phase 4 | Complete |
| REFR-03 | Phase 4 | Complete |
| REFR-04 | Phase 4 | Complete |
| REFR-05 | Phase 4 | Complete |
| REFR-06 | Phase 4 | Complete |
| REFR-07 | Phase 4 | Complete |
| REFR-08 | Phase 4 | Complete |
| REFR-09 | Phase 4 | Complete |
| REFR-10 | Phase 4 | Complete |
| REFR-11 | Phase 4 | Complete |
| REFR-12 | Phase 4 | Complete |
| REFR-13 | Phase 4 | Complete |
| FEAT-01 | Phase 5 | Pending |
| FEAT-02 | Phase 5 | Pending |
| FEAT-03 | Phase 5 | Pending |
| FEAT-04 | Phase 5 | Pending |
| FEAT-05 | Phase 5 | Pending |
| FEAT-06 | Phase 5 | Pending |
| FEAT-07 | Phase 5 | Pending |
| FEAT-08 | Phase 5 | Pending |
| DOCS-01 | Phase 5 | Pending |
| DOCS-02 | Phase 5 | Pending |
| DOCS-03 | Phase 5 | Pending |

**Coverage:**
- v1 requirements: 44 total (note: count in header was 43 — INFR through DOCS per category counts)
- Mapped to phases: 44
- Unmapped: 0

---
*Requirements defined: 2026-02-18*
*Last updated: 2026-02-18 after roadmap creation — all requirements mapped*
