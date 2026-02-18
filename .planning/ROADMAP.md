# Roadmap: hardnormly

## Milestones

- 🚧 **v0.7.0 — Code Quality, Testing & Features** - Phases 1-5 (in progress)

## Overview

Transform hardnormly from a 463-line monolithic bash script with no tests or CI into a
well-structured, tested, linted, and documented tool with new capabilities. The path runs
infrastructure first (linting, CI), then test data, then tests that make refactoring safe,
then modularization, then features that depend on the modular structure.

## Phases

**Phase Numbering:**
- Integer phases (1, 2, 3): Planned milestone work
- Decimal phases (2.1, 2.2): Urgent insertions (marked with INSERTED)

Decimal phases appear between their surrounding integers in numeric order.

- [x] **Phase 1: Infrastructure** - Linting, CI, strict mode — quality foundation that unblocks everything
- [ ] **Phase 2: Test Data** - Synthetic and real data subsets that make automated tests possible
- [ ] **Phase 3: Test Framework** - BATS smoke, filter, integration, and regression tests
- [ ] **Phase 4: Refactoring** - Modularize into lib/, replace eval, unify error handling
- [ ] **Phase 5: Features & Docs** - Subcommands, --caller flag, exclusion BED generation, help text

## Phase Details

### Phase 1: Infrastructure

**Goal**: The project enforces consistent code quality automatically — locally and in CI
**Depends on**: Nothing (first phase)
**Requirements**: INFR-01, INFR-02, INFR-03, INFR-04, INFR-05, INFR-06
**Plans:** 3 plans
**Success Criteria** (what must be TRUE):
  1. Running `shellcheck hardnormly.sh` exits 0 with no warnings
  2. Running `shfmt -d hardnormly.sh` exits 0 (no diffs)
  3. Pushing a commit triggers GitHub Actions and shows green checks for lint and test jobs
  4. Running `git commit` with a ShellCheck violation causes the pre-commit hook to reject it
  5. The script exits non-zero immediately when any piped command fails

Plans:
- [x] 01-01-PLAN.md — ShellCheck config, fix findings, apply shfmt formatting (INFR-01, INFR-02)
- [x] 01-02-PLAN.md — Strict mode, traps, exception handling, eval replacement (INFR-03)
- [x] 01-03-PLAN.md — Pre-commit hooks, GitHub Actions CI, Makefile (INFR-04, INFR-05, INFR-06)

### Phase 2: Test Data

**Goal**: All the data needed to run automated tests exists in the repo and can be regenerated
**Depends on**: Phase 1
**Requirements**: TDAT-01, TDAT-02, TDAT-03, TDAT-04, TDAT-05, TDAT-06, TDAT-07, TDAT-08
**Plans:** 4 plans
**Success Criteria** (what must be TRUE):
  1. `tests/data/` contains synthetic VCFs that include variants triggering every GATK and Freebayes filter
  2. A minimal reference FASTA and matching BED files exist for the chr22 test region
  3. Running `tests/generate_test_data.sh` from a clean checkout reproduces all real data subsets
  4. Pre-computed expected output VCFs exist that the regression tests can diff against

Plans:
- [ ] 02-01-PLAN.md — Fix gitignore, create directory structure, reference FASTA, BED and genome files (TDAT-04)
- [ ] 02-02-PLAN.md — Create synthetic GATK and Freebayes VCFs with all filter triggers (TDAT-01, TDAT-02)
- [ ] 02-03-PLAN.md — Utility VCFs, real data generation script, download real data subsets (TDAT-03, TDAT-05, TDAT-06)
- [ ] 02-04-PLAN.md — Generate expected output files and Snakemake test configs (TDAT-07, TDAT-08)

### Phase 3: Test Framework

**Goal**: Automated tests run against the existing script and catch regressions before any refactoring
**Depends on**: Phase 2
**Requirements**: TFWK-01, TFWK-02, TFWK-03, TFWK-04, TFWK-05, TFWK-06, TFWK-07
**Success Criteria** (what must be TRUE):
  1. Running `bats tests/` passes all tests against the current (pre-refactor) script
  2. Each GATK filter is verified to tag exactly the variants it should (and no others)
  3. Each Freebayes filter is verified to tag exactly the variants it should
  4. Running the full pipeline on real data produces a valid, non-empty VCF
  5. Passing an empty VCF, omitting filters, or omitting BED files does not crash the script
**Plans:** 3 plans

Plans:
- [ ] 03-01-PLAN.md — BATS infrastructure, shared helper, smoke tests, CI/Makefile (TFWK-01)
- [ ] 03-02-PLAN.md — GATK and Freebayes filter unit tests (TFWK-02, TFWK-03)
- [ ] 03-03-PLAN.md — Integration, regression, genome flag, and edge case tests (TFWK-04, TFWK-05, TFWK-06, TFWK-07)

### Phase 4: Refactoring

**Goal**: The script is split into focused lib/ modules with safe error handling — tests still pass
**Depends on**: Phase 3
**Requirements**: REFR-01, REFR-02, REFR-03, REFR-04, REFR-05, REFR-06, REFR-07, REFR-08,
                 REFR-09, REFR-10, REFR-11, REFR-12, REFR-13
**Success Criteria** (what must be TRUE):
  1. All Phase 3 tests still pass after refactoring (no regressions)
  2. `lib/` contains 8 focused modules (logging, cli, genome, bed, annotate, normalize, filter, stats)
  3. The filter pipeline runs without eval on dynamically built strings
  4. Every external command failure is caught by a consistent run_cmd or equivalent pattern
  5. Unit tests exist for key functions in each lib/ module and pass
**Plans**: TBD

Plans:
- [ ] 04-01: Create lib/logging.sh and lib/cli.sh, update main script to source them
- [ ] 04-02: Create lib/genome.sh, lib/bed.sh, lib/annotate.sh
- [ ] 04-03: Create lib/normalize.sh, lib/filter.sh (replace eval), lib/stats.sh
- [ ] 04-04: Finalize hardnormly.sh as orchestrator, fix cleanup/trap pattern
- [ ] 04-05: Write unit tests for each lib/ module

### Phase 5: Features & Docs

**Goal**: New user-facing capabilities work and are discoverable via help text
**Depends on**: Phase 4
**Requirements**: FEAT-01, FEAT-02, FEAT-03, FEAT-04, FEAT-05, FEAT-06, FEAT-07,
                 DOCS-01, DOCS-02, DOCS-03
**Success Criteria** (what must be TRUE):
  1. Running `hardnormly.sh run-pipeline [args]` and `hardnormly.sh [args]` produce identical output
  2. Running `hardnormly.sh generate-inclusion-bed` and `generate-exclusion-bed` produce valid merged BED files
  3. Running `hardnormly.sh --caller gatk` uses gatk_filters.txt without requiring --filters-file
  4. A plot-vcfstats failure is logged but the pipeline continues and exits 0
  5. `hardnormly.sh --help` shows usage examples and filter file format description
**Plans**: TBD

Plans:
- [ ] 05-01: Make plot-vcfstats errors non-fatal (FEAT-01)
- [ ] 05-02: Add subcommand dispatcher and generate-inclusion-bed, generate-exclusion-bed (FEAT-02 to FEAT-04, FEAT-07)
- [ ] 05-03: Add --caller flag and exclusion BED generation helper script (FEAT-05, FEAT-06)
- [ ] 05-04: Enhance --help with examples, filter format, BED source docs (DOCS-01, DOCS-02, DOCS-03)

## Progress

**Execution Order:**
Phases execute in numeric order: 1 → 2 → 3 → 4 → 5

| Phase | Plans Complete | Status | Completed |
|-------|----------------|--------|-----------|
| 1. Infrastructure | 3/3 | Complete | 2026-02-18 |
| 2. Test Data | 0/4 | Planned | - |
| 3. Test Framework | 0/3 | Planned | - |
| 4. Refactoring | 0/5 | Not started | - |
| 5. Features & Docs | 0/4 | Not started | - |
