# State: hardnormly

## Project Reference

See: .planning/PROJECT.md (updated 2026-02-18)

**Core value:** Reliably normalize and filter VCF files for clean variant calls
**Current focus:** Phase 1 — Infrastructure

## Current Position

Phase: 1 of 5 (Infrastructure)
Plan: 3 of 4 in current phase
Status: In progress
Last activity: 2026-02-18 — Completed 01-03-PLAN.md (pre-commit hook, CI workflow, Makefile)

Progress: [███░░░░░░░] 20% (4/20 plans complete)

## Performance Metrics

**Velocity:**
- Total plans completed: 4
- Average duration: ~20 minutes
- Total execution time: ~80 minutes

**By Phase:**

| Phase | Plans | Total | Avg/Plan |
|-------|-------|-------|----------|
| 01-infrastructure | 3/4 | ~65 min | ~22 min |
| 02-test-data | 1/4 | ~25 min | ~25 min |

**Recent Trend:**
- Last 5 plans: 01-01 (15 min), 02-01 (25 min), 01-02 (25 min), 01-03 (25 min)
- Trend: stable ~25 min

*Updated after each plan completion*

## Accumulated Context

### Decisions

- Focus on hardnormly.sh main script, not Snakemake workflow (Snakemake needs cluster testing)
- Replace real test data with synthetic + public data (no clinical data in repo)
- Issue #13 (genome file option) already implemented — verify with TFWK-06 test and close
- BATS chosen for testing (bash-native, well-supported)
- lib/ modularization pattern: 8 focused modules sourced by main script
- SC2294 (eval) temporary disable in 01-01 — RESOLVED in 01-02 with array-based pipeline (no eval)
- shfmt flags chosen: -i 0 (tabs) -bn (binary ops at line start) -ci (case body indent)
- cleanup_handler uses trap - EXIT guard before exit to prevent recursive EXIT trap invocation
- filter_stages array uses | delimiter (pipe cannot appear in bcftools filter expressions)
- shfmt requires || { multiline } blocks — inline single-line style rejected by formatter
- Gitignore negations must come AFTER the pattern they override (git last-match-wins rule)
- chr22:1-100001 extracted from hs37d5.fa as mini_ref.fa for synthetic VCF testing
- bioinformatics tools (samtools, bcftools, bedtools) available via WSL on Windows
- scripts/run_snakemake.sh excluded from CI lint and Makefile SH_FILES (Snakemake out of scope)
- hooks stored in .githooks/ (tracked in git, activated per-clone via setup-hooks.sh)
- Makefile omits SHELL := override — Windows GNU Make (Windows32) doesn't resolve /bin/bash

### Pending Todos

None yet.

### Blockers/Concerns

- `make help` fails on Windows dev environment (sh.exe @echo issue) — works on Linux CI. Not blocking.
- New .sh files in future phases must be manually added to Makefile SH_FILES and CI lint steps.

## Session Continuity

Last session: 2026-02-18
Stopped at: Completed 01-03-PLAN.md — pre-commit hook, GitHub Actions CI, Makefile
Resume file: None
