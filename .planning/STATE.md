# State: hardnormly

## Project Reference

See: .planning/PROJECT.md (updated 2026-02-18)

**Core value:** Reliably normalize and filter VCF files for clean variant calls
**Current focus:** Phase 2 — Test Data

## Current Position

Phase: 2 of 5 (Test Data)
Plan: 1 of 4 in current phase
Status: In progress
Last activity: 2026-02-18 — Completed 02-01-PLAN.md (test data foundation)

Progress: [█░░░░░░░░░] 10% (2/20 plans complete)

## Performance Metrics

**Velocity:**
- Total plans completed: 2
- Average duration: ~20 minutes
- Total execution time: ~40 minutes

**By Phase:**

| Phase | Plans | Total | Avg/Plan |
|-------|-------|-------|----------|
| 01-infrastructure | 1/4 | ~15 min | ~15 min |
| 02-test-data | 1/4 | ~25 min | ~25 min |

**Recent Trend:**
- Last 5 plans: 01-01 (15 min), 02-01 (25 min)
- Trend: --

*Updated after each plan completion*

## Accumulated Context

### Decisions

- Focus on hardnormly.sh main script, not Snakemake workflow (Snakemake needs cluster testing)
- Replace real test data with synthetic + public data (no clinical data in repo)
- Issue #13 (genome file option) already implemented — verify with TFWK-06 test and close
- BATS chosen for testing (bash-native, well-supported)
- lib/ modularization pattern: 8 focused modules sourced by main script
- SC2294 (eval) gets a documented temporary disable in 01-01; Plan 02 replaces eval with array-based pipeline
- shfmt flags chosen: -i 0 (tabs) -bn (binary ops at line start) -ci (case body indent)
- Gitignore negations must come AFTER the pattern they override (git last-match-wins rule)
- chr22:1-100001 extracted from hs37d5.fa as mini_ref.fa for synthetic VCF testing
- bioinformatics tools (samtools, bcftools, bedtools) available via WSL on Windows

### Pending Todos

None yet.

### Blockers/Concerns

None yet.

## Session Continuity

Last session: 2026-02-18
Stopped at: Completed 02-01-PLAN.md — test data foundation (gitignore, mini_ref.fa, BED files, genome file)
Resume file: None
