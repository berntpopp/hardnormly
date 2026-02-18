# State: hardnormly

## Project Reference

See: .planning/PROJECT.md (updated 2026-02-18)

**Core value:** Reliably normalize and filter VCF files for clean variant calls
**Current focus:** Phase 2 — Test Data

## Current Position

Phase: 2 of 5 (Test Data)
Plan: 3 of 4 in current phase
Status: In progress
Last activity: 2026-02-18 — Completed 02-03-PLAN.md (utility VCFs and real data subsets)

Progress: [███░░░░░░░] 30% (6/20 plans complete)

## Performance Metrics

**Velocity:**
- Total plans completed: 6
- Average duration: ~28 minutes
- Total execution time: ~170 minutes

**By Phase:**

| Phase | Plans | Total | Avg/Plan |
|-------|-------|-------|----------|
| 01-infrastructure | 3/3 ✓ | ~65 min | ~22 min |
| 02-test-data | 3/4 | ~115 min | ~38 min |

**Recent Trend:**
- Last 5 plans: 02-01 (25 min), 01-02 (25 min), 01-03 (25 min), 02-02 (45 min), 02-03 (45 min)
- Trend: stable ~25-45 min

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
- chr22:1-100001 extracted from hs37d5.fa as mini_ref.fa for synthetic VCF testing (all-N telomeric region)
- bioinformatics tools (samtools, bcftools, bedtools) available via WSL on Windows (/home/bernt/miniconda3/envs/hardnormly/bin)
- scripts/run_snakemake.sh excluded from CI lint and Makefile SH_FILES (Snakemake out of scope)
- hooks stored in .githooks/ (tracked in git, activated per-clone via setup-hooks.sh)
- Makefile omits SHELL := override — Windows GNU Make (Windows32) doesn't resolve /bin/bash
- Synthetic VCFs use REF=N throughout (mini_ref.fa is all-N telomeric sequence; bcftools TYPE() correctly classifies N>A as SNP, N>NA as INDEL)
- Freebayes AO/SAF/SAR/RPR/RPL declared as INFO fields (site-level Type=Integer) to match filter expression INFO/SAF==0 syntax
- NA19247 (YRI/AFR) used instead of NA19240 — NA19240 not present in 1000G Phase 3 chr22 dataset
- WSL tool invocations: use env -i with explicit paths to avoid MSYS path rewriting in Git Bash environment
- 1000G samples acquired: NA12878 (CEU/EUR, 20 variants), NA19247 (YRI/AFR, 30 variants), HG00096 (GBR/EUR, 17 variants)

### Pending Todos

None yet.

### Blockers/Concerns

- `make help` fails on Windows dev environment (sh.exe @echo issue) — works on Linux CI. Not blocking.
- New .sh files in future phases must be manually added to Makefile SH_FILES and CI lint steps.
- Plan 04 BATS tests should reference `1kg_NA19247` not `1kg_NA19240` (sample substitution from 02-03)

## Session Continuity

Last session: 2026-02-18
Stopped at: Completed 02-03-PLAN.md — utility VCFs and real data subsets (22 files committed)
Resume file: None
