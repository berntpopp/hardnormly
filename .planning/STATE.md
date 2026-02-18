# State: hardnormly

## Project Reference

See: .planning/PROJECT.md (updated 2026-02-18)

**Core value:** Reliably normalize and filter VCF files for clean variant calls
**Current focus:** Phase 3 — Test Framework

## Current Position

Phase: 3 of 5 (Test Framework)
Plan: 1 of 3 in current phase
Status: In progress
Last activity: 2026-02-18 — Completed 03-01-PLAN.md (BATS infrastructure and smoke tests)

Progress: [████░░░░░░] 40% (8/20 plans complete)

## Performance Metrics

**Velocity:**
- Total plans completed: 8
- Average duration: ~26 minutes
- Total execution time: ~205 minutes

**By Phase:**

| Phase | Plans | Total | Avg/Plan |
|-------|-------|-------|----------|
| 01-infrastructure | 3/3 ✓ | ~65 min | ~22 min |
| 02-test-data | 4/4 ✓ | ~125 min | ~31 min |
| 03-test-framework | 1/3 | ~25 min | ~25 min |

**Recent Trend:**
- Last 5 plans: 02-02 (45 min), 02-03 (45 min), 02-04 (10 min), 03-01 (25 min)
- Trend: stable ~10-45 min

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
- hardnormly.sh had 3 bugs in --auto-index and exclude BED annotation (fixed in 02-04): exclude annotation "exclude"→"1", --write-index=tbi flag syntax (was -W tbi)
- Test configs include ref.genome_file pointing to tests/data/hg19_chr22.genome to avoid UCSC MySQL queries during tests
- Expected output baseline records commit at time of generation (not task commit) for accurate provenance
- BATS installed via git submodules (not bats-action) — uses load not bats_load_library, requires submodules:recursive in CI checkout
- _get_filter helper takes 3 args (vcf, chrom, pos) for forward-compatibility with hg38 chr-prefixed VCFs
- .gitattributes eol=lf added for *.sh *.bash *.bats — fixes CRLF corruption on Windows with core.autocrlf=true
- CI test job uses apt-get bcftools/bedtools (not conda) — sufficient for BATS tests

### Pending Todos

None yet.

### Blockers/Concerns

- `make help` fails on Windows dev environment (sh.exe @echo issue) — works on Linux CI. Not blocking.
- New .sh files in future phases must be manually added to Makefile SH_FILES and CI lint steps.
- hardnormly.sh 3 bug fixes in 02-04 may affect pre-existing test expectations — BATS tests in Phase 3 should use the fixed expected outputs in tests/data/expected/
- Plan 03-02/03-03 must call _require_tools in setup() since filter tests need bcftools/bedtools
- Windows dev: BATS submodule files need dos2unix after clone (CRLF from autocrlf=true) — .gitattributes eol=lf prevents this going forward but existing clones need manual fix

## Session Continuity

Last session: 2026-02-18
Stopped at: Completed 03-01-PLAN.md — BATS infrastructure and smoke tests (10 files created/modified)
Resume file: None
