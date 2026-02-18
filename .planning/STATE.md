# State: hardnormly

## Project Reference

See: .planning/PROJECT.md (updated 2026-02-18)

**Core value:** Reliably normalize and filter VCF files for clean variant calls
**Current focus:** Phase 4 — Refactoring

## Current Position

Phase: 4 of 5 (Refactoring) — In progress
Plan: 1 of 5 in phase 04 (done)
Status: In progress
Last activity: 2026-02-18 — Completed 04-01-PLAN.md (logging and CLI extraction into lib/)

Progress: [███████░░░] 55% (11/20 plans complete)

## Performance Metrics

**Velocity:**
- Total plans completed: 11
- Average duration: ~21 minutes
- Total execution time: ~237 minutes

**By Phase:**

| Phase | Plans | Total | Avg/Plan |
|-------|-------|-------|----------|
| 01-infrastructure | 3/3 ✓ | ~65 min | ~22 min |
| 02-test-data | 4/4 ✓ | ~125 min | ~31 min |
| 03-test-framework | 3/3 ✓ | ~43 min | ~14 min |
| 04-refactoring | 1/5 | ~4 min | ~4 min |

**Recent Trend:**
- Last 5 plans: 03-01 (25 min), 03-02 (3 min), 03-03 (25 min), 04-01 (4 min)
- Trend: refactoring plans fast (pure code extraction, no new logic)

*Updated after each plan completion*

## Accumulated Context

### Decisions

- Focus on hardnormly.sh main script, not Snakemake workflow (Snakemake needs cluster testing)
- Replace real test data with synthetic + public data (no clinical data in repo)
- Issue #13 (genome file option) already implemented — VERIFIED CLOSED by TFWK-06 tests
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
- BATS installed via setup script (tests/setup_bats.sh) — clones pinned versions with core.autocrlf=false, gitignored dirs, CI runs setup script before tests
- _get_filter helper takes 3 args (vcf, chrom, pos) for forward-compatibility with hg38 chr-prefixed VCFs
- .gitattributes eol=lf added for *.sh *.bash *.bats — fixes CRLF corruption on Windows with core.autocrlf=true
- CI test job uses apt-get bcftools/bedtools (not conda) — sufficient for BATS tests
- Filter test assertions use exact string equality — catches filter ordering regressions as well as tag presence
- setup_file() pattern: pipeline runs once per bats file, all @test blocks query BATS_FILE_TMPDIR output
- Integration tests use setup() not setup_file() — each test is an independent pipeline run
- Real-data integration tests require ref/hs37d5.fa (not mini_ref.fa) — 1000G variants at 16M+ exceed mini_ref coverage
- HAVE_FULL_REF guard allows TFWK-04 tests to skip when ref/hs37d5.fa absent (CI environments)
- Regression comparison uses bcftools query CHROM/POS/FILTER (not binary diff) — avoids header timestamp noise
- 1000G VCFs have only GT FORMAT field — GATK/Freebayes filters not applicable; integration tests run without --filters-file
- Include guard pattern: [[ -n "${_LIB_X_LOADED:-}" ]] && return 0; readonly _LIB_X_LOADED=1 (all lib/ modules)
- parse_args sets caller-scope globals directly (not namerefs) — simpler for 20+ variables, only called once
- parse_filter_args uses nameref (local -n) for caller-controlled array (testability)
- Boolean comparisons changed from bare $var to [[ "$var" == "true" ]] in hardnormly.sh for string variables
- run_cmd uses ${_TMP_DIR:+${_TMP_DIR}/} expansion — falls back to /tmp before set_tmp_dir is called
- ERR trap kept in hardnormly.sh for transition period (will be removed in 04-04 after all commands wrapped in run_cmd)
- lib/ modules sourced via _SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)" for robust path resolution
- New .sh files added to Makefile SH_FILES needed: lib/logging.sh, lib/cli.sh (blocker for next plans too)

### Pending Todos

- Add lib/logging.sh and lib/cli.sh (and future lib/ files) to Makefile SH_FILES and CI lint steps

### Blockers/Concerns

- `make help` fails on Windows dev environment (sh.exe @echo issue) — works on Linux CI. Not blocking.
- New .sh files in future phases must be manually added to Makefile SH_FILES and CI lint steps.
- Windows dev: BATS setup script uses -c core.autocrlf=false to prevent CRLF corruption on clone
- TFWK-04 real-data tests require ref/hs37d5.fa — CI must either provide this file or accept that tests skip
- run_cmd_with_retry not yet wired to UCSC MySQL call in hardnormly.sh — to be handled in 04-02 or 04-03 (genome.sh extraction)

## Session Continuity

Last session: 2026-02-18T19:27:45Z
Stopped at: Completed 04-01-PLAN.md — logging and CLI extraction into lib/logging.sh and lib/cli.sh
Resume file: None
