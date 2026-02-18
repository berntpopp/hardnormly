# State: hardnormly

## Project Reference

See: .planning/PROJECT.md (updated 2026-02-18)

**Core value:** Reliably normalize and filter VCF files for clean variant calls
**Current focus:** Phase 4 — Refactoring

## Current Position

Phase: 4 of 5 (Refactoring) — In progress
Plan: 3 of 5 in phase 04 (done)
Status: In progress
Last activity: 2026-02-18 — Completed 04-03-PLAN.md (normalize/filter/stats extraction into lib/)

Progress: [████████░░] 65% (13/20 plans complete)

## Performance Metrics

**Velocity:**
- Total plans completed: 12
- Average duration: ~20 minutes
- Total execution time: ~242 minutes

**By Phase:**

| Phase | Plans | Total | Avg/Plan |
|-------|-------|-------|----------|
| 01-infrastructure | 3/3 ✓ | ~65 min | ~22 min |
| 02-test-data | 4/4 ✓ | ~125 min | ~31 min |
| 03-test-framework | 3/3 ✓ | ~43 min | ~14 min |
| 04-refactoring | 3/5 | ~14 min | ~4.7 min |

**Recent Trend:**
- Last 5 plans: 03-03 (25 min), 04-01 (4 min), 04-02 (5 min), 04-03 (5 min)
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
- SC2310/SC2311 pattern: avoid || for function calls in set -Eeuo pipefail scripts; rely on propagation; use $(set -e; fn) form for command substitutions
- Makefile SH_FILES now includes all 8 lib/ modules (lib/logging.sh, lib/cli.sh, lib/genome.sh, lib/bed.sh, lib/annotate.sh, lib/normalize.sh, lib/filter.sh, lib/stats.sh)
- create_genome_file uses manual retry loop (not run_cmd_with_retry) because it needs stdout capture; run_cmd_with_retry resolved for genome in 04-02 via retry loop in lib/genome.sh
- normalize_vcf does NOT use run_cmd — bcftools norm emits warnings on stderr even on success; manual mktemp capture used to preserve Warning/Lines logging
- apply_filter_stages takes stages as positional args after tmp_dir — clean varargs signature avoiding nameref complexity
- All 8 lib/ modules complete; hardnormly.sh is now a slim orchestrator (~227 lines)

### Pending Todos

- None (Makefile SH_FILES pending todo resolved in 04-02)

### Blockers/Concerns

- `make help` fails on Windows dev environment (sh.exe @echo issue) — works on Linux CI. Not blocking.
- Windows dev: BATS setup script uses -c core.autocrlf=false to prevent CRLF corruption on clone
- TFWK-04 real-data tests require ref/hs37d5.fa — CI must either provide this file or accept that tests skip

## Session Continuity

Last session: 2026-02-18T19:44:35Z
Stopped at: Completed 04-03-PLAN.md — normalize/filter/stats extraction into lib/normalize.sh, lib/filter.sh, lib/stats.sh
Resume file: None
