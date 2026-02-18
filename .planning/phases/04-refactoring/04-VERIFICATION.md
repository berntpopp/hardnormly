---
phase: 04-refactoring
verified: 2026-02-18T20:10:18Z
status: gaps_found
score: 4/5 must-haves verified
gaps:
  - truth: Every external command failure is caught by a consistent run_cmd or equivalent pattern
    status: failed
    reason: lib/bed.sh calls bedtools, bgzip, and tabix without run_cmd and without error capture or logging. lib/genome.sh implements its own retry loop without using run_cmd, despite the Requires comment claiming otherwise.
    artifacts:
      - path: lib/bed.sh
        issue: normalize_bed, merge_include_beds, merge_exclude_beds, compress_index_bed all call bedtools/bgzip/tabix with no run_cmd and no error handling. Failures propagate via pipefail with no log message naming the failed command.
      - path: lib/genome.sh
        issue: create_genome_file implements a manual retry loop calling mysql directly, not via run_cmd. MySQL failures do not produce the standard run_cmd error message. Requires comment falsely claims run_cmd_with_retry dependency.
    missing:
      - run_cmd or error-capture guards wrapping bedtools/bgzip/tabix calls in bed.sh
      - Either use run_cmd_with_retry in genome.sh OR document manual retry as intentional and correct the Requires comment
---

# Phase 4: Refactoring - Verification Report

**Phase Goal:** The script is split into focused lib/ modules with safe error handling -- tests still pass
**Verified:** 2026-02-18T20:10:18Z
**Status:** gaps_found
**Re-verification:** No -- initial verification

## Goal Achievement

### Observable Truths

| # | Criterion | Status | Evidence |
|---|-----------|--------|----------|
| 1 | All Phase 3 tests still pass after refactoring | VERIFIED | 112 tests run, 0 failures. 55 skipped (tools not available in env). |
| 2 | lib/ contains 8 focused modules | VERIFIED | annotate.sh, bed.sh, cli.sh, filter.sh, genome.sh, logging.sh, normalize.sh, stats.sh -- all present. |
| 3 | Filter pipeline runs without eval on dynamically built strings | VERIFIED | grep for eval in hardnormly.sh and lib/*.sh returns zero matches. Filter stages encoded as array elements, applied iteratively via IFS split. |
| 4 | Every external command failure caught by consistent run_cmd or equivalent | FAILED | bed.sh calls bedtools/bgzip/tabix with no run_cmd and no error capture. genome.sh has manual retry loop bypassing run_cmd. See Gaps section. |
| 5 | Unit tests exist for key functions in each lib/ module and pass | VERIFIED | tests/05-lib-logging.bats (152 lines), tests/06-lib-cli.bats (138 lines), tests/07-lib-modules.bats (300 lines) -- 57 unit tests total, 0 failures. |

**Score:** 4/5 truths verified

### Required Artifacts

| Artifact | Expected | Status | Details |
|----------|----------|--------|---------|
| lib/logging.sh | log_msg, debug_msg, error_msg, run_cmd, run_cmd_with_retry | VERIFIED | 98 lines, substantive, sourced line 12. run_cmd defined lines 53-76, run_cmd_with_retry lines 80-98. |
| lib/cli.sh | Argument parsing, validation, show_help | VERIFIED | 236 lines, substantive, sourced line 13. parse_args, validate_args, show_help, show_version, parse_filter_args all present. |
| lib/genome.sh | Genome file creation via UCSC MySQL | VERIFIED | 35 lines, sourced line 14. Manual retry loop present (not run_cmd_with_retry despite Requires comment). |
| lib/bed.sh | BED normalization, merge, compress, index | PARTIAL | 82 lines, sourced line 15. All 5 functions exist but external commands not wrapped in run_cmd -- REFR-11 gap. |
| lib/annotate.sh | VCF annotation with BED regions | VERIFIED | 24 lines, uses run_cmd correctly, sourced line 16. |
| lib/normalize.sh | bcftools norm wrapper | VERIFIED | 46 lines, uses manual stderr capture with error handling -- documented as intentional. |
| lib/filter.sh | Filter pipeline construction and execution | VERIFIED | 89 lines, consistent error return pattern throughout, sourced line 17. |
| lib/stats.sh | Stats generation and plotting | VERIFIED | 46 lines, uses run_cmd for bcftools stats, sourced line 18. |
| tests/05-lib-logging.bats | Unit tests for logging module | VERIFIED | 152 lines, 17 tests, all pass. |
| tests/06-lib-cli.bats | Unit tests for cli module | VERIFIED | 138 lines, 16 tests, all pass. |
| tests/07-lib-modules.bats | Unit tests for remaining modules | VERIFIED | 300 lines, 24 tests, all pass (tool-requiring tests skip gracefully). |

### Key Link Verification

| From | To | Via | Status | Details |
|------|----|-----|--------|---------|
| hardnormly.sh | lib/logging.sh | source line 12 | WIRED | Sourced first; all other modules depend on it. |
| hardnormly.sh | lib/cli.sh | source line 13 | WIRED | parse_args called line 59, validate_args called line 66. |
| hardnormly.sh | lib/genome.sh | source line 14 | WIRED | create_genome_file called line 83. |
| hardnormly.sh | lib/bed.sh | source line 15 | WIRED | normalize_bed, merge_include_beds, merge_exclude_beds, compress_index_bed, create_header_file all called in steps 2-3. |
| hardnormly.sh | lib/annotate.sh | source line 16 | WIRED | annotate_vcf_with_regions called lines 144, 158. |
| hardnormly.sh | lib/normalize.sh | source line 17 | WIRED | normalize_vcf called line 172. |
| hardnormly.sh | lib/filter.sh | source line 18 | WIRED | init_filter_pipeline, apply_filter_stages, write_filtered_output all called in step 6. |
| hardnormly.sh | lib/stats.sh | source line 19 | WIRED | generate_stats called line 204, plot_stats_output called line 210. |
| filter pipeline | no eval | array stages | WIRED | filter_stages array holds name-action-expr tuples; apply_filter_stages uses IFS split -- no eval. |
| cleanup | single EXIT trap | cleanup_handler | WIRED | Single trap cleanup_handler EXIT at line 56; no duplicate cleanup in any lib module. |

### Requirements Coverage

| Requirement | Status | Notes |
|-------------|--------|-------|
| REFR-01: lib/logging.sh | SATISFIED | Provides log_msg, debug_msg, error_msg, run_cmd, run_cmd_with_retry. |
| REFR-02: lib/cli.sh | SATISFIED | parse_args, validate_args, show_help, show_version, parse_filter_args. |
| REFR-03: lib/genome.sh | SATISFIED | create_genome_file with retry logic and error return. |
| REFR-04: lib/bed.sh | SATISFIED | All 5 BED functions exist and are callable. |
| REFR-05: lib/annotate.sh | SATISFIED | annotate_vcf_with_regions via run_cmd. |
| REFR-06: lib/normalize.sh | SATISFIED | normalize_vcf with intentional stderr capture pattern. |
| REFR-07: lib/filter.sh | SATISFIED | init_filter_pipeline, apply_filter_stages, write_filtered_output with consistent error handling. |
| REFR-08: lib/stats.sh | SATISFIED | generate_stats (run_cmd), plot_stats_output (non-fatal pattern). |
| REFR-09: hardnormly.sh sources lib/ and orchestrates | SATISFIED | 8 source lines, clean pipeline orchestration. |
| REFR-10: No eval on dynamically built filter strings | SATISFIED | Zero eval matches in hardnormly.sh and all lib/*.sh. |
| REFR-11: Consistent error handling for external commands | BLOCKED | bed.sh (bedtools/bgzip/tabix) and genome.sh (mysql) bypass run_cmd with no error capture. |
| REFR-12: Duplicate cleanup eliminated, single trap | SATISFIED | Single trap cleanup_handler EXIT in hardnormly.sh; no cleanup code in lib modules. |
| REFR-13: Unit tests for each lib/ module | SATISFIED | 05, 06, 07 bats files cover logging, cli, and all other modules. |

### Anti-Patterns Found

| File | Line | Pattern | Severity | Impact |
|------|------|---------|----------|--------|
| lib/bed.sh | 18, 34-36, 40, 55-57 | External commands (bedtools) with no error capture | Warning | Failures produce no error message -- only pipefail-triggered exit. |
| lib/bed.sh | 70-71 | bgzip/tabix with no error handling | Warning | bgzip or tabix failure produces no log entry. |
| lib/bed.sh | 5 | Comment says Requires run_cmd but run_cmd is never called | Info | Misleading comment -- bed.sh relies on pipefail propagation. |
| lib/genome.sh | 4 | Comment says Requires run_cmd_with_retry but manual retry loop used instead | Info | run_cmd_with_retry exists in logging.sh but genome.sh does not call it. |

### Human Verification Required

None -- all checks automated.

## Gaps Summary

One gap blocks full REFR-11 compliance.

### REFR-11 partial gap: bed.sh and genome.sh do not use run_cmd

lib/bed.sh: All five functions (normalize_bed, merge_include_beds, merge_exclude_beds,
compress_index_bed, create_header_file) call external tools (bedtools, bgzip, tabix) without
run_cmd or any error guard. When bedtools fails, the pipeline exits via set -Eeuo pipefail
without any log message identifying which command failed or showing its stderr output.

lib/genome.sh: create_genome_file implements its own retry loop (lines 20-32) calling mysql
directly. While it retries and returns 1 on exhaustion, it does not produce the standard
run_cmd-style error message, and the Requires: run_cmd_with_retry comment is incorrect.

Impact assessment: The pipeline does fail-fast on errors (pipefail + set -Eeuo enforces this).
However, operators receive no diagnostic identifying which bedtools command failed or its stderr
output. This weakens the safe error handling claim versus modules like annotate.sh, filter.sh,
and stats.sh which do log failure details via run_cmd or explicit error handling blocks.

### All other criteria fully met

- 112 tests pass, 0 failures (criterion 1)
- 8 lib/ modules exist and are substantive (criterion 2)
- No eval anywhere in hardnormly.sh or lib/*.sh (criterion 3)
- 57 unit tests in 3 files, all pass (criterion 5)
- make lint passes cleanly (shellcheck + shfmt produce no output)
- Single EXIT trap for cleanup, no duplicate cleanup logic (REFR-12)
- All 8 modules sourced and wired to actual pipeline calls (REFR-09)

---

_Verified: 2026-02-18T20:10:18Z_
_Verifier: Claude (gsd-verifier)_
