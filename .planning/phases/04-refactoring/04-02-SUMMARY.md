---
phase: 04-refactoring
plan: 02
subsystem: pipeline
tags: [bash, bcftools, bedtools, bed-processing, genome-file, vcf-annotation, lib-modules]

requires:
  - phase: 04-01
    provides: lib/logging.sh and lib/cli.sh modules, module sourcing pattern, shellcheck/shfmt baseline

provides:
  - lib/genome.sh with create_genome_file (UCSC MySQL retry logic)
  - lib/bed.sh with normalize_bed, merge_include_beds, merge_exclude_beds, compress_index_bed, create_header_file
  - lib/annotate.sh with annotate_vcf_with_regions
  - hardnormly.sh Steps 1-4 delegated to lib modules
  - Makefile SH_FILES updated to include all lib/ modules

affects: [04-03, 04-04, 04-05]

tech-stack:
  added: []
  patterns:
    - "SC2311 fix: use $(set -e; fn) multiline form for command substitution that must preserve set -e"
    - "SC2310 fix: avoid || for function calls; rely on set -Eeuo pipefail to propagate failures"
    - "Lib module functions use return 1 not exit 1 (orchestrator handles exit)"
    - "annotate_vcf_with_regions uses run_cmd for bcftools annotate (error reporting included)"

key-files:
  created:
    - lib/genome.sh
    - lib/bed.sh
    - lib/annotate.sh
  modified:
    - hardnormly.sh
    - Makefile

key-decisions:
  - "SC2310/SC2311: removed || error handling from annotate_vcf_with_regions calls — rely on set -Eeuo pipefail to propagate failures from run_cmd"
  - "genome_file command substitution uses $(set -e; create_genome_file ...) multiline form per shfmt/shellcheck requirements"
  - "Makefile SH_FILES now includes all lib/ modules (lib/logging.sh, lib/cli.sh, lib/genome.sh, lib/bed.sh, lib/annotate.sh)"

patterns-established:
  - "Function calls that could fail: avoid || condition; let set -e propagate the failure"
  - "Command substitution preserving set -e: $(\\n\\tset -e\\n\\tcmd\\n) multiline form"

duration: 5min
completed: 2026-02-18
---

# Phase 4 Plan 02: Genome/BED/Annotate Extraction Summary

**Three lib modules extract ~100 lines from hardnormly.sh: lib/genome.sh (UCSC MySQL retry), lib/bed.sh (BED normalize/merge/compress/header), lib/annotate.sh (bcftools annotate via run_cmd)**

## Performance

- **Duration:** 5 min
- **Started:** 2026-02-18T19:31:15Z
- **Completed:** 2026-02-18T19:36:25Z
- **Tasks:** 2 (+1 deviation fix)
- **Files modified:** 5

## Accomplishments

- lib/genome.sh provides `create_genome_file` with 3-attempt retry loop replacing inline MySQL query
- lib/bed.sh provides 5 functions covering all BED normalization, merging, compression, and header creation
- lib/annotate.sh provides `annotate_vcf_with_regions` using `run_cmd` for proper error reporting
- hardnormly.sh Steps 1-4 reduced from ~90 inline lines to clean module function calls
- All 55 Phase 3 BATS tests pass without modification
- Makefile SH_FILES updated so `make lint` covers all lib modules

## Task Commits

Each task was committed atomically:

1. **Task 1: Create lib/genome.sh, lib/bed.sh, and lib/annotate.sh** - `906d731` (feat)
2. **Task 2: Update hardnormly.sh to delegate Steps 1-4** - `6eef3b7` (refactor)
3. **Deviation: Add lib/ modules to Makefile SH_FILES** - `d39088c` (chore)

**Plan metadata:** (docs commit follows)

## Files Created/Modified

- `lib/genome.sh` - `create_genome_file` with 3-attempt UCSC MySQL retry loop
- `lib/bed.sh` - `normalize_bed`, `merge_include_beds`, `merge_exclude_beds`, `compress_index_bed`, `create_header_file`
- `lib/annotate.sh` - `annotate_vcf_with_regions` wrapping bcftools annotate via run_cmd
- `hardnormly.sh` - Sources 3 new modules; Steps 1-4 replaced with module calls; inline `normalize_bed()` removed
- `Makefile` - All 5 lib/ modules added to SH_FILES

## Decisions Made

- SC2310/SC2311 shellcheck warnings for annotate_vcf_with_regions: removed `||` error handling and rely on `set -Eeuo pipefail` to propagate failure from `run_cmd`'s `return "$exit_code"`. Error messages already handled inside `run_cmd`.
- genome_file command substitution: uses shfmt-required multiline form `$(\n\tset -e\n\tcreate_genome_file ...\n)` to satisfy SC2311 (preserve set -e inside command substitution).
- Makefile SH_FILES: resolved STATE.md pending todo — all lib/ modules now covered by `make lint`.

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 2 - Missing Critical] Added lib/ modules to Makefile SH_FILES**

- **Found during:** Post-task verification
- **Issue:** STATE.md noted "New .sh files added to Makefile SH_FILES needed" as a pending todo; lib/genome.sh, lib/bed.sh, lib/annotate.sh were not in SH_FILES so `make lint` would not check them
- **Fix:** Added all 5 lib/ modules to SH_FILES in Makefile; verified `make lint` passes
- **Files modified:** Makefile
- **Verification:** `make lint` exits 0 with all lib files checked
- **Committed in:** `d39088c`

---

**Total deviations:** 1 auto-fixed (Rule 2 - missing lint coverage)
**Impact on plan:** Lint coverage fix required for correctness of CI pipeline. No scope creep.

## Issues Encountered

- shellcheck SC2310/SC2311 warnings triggered by `check-set-e-suppressed` optional check in `.shellcheckrc`. Required pattern adjustment: removed `||` from annotate calls (rely on set -e propagation) and used `$(set -e; fn)` form for genome file command substitution. Both patterns now documented as established conventions.

## Next Phase Readiness

- lib/genome.sh, lib/bed.sh, lib/annotate.sh ready for 04-03 (normalization extraction)
- 5 source lines in hardnormly.sh (logging, cli, genome, bed, annotate); 3 more planned (norm, filter, stats)
- All BATS tests passing, `make lint` clean

---
*Phase: 04-refactoring*
*Completed: 2026-02-18*
