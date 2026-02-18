---
phase: 05-features-docs
plan: 01
subsystem: pipeline
tags: [bash, bcftools, shellcheck, shfmt, cli, vcf, filtering]

# Dependency graph
requires:
  - phase: 04-refactoring
    provides: "8 lib/ modules; hardnormly.sh clean orchestrator; lib/cli.sh parse_args; lib/annotate.sh annotate_vcf_with_regions"
provides:
  - "--caller flag: auto-selects gatk or freebayes default filter file"
  - "--strip-annotations flag: removes INFO fields before normalization"
  - "strip_vcf_annotations function in lib/annotate.sh"
  - "Non-fatal plot-vcfstats handling with if ! + SC2310 disable"
affects: [05-03-subcommand-dispatcher, testing]

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "Non-fatal function call: # shellcheck disable=SC2310 then if ! fn; then _rc=1; fi pattern (check-set-e-suppressed fires on if ! under .shellcheckrc)"
    - "--caller resolution happens in hardnormly.sh (has _SCRIPT_DIR), not in lib/cli.sh (parse_args only stores raw value)"
    - "Python script files used for in-place substitution when tab-exact string replacement is needed in WSL"

key-files:
  created: []
  modified:
    - lib/cli.sh
    - lib/annotate.sh
    - hardnormly.sh

key-decisions:
  - "SC2310 disable required: check-set-e-suppressed in .shellcheckrc fires on if ! fn even without ||; # shellcheck disable=SC2310 used on the if ! line"
  - "--caller path resolution placed in hardnormly.sh (after set_log_file) because _SCRIPT_DIR is not available in lib/cli.sh"
  - "strip_vcf_annotations inserted as step 4.5 between exclusion annotation and normalization"
  - "if ! pattern used for non-fatal plot: _plot_rc captured separately to avoid re-triggering SC2310 on compound condition"

patterns-established:
  - "Non-fatal function: # shellcheck disable=SC2310 + if ! fn; then _rc=1; fi + if [[ _rc -ne 0 ]]; then log; fi"

# Metrics
duration: 7min
completed: 2026-02-18
---

# Phase 5 Plan 01: New Pipeline Flags (--caller, --strip-annotations, non-fatal plot) Summary

**--caller auto-selects gatk/freebayes filter file; --strip-annotations removes INFO fields before normalization; plot-vcfstats failures non-fatal via shellcheck-safe if ! pattern**

## Performance

- **Duration:** 7 min
- **Started:** 2026-02-18T20:26:33Z
- **Completed:** 2026-02-18T20:33:38Z
- **Tasks:** 3
- **Files modified:** 3

## Accomplishments

- parse_args in lib/cli.sh handles --caller and --strip-annotations flags
- lib/annotate.sh has strip_vcf_annotations (bcftools annotate -x)
- hardnormly.sh: --caller resolves to gatk/freebayes filter file after logging config; --filters-file wins with warning; unknown --caller exits non-zero
- Step 4.5 in pipeline applies strip_vcf_annotations before normalization when --strip-annotations is set
- plot-vcfstats failure is caught and logged as warning; pipeline continues

## Task Commits

Each task was committed atomically:

1. **Task 1: Add --caller and --strip-annotations parsing to lib/cli.sh** - `150c75e` (feat)
2. **Task 2: Add strip_vcf_annotations to lib/annotate.sh** - `91dc2cc` (feat)
3. **Task 3: Wire non-fatal plot, --caller resolution, and --strip-annotations into hardnormly.sh** - `e7eb7f5` (feat)

## Files Created/Modified

- `lib/cli.sh` - Added --caller and --strip-annotations cases to parse_args; updated SC2034 disable comment
- `lib/annotate.sh` - Added strip_vcf_annotations function; updated Provides: comment
- `hardnormly.sh` - Added caller/strip_annotations defaults; --caller resolution block; step 4.5; non-fatal plot pattern

## Decisions Made

- **SC2310 on `if !`:** The project's `.shellcheckrc` enables `check-set-e-suppressed`, which causes SC2310 to fire on `if ! fn` (not just `||`). Used `# shellcheck disable=SC2310` inline before the `if !` line. The plan said "if ! avoids SC2310" but that was incorrect under check-set-e-suppressed.
- **--caller resolution in hardnormly.sh, not cli.sh:** `_SCRIPT_DIR` (needed for absolute path to `defaults/`) is only defined in hardnormly.sh. parse_args stores the raw value; resolution happens after `set_log_file`/`set_debug` so `log_msg` is available for the warning.
- **Step ordering:** --caller resolution block placed between `set_log_file`/`set_debug` and `validate_args` to enable `log_msg` for the conflict warning.

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 1 - Bug] SC2310 fires on `if !` under check-set-e-suppressed**

- **Found during:** Task 3 (wiring hardnormly.sh)
- **Issue:** Plan specified `if ! plot_stats_output ...` as SC2310-safe, but `.shellcheckrc` enables `check-set-e-suppressed` which triggers SC2310 even for `if !` on function calls
- **Fix:** Added `# shellcheck disable=SC2310` on the `if !` line; kept the `if !` pattern (correct intent) with explicit disable
- **Files modified:** hardnormly.sh
- **Verification:** `shellcheck hardnormly.sh` exits 0
- **Committed in:** e7eb7f5 (Task 3 commit)

---

**Total deviations:** 1 auto-fixed (1 bug — shellcheck SC2310 under check-set-e-suppressed)
**Impact on plan:** Minor — code intent preserved, disable comment added for transparency.

## Issues Encountered

- Python heredoc dollar-sign expansion: `$2`/`$1` in heredoc body were expanded to empty strings even with single-quoted PYEOF delimiter when running via `wsl bash -c "..."`. Worked around by writing a Python file to disk (`fix_cli.py`) and executing it separately.

## User Setup Required

None - no external service configuration required.

## Next Phase Readiness

- All three features (FEAT-01, FEAT-06, FEAT-08) ready for use by the pipeline
- 05-03 subcommand dispatcher plan can route to run-pipeline knowing --caller and --strip-annotations are already handled
- No blockers

---
*Phase: 05-features-docs*
*Completed: 2026-02-18*
