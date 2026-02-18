---
phase: 01-infrastructure
plan: 02
subsystem: infra
tags: [bash, shellcheck, shfmt, strict-mode, error-handling, bcftools, filter-pipeline]

# Dependency graph
requires:
  - phase: 01-01
    provides: ShellCheck + shfmt baseline with zero findings on hardnormly.sh
provides:
  - set -Eeuo pipefail strict mode with -E errtrace for function ERR propagation
  - err_handler function reporting command, exit code, and line number on failure
  - cleanup_handler EXIT trap with recursive-entry guard (trap - EXIT)
  - Array-based sequential BCF filter pipeline replacing eval string-building
  - Zero shellcheck findings, zero shfmt findings after all changes
affects:
  - 01-03 (lib/ modularization — strict mode must propagate into sourced modules)
  - 01-04 (BATS tests — tests can now assert on err_handler stderr output)
  - 02-test-data (any test invocations must handle strict mode exits correctly)

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "err_handler: capture $? before local, print BASH_COMMAND and LINENO to stderr"
    - "cleanup_handler: trap - EXIT before exit to prevent recursive invocation"
    - "filter_stages array: encode name|action|expression, apply sequentially via BCF temp files"
    - "|| { multiline } over || { inline } for shfmt compliance"

key-files:
  created: []
  modified:
    - hardnormly.sh

key-decisions:
  - "Strict mode set -Eeuo pipefail placed immediately after version declaration (before default values), so variable setup itself is covered"
  - "cleanup_handler uses trap - EXIT guard before calling exit to prevent infinite recursion"
  - "filter_stages uses | as delimiter since | cannot appear in bcftools filter expressions"
  - "Both Task 1 and Task 2 committed atomically in one commit (81bff93) since changes are inseparable in the same file"
  - "shfmt requires || { multiline } blocks — inline single-line || { }; style was rejected"

patterns-established:
  - "ERR trap pattern: trap 'err_handler ${LINENO}' ERR with single quotes so LINENO expands at fire time"
  - "EXIT trap pattern: cleanup_handler resets EXIT trap before calling exit to prevent re-entry"
  - "Sequential BCF pipeline: initialize to filter_current.bcf, loop stages renaming filter_next.bcf to filter_current.bcf"
  - "Empty array loops: safe under set -u because arrays are declared with =() before use"

# Metrics
duration: 25min
completed: 2026-02-18
---

# Phase 1 Plan 02: Strict Mode and Eval-Free Filter Pipeline Summary

**`set -Eeuo pipefail` strict mode, ERR/EXIT traps with line-number error reporting, and array-based sequential BCF filter pipeline replacing the eval string builder — hardnormly.sh fails fast on errors and passes shellcheck + shfmt clean**

## Performance

- **Duration:** ~25 min
- **Started:** 2026-02-18T00:00:00Z (approx)
- **Completed:** 2026-02-18
- **Tasks:** 2 (implemented atomically)
- **Files modified:** 1 (hardnormly.sh)

## Accomplishments

- `set -Eeuo pipefail` now covers the entire script including functions (via `-E` errtrace), preventing silent failures in VCF processing pipelines
- ERR trap fires with the exact failed command, exit code, and line number — eliminates silent error debugging
- EXIT trap always cleans up `$tmp_dir` on both normal and error exits, with recursive-entry guard
- Filter pipeline is now eval-free: `filter_stages` array encodes each filter as `name|action|expression`, applied sequentially via BCF temp files
- All deprecated patterns removed: `eval`, `pipeline_cmd` string building, `norm_exit_code`/`plot_exit_code` capture patterns, SC2294 shellcheck disable, duplicate cleanup block

## Task Commits

Both tasks were implemented atomically in a single rewrite (changes are inseparable in the same file):

1. **Task 1: Add strict mode, traps, and fix exception handling patterns** - `81bff93` (feat)
2. **Task 2: Replace eval-based filter pipeline with array-based sequential approach** - `81bff93` (feat)

**Plan metadata:** (see docs commit below)

## Files Created/Modified

- `C:\development\hardnormly\hardnormly.sh` — Main pipeline script with strict mode, traps, safe exception patterns, and array-based filter pipeline

## Decisions Made

- **Single atomic commit**: Tasks 1 and 2 were implemented as a complete file rewrite rather than two incremental edits. The changes are logically separable but file-inseparable — committing them together in `81bff93` is the accurate representation.
- **`trap - EXIT` guard**: The `cleanup_handler` resets the EXIT trap before calling `exit "$exit_code"` to prevent infinite recursive invocation. This is a correctness fix not explicitly in the plan, discovered during implementation (Rule 1 — bug fix).
- **`filter_stages` `|` delimiter**: Pipe character chosen as field separator because it cannot appear in bcftools filter expressions (which use operators like `!=`, `=`, `<`, `>`).
- **shfmt multiline `|| {}`**: shfmt 3.x requires `|| {` blocks to be multi-line. All five inline `|| { ...; exit 1; }` patterns were expanded to 3-line blocks.

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 1 - Bug] Added `trap - EXIT` guard in cleanup_handler to prevent recursive invocation**

- **Found during:** Task 1 (cleanup_handler implementation)
- **Issue:** When `cleanup_handler` (the EXIT trap) calls `exit "$exit_code"`, bash fires the EXIT trap again, calling `cleanup_handler` recursively in an infinite loop. This is a real correctness bug that would hang the shell on any exit.
- **Fix:** Added `trap - EXIT` as the first statement in `cleanup_handler`, resetting the EXIT trap before performing cleanup or calling `exit`. Standard pattern for EXIT trap handlers.
- **Files modified:** hardnormly.sh
- **Verification:** Pattern visible at line 40: `trap - EXIT # Prevent recursive re-entry`
- **Committed in:** 81bff93 (combined task commit)

**2. [Rule 1 - Bug] Expanded inline `|| { }` blocks to multi-line format for shfmt compliance**

- **Found during:** Task 2 verification (shfmt -d check)
- **Issue:** Five inline `|| { log_msg "..."; exit 1; }` patterns failed shfmt formatting check — shfmt requires multi-line block syntax for `|| {` compound commands.
- **Fix:** Expanded all five occurrences to 3-line blocks with proper indentation.
- **Files modified:** hardnormly.sh
- **Verification:** `shfmt -d -i 0 -bn -ci hardnormly.sh` exits 0
- **Committed in:** 81bff93 (combined task commit)

---

**Total deviations:** 2 auto-fixed (both Rule 1 — bugs)
**Impact on plan:** Both fixes necessary for correctness and tool compliance. No scope creep.

## Issues Encountered

None — shellcheck and shfmt both pass cleanly after the two auto-fixes above.

## User Setup Required

None — no external service configuration required.

## Next Phase Readiness

- `hardnormly.sh` now has strict mode, proper traps, and eval-free filtering — ready for Plan 01-03 (lib/ modularization)
- Future sourced modules must work under `set -Eeuo pipefail`; the `err_handler` and `cleanup_handler` will propagate correctly via `-E` (errtrace)
- BATS tests (Plan 01-04) can now assert on stderr output from err_handler for negative test cases

---
*Phase: 01-infrastructure*
*Completed: 2026-02-18*
