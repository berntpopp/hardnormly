---
phase: 04-refactoring
plan: 01
subsystem: infra
tags: [bash, modularization, logging, cli, shellcheck, shfmt, bats]

# Dependency graph
requires:
  - phase: 03-test-framework
    provides: 55 BATS tests validating hardnormly.sh behavior (regression safety net for extraction)
provides:
  - lib/logging.sh with log_msg, debug_msg, error_msg, run_cmd, run_cmd_with_retry, set_log_file, set_debug, set_tmp_dir
  - lib/cli.sh with show_help, show_version, parse_args, validate_args, parse_filter_args
  - hardnormly.sh updated to source both modules (no duplicate logging/CLI code)
  - Include guard pattern (_LIB_LOGGING_LOADED, _LIB_CLI_LOADED) established for all lib/ modules
affects: [04-02, 04-03, 04-04, 04-05]

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "Include guard: [[ -n \"${_LIB_X_LOADED:-}\" ]] && return 0; readonly _LIB_X_LOADED=1"
    - "Setter pattern for logging config (set_log_file, set_debug, set_tmp_dir)"
    - "parse_args sets caller-scope globals (standard bash — no local declarations)"
    - "parse_filter_args uses nameref (local -n) for caller-controlled array"
    - "SC2034/SC2154 disable comment on parse_args (variables used in caller scope)"
    - "_SCRIPT_DIR via BASH_SOURCE[0] for robust lib/ sourcing"

key-files:
  created:
    - lib/logging.sh
    - lib/cli.sh
  modified:
    - hardnormly.sh

key-decisions:
  - "parse_args sets caller-scope globals rather than using namerefs (simpler for 20+ variables, only called once)"
  - "SC2034/SC2154 suppressed at function level only (not file level) for parse_args"
  - "Boolean comparisons changed from bare $var to [[ \"$var\" == \"true\" ]] in hardnormly.sh for string variables"
  - "run_cmd uses _TMP_DIR:+ parameter expansion for fallback to /tmp when _TMP_DIR not set"
  - "ERR trap kept in hardnormly.sh for transition period (will be removed in 04-04)"

patterns-established:
  - "Include guard: [[ -n \"${_LIB_X_LOADED:-}\" ]] && return 0; readonly _LIB_X_LOADED=1"
  - "Source order: logging.sh before cli.sh (cli.sh calls log_msg)"
  - "Setter functions called immediately after parse_args (set_log_file, set_debug)"
  - "set_tmp_dir called after mkdir -p (tmp_dir confirmed to exist before registering)"

# Metrics
duration: 4min
completed: 2026-02-18
---

# Phase 4 Plan 1: Logging and CLI Extraction Summary

**Extracted logging (log_msg/debug_msg/error_msg/run_cmd) and CLI (parse_args/validate_args/parse_filter_args) into lib/logging.sh and lib/cli.sh with include guard pattern, all 55 BATS tests pass**

## Performance

- **Duration:** ~4 min
- **Started:** 2026-02-18T19:24:00Z
- **Completed:** 2026-02-18T19:27:45Z
- **Tasks:** 2
- **Files modified:** 3 (lib/logging.sh created, lib/cli.sh created, hardnormly.sh modified)

## Accomplishments

- Created `lib/logging.sh` with 8 functions: `set_log_file`, `set_debug`, `set_tmp_dir`, `log_msg`, `debug_msg`, `error_msg`, `run_cmd`, `run_cmd_with_retry`
- Created `lib/cli.sh` with 5 functions: `show_help`, `show_version`, `parse_args`, `validate_args`, `parse_filter_args`
- Updated `hardnormly.sh` to source both modules via `_SCRIPT_DIR` resolution, removing all duplicated logging and CLI code
- All 55 Phase 3 BATS tests pass with zero regressions after extraction

## Task Commits

Each task was committed atomically:

1. **Task 1: Create lib/logging.sh** - `df5def3` (feat)
2. **Task 2: Create lib/cli.sh and update hardnormly.sh** - `b5336dd` (feat)

**Plan metadata:** (docs commit — see below)

## Files Created/Modified

- `lib/logging.sh` — Include guard, global state vars, setter functions, log_msg, debug_msg, error_msg, run_cmd (stderr capture + error reporting), run_cmd_with_retry (UCSC MySQL retry)
- `lib/cli.sh` — Include guard, show_help, show_version, parse_args (caller-scope globals pattern), validate_args (uses log_msg from logging.sh), parse_filter_args (nameref for caller array)
- `hardnormly.sh` — Added `_SCRIPT_DIR`, source lines for both modules; removed duplicated log_msg, debug_msg, show_help, arg-parsing loop, validation block; boolean tests changed to string comparison

## Decisions Made

- `parse_args` sets caller-scope globals directly (not namerefs): for 20+ variables called exactly once from orchestrator, this is the simplest correct approach
- `SC2034`/`SC2154` disabled at function level only (single comment above `parse_args`) — not file-wide
- Boolean comparisons changed from `if $var` to `if [[ "$var" == "true" ]]` — required when `debug`, `only_pass` etc. are strings (set by parse_args without `local`, received as strings)
- `run_cmd` uses `${_TMP_DIR:+${_TMP_DIR}/}` expansion so stderr temp files go to `_TMP_DIR` when set, otherwise mktemp default `/tmp` — handles early calls before tmp_dir creation
- ERR trap left in hardnormly.sh during this transition plan — will be removed in 04-04 after all command calls are wrapped in run_cmd

## Deviations from Plan

None — plan executed exactly as written.

## Issues Encountered

- shellcheck SC2034 warnings on `parse_args` variables (set in callers scope, not used in cli.sh). Resolved with targeted `# shellcheck disable=SC2034,SC2154` comment above the function definition.

## User Setup Required

None — no external service configuration required.

## Next Phase Readiness

- `lib/logging.sh` and `lib/cli.sh` are in place — Plans 04-02 through 04-05 can source them
- Include guard pattern established — all subsequent lib/ modules must follow it
- `run_cmd` available for pipeline steps to use in 04-02 (bed.sh) and beyond
- Blocker: `run_cmd_with_retry` not yet used in hardnormly.sh (UCSC MySQL still uses inline `|| {}`) — will be wired in 04-02 or 04-03 when genome.sh is extracted

---
*Phase: 04-refactoring*
*Completed: 2026-02-18*
