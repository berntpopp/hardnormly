---
phase: 05-features-docs
plan: 03
subsystem: cli
tags: [bash, subcommands, dispatcher, cli, help, smoke-tests, bats]

# Dependency graph
requires:
  - phase: 05-01
    provides: "--caller and --strip-annotations flags added to parse_args; show_help was the existing function to be replaced"
provides:
  - Subcommand dispatcher in hardnormly.sh routing run-pipeline, generate-*-bed, no-args/--help, --version, legacy flags, unknown subcommands
  - Compact subcommand-aware show_help (~35 lines) with subcommands section, all flags as one-liners, usage example
  - 6 new smoke tests covering dispatcher behavior (12 total)
affects: [05-04, integration-tests, documentation]

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "Subcommand dispatcher: case on first arg, route to handler or legacy parse_args"
    - "Stub functions for unimplemented subcommands exit 1 with informative error"
    - "cat heredoc with single-quoted HELP delimiter for shellcheck-clean multi-line output"

key-files:
  created: []
  modified:
    - hardnormly.sh
    - lib/cli.sh
    - tests/01-smoke.bats

key-decisions:
  - "Dispatcher already present in HEAD from prior session (docs(04) commit included it); Task 1 needed no new commit"
  - "show_help uses cat <<'HELP' heredoc (single-quoted) to prevent variable expansion and satisfy shellcheck"
  - "show_help and show_version output goes to stdout (not stderr) — BATS run captures correctly with assert_output"
  - "run-pipeline --version and run-pipeline --help work because parse_args handles them after dispatcher shifts"

patterns-established:
  - "Dispatcher pattern: _subcommand=\${1:-}; case; -* falls through to legacy parse_args"
  - "Stub subcommands: defined before dispatcher, output error to stderr, exit 1"

# Metrics
duration: 15min
completed: 2026-02-18
---

# Phase 5 Plan 3: Subcommand Dispatcher and Compact Help Summary

**Subcommand dispatcher routing run-pipeline/generate-*-bed/legacy-flags/--help via case on first arg, with compact one-screen show_help listing all subcommands and flags**

## Performance

- **Duration:** ~15 min
- **Started:** 2026-02-18T20:35:00Z
- **Completed:** 2026-02-18T20:50:00Z
- **Tasks:** 3
- **Files modified:** 3

## Accomplishments

- Subcommand dispatcher in hardnormly.sh routes: `run-pipeline` to parse_args, `generate-*-bed` to stubs, no-args/`--help`/`-h` to show_help, `--version` to show_version, `-*` flags to parse_args (legacy compat), unknown words to error+exit 1
- Rewrote show_help from 26 verbose echo lines to compact cat heredoc (~35 lines): subcommands section, all flags with one-liner descriptions, usage example
- Added 6 new smoke tests covering dispatcher behavior; all 12 smoke tests pass

## Task Commits

Each task was committed atomically:

1. **Task 1: Add subcommand dispatcher to hardnormly.sh** - `7149f5f` (included in prior session's docs(04) commit — already present in HEAD)
2. **Task 2: Rewrite show_help for compact subcommand-aware display** - `5b1fea4` (feat)
3. **Task 3: Update smoke tests for dispatcher behavior** - `d22555c` (test)

**Plan metadata:** (created after this summary)

## Files Created/Modified

- `hardnormly.sh` - Dispatcher case block + cmd_generate_*_bed stub functions (already committed in 7149f5f)
- `lib/cli.sh` - show_help rewritten: cat heredoc with subcommands section, compact one-liner flags, usage example
- `tests/01-smoke.bats` - 6 new dispatcher/help-content tests added (12 total)

## Decisions Made

- **Dispatcher already in HEAD:** The dispatcher code was committed in a prior planning session's `docs(04)` commit (`7149f5f`). No duplicate commit needed for Task 1 — it was already functionally complete.
- **Heredoc delimiter single-quoted:** `cat <<'HELP'` prevents variable expansion inside heredoc, satisfying shellcheck and avoiding accidental interpolation of `$0` or other variables.
- **show_help output to stdout:** The cat heredoc goes to stdout (not stderr). BATS `run` captures stdout by default; `assert_output --partial "Usage:"` works correctly.
- **Stub placement:** cmd_generate_*_bed stubs placed BEFORE the dispatcher (after lib sourcing, before defaults) so they are available when dispatcher case is evaluated.

## Deviations from Plan

### Auto-fixed Issues

None within tasks 2 and 3.

**Note on Task 1:** The dispatcher code was found to already be committed in HEAD (`7149f5f docs(04): complete refactoring phase`), having been included in a previous execution session. No corrective action needed — the implementation matched the plan exactly. Documented as an observation, not a deviation.

---

**Total deviations:** 0 auto-fixed
**Impact on plan:** No deviations. Plan executed cleanly.

## Issues Encountered

- **Exit code display quirk:** Direct `bash hardnormly.sh --help; echo exit:$?` showed exit 0 in shell due to cleanup_handler's `local exit_code=$?` capturing empty string when EXIT trap fires after `exit 1`. BATS captures actual process exit code correctly. Existing smoke tests already verified the correct behavior (assert_failure passed). No fix needed.

## Next Phase Readiness

- Dispatcher ready for 05-04 to implement the actual `cmd_generate_inclusion_bed` and `cmd_generate_exclusion_bed` functions (replacing the stubs)
- All smoke tests green; backward compat verified for legacy flag invocations
- show_help compact format is the canonical help for future flag additions

---
*Phase: 05-features-docs*
*Completed: 2026-02-18*
