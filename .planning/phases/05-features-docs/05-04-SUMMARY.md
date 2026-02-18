---
phase: 05-features-docs
plan: 04
subsystem: cli
tags: [bash, bedtools, subcommands, bed-merging, argument-parsing]

# Dependency graph
requires:
  - phase: 05-03
    provides: subcommand dispatcher and stubs for generate-inclusion-bed and generate-exclusion-bed
  - phase: 04-refactoring
    provides: lib/bed.sh with normalize_bed, merge_include_beds, merge_exclude_beds
  - phase: 04-refactoring
    provides: lib/cli.sh with show_help, parse_args, validate_args
provides:
  - cmd_generate_inclusion_bed: fully implemented in hardnormly.sh
  - cmd_generate_exclusion_bed: fully implemented in hardnormly.sh
  - show_help_generate_inclusion_bed: per-subcommand help in lib/cli.sh (exits 0)
  - show_help_generate_exclusion_bed: per-subcommand help in lib/cli.sh (exits 0)
affects: []

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "Per-subcommand argument parsers: each subcommand has own while/case loop (not reusing parse_args)"
    - "Subcommand help exits 0 (informational), main show_help exits 1 (error path)"
    - "SC2064 disable inline for trap with captured variable: trap \"rm -rf '$tmp_dir'\" EXIT"
    - "Silent by default pattern: log_msg calls guarded by [[ \"$verbose\" == true ]]"

key-files:
  created: []
  modified:
    - lib/cli.sh
    - hardnormly.sh

key-decisions:
  - "Per-subcommand help exits 0 (informational), unlike main show_help which exits 1"
  - "Subcommand arg parsers are self-contained — not reusing parse_args (different flag sets)"
  - "genome_file required for generate-inclusion-bed (slop needs it) but not for generate-exclusion-bed"
  - "No bgzip/tabix step in generate-*-bed subcommands — output is plain BED, users compress separately"
  - "Local tmp_dir per subcommand with inline trap cleanup — EXIT trap fires on subcommand exit"

patterns-established:
  - "Per-subcommand help: show_help_<subcommand> exits 0 (not 1)"
  - "Subcommand handlers: local vars, own arg parser, validate required args, call lib/ functions"

# Metrics
duration: 4min
completed: 2026-02-18
---

# Phase 5 Plan 4: Generate BED Subcommands Summary

**generate-inclusion-bed and generate-exclusion-bed subcommands expose BED merging as standalone utilities via per-subcommand arg parsers that wrap lib/bed.sh normalize_bed + merge_include/exclude_beds**

## Performance

- **Duration:** ~4 min
- **Started:** 2026-02-18T20:45:39Z
- **Completed:** 2026-02-18T20:49:11Z
- **Tasks:** 2
- **Files modified:** 2

## Accomplishments
- Added `show_help_generate_inclusion_bed` and `show_help_generate_exclusion_bed` to lib/cli.sh (both exit 0)
- Replaced stubs in hardnormly.sh with full `cmd_generate_inclusion_bed` and `cmd_generate_exclusion_bed` implementations
- Each subcommand has its own argument parser, required-arg validation, and calls lib/bed.sh merging functions
- Silent by default; `-v/--verbose` enables progress messages
- All 12 existing smoke tests pass without regression

## Task Commits

Each task was committed atomically:

1. **Task 1: Add per-subcommand help functions to lib/cli.sh** - `705b881` (feat)
2. **Task 2: Implement cmd_generate_inclusion_bed and cmd_generate_exclusion_bed** - `a0437f7` (feat)

## Files Created/Modified
- `lib/cli.sh` - Added `show_help_generate_inclusion_bed` and `show_help_generate_exclusion_bed`; updated Provides: comment
- `hardnormly.sh` - Replaced stubs with full implementations of both subcommand handlers

## Decisions Made
- Per-subcommand help exits 0 (informational) unlike main show_help which exits 1 (error path)
- Each subcommand has its own while/case argument parser rather than reusing parse_args — the flag sets differ (-b/-g vs -e, neither uses VCF/FASTA)
- genome_file is REQUIRED for generate-inclusion-bed (slop operation needs it) but NOT required for generate-exclusion-bed (no slop)
- No bgzip/tabix step in subcommands — output is plain BED; users can compress separately if needed
- SC2064 disabled inline for `trap "rm -rf '$tmp_dir'" EXIT` — double quotes intentionally capture current value of $tmp_dir at trap definition time

## Deviations from Plan

None - plan executed exactly as written.

## Issues Encountered

Minor: Initial exit-code verification appeared to show exit 0 due to WSL bash -c test artifact (the `echo exit:\$?` in the same -c string ran as last command, returning its own exit 0). Confirmed actual exit code is 1 using WSL heredoc invocation. No code change needed.

## User Setup Required

None - no external service configuration required.

## Next Phase Readiness
- Phase 5 complete: all 4 plans finished (05-01 --caller/--strip-annotations, 05-02 generate_exclusion_bed.sh, 05-03 dispatcher, 05-04 subcommand handlers)
- Full refactor/v0.7.0 feature set implemented
- No blockers for final PR/merge

---
*Phase: 05-features-docs*
*Completed: 2026-02-18*
