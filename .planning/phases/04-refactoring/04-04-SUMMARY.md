---
phase: 04-refactoring
plan: 04
subsystem: infra
tags: [bash, shellcheck, shfmt, makefile, ci, github-actions, cleanup, traps]

# Dependency graph
requires:
  - phase: 04-03
    provides: All 8 lib/ modules extracted; hardnormly.sh slim orchestrator with run_cmd wrapping all external commands
provides:
  - ERR trap removed from hardnormly.sh — cleanup_handler is sole cleanup mechanism
  - Makefile uses $(wildcard lib/*.sh) — auto-discovers new lib/ modules
  - CI shellcheck and shfmt steps use lib/*.sh glob — auto-discovers new lib/ modules
affects: [04-05, release]

# Tech tracking
tech-stack:
  added: []
  patterns: [single-cleanup-handler, wildcard-lib-discovery]

key-files:
  created: []
  modified:
    - hardnormly.sh
    - Makefile
    - .github/workflows/ci.yml

key-decisions:
  - "ERR trap removed; set -Eeuo pipefail kept (the -E flag is harmless without ERR trap)"
  - "cleanup_handler with trap - EXIT guard is the single cleanup path (no inline rm -f anywhere)"
  - "Makefile SH_FILES uses $(wildcard lib/*.sh) — new lib/ files auto-included without manual editing"
  - "CI uses lib/*.sh glob — matches Makefile coverage automatically"

patterns-established:
  - "Single cleanup path: only cleanup_handler on EXIT, no inline cleanup at error sites"
  - "Wildcard discovery: $(wildcard lib/*.sh) in Makefile covers all current and future lib/ modules"

# Metrics
duration: 2min
completed: 2026-02-18
---

# Phase 4 Plan 4: Orchestrator Finalization Summary

**ERR trap removed from hardnormly.sh and Makefile/CI updated to auto-discover lib/*.sh via wildcard/glob**

## Performance

- **Duration:** ~2 min
- **Started:** 2026-02-18T19:47:43Z
- **Completed:** 2026-02-18T19:49:53Z
- **Tasks:** 2
- **Files modified:** 3

## Accomplishments
- Removed `err_handler()` function and `trap 'err_handler ${LINENO}' ERR` from hardnormly.sh — `cleanup_handler` is now the sole cleanup mechanism
- Makefile `SH_FILES` switched from explicit file list to `$(wildcard lib/*.sh)` — no manual updates needed when adding new lib/ modules
- CI `shellcheck` and `shfmt` steps updated to `lib/*.sh` glob — matches Makefile coverage automatically, resolves STATE.md blocker

## Task Commits

Each task was committed atomically:

1. **Task 1: Remove ERR trap, consolidate cleanup handler** - `a331da3` (refactor)
2. **Task 2: Update Makefile and CI to include lib/*.sh** - `7c0da78` (chore)

**Plan metadata:** `[pending]` (docs: complete orchestrator finalization plan)

## Files Created/Modified
- `hardnormly.sh` - Removed err_handler() and ERR trap; cleanup_handler is sole cleanup path; now ~219 lines
- `Makefile` - SH_FILES uses `$(wildcard lib/*.sh)` wildcard instead of explicit list
- `.github/workflows/ci.yml` - shellcheck and shfmt steps use `lib/*.sh` glob

## Decisions Made
- Kept `set -Eeuo pipefail` including `-E` (errtrace) — harmless without ERR trap, provides defense-in-depth if trap is added later
- cleanup_handler already had the `trap - EXIT` guard to prevent recursive invocation; no changes needed there

## Deviations from Plan

None - plan executed exactly as written.

## Issues Encountered
None

## User Setup Required
None - no external service configuration required.

## Next Phase Readiness
- hardnormly.sh is a clean orchestrator: module loading, cleanup handler, arg parsing, 7 pipeline step calls (~219 lines)
- Makefile and CI now auto-discover lib/ modules — the "must add new files manually" blocker is resolved
- All 55 BATS tests pass
- Ready for plan 04-05 (final cleanup / version bump / release prep)

---
*Phase: 04-refactoring*
*Completed: 2026-02-18*
