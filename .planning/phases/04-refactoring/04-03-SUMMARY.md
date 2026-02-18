---
phase: 04-refactoring
plan: 03
subsystem: refactoring
tags: [bash, bcftools, lib-modules, vcf-normalization, filter-pipeline, stats]

# Dependency graph
requires:
  - phase: 04-02
    provides: lib/genome.sh, lib/bed.sh, lib/annotate.sh (Steps 1-4 extracted)
provides:
  - lib/normalize.sh with normalize_vcf (special stderr handling, no run_cmd)
  - lib/filter.sh with init_filter_pipeline, apply_filter_stages, write_filtered_output
  - lib/stats.sh with generate_stats, plot_stats_output
  - All 8 lib/ modules exist; hardnormly.sh is a slim orchestrator
affects: [04-04, 04-05]

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "normalize_vcf avoids run_cmd — bcftools norm emits warnings on stderr even on success; captured manually"
    - "apply_filter_stages receives filter_stages as positional args (not nameref) for clean function signature"
    - "plot_stats_output creates its own temp file in tmp_dir (not caller-managed) for clean resource ownership"

key-files:
  created:
    - lib/normalize.sh
    - lib/filter.sh
    - lib/stats.sh
  modified:
    - hardnormly.sh
    - Makefile

key-decisions:
  - "normalize_vcf does NOT use run_cmd — bcftools norm emits warnings on stderr even on success; manual capture needed"
  - "apply_filter_stages takes stages as positional args after tmp_dir — clean varargs signature"
  - "Makefile SH_FILES updated to include all 8 lib/ modules"

patterns-established:
  - "normalize_vcf pattern: mktemp in tmp_dir for stderr/stdout, grep for Warning/Lines after success"
  - "filter.sh pipeline: init (fill-tags) → apply stages → write output as 3 separate functions"
  - "stats.sh: generate_stats uses run_cmd (redirecting stdout), plot_stats_output creates own temp file"

# Metrics
duration: 5min
completed: 2026-02-18
---

# Phase 4 Plan 03: Normalize/Filter/Stats Extraction Summary

**VCF normalization, filter pipeline, and stats generation extracted into lib/normalize.sh, lib/filter.sh, and lib/stats.sh — hardnormly.sh reduced to slim orchestrator sourcing 8 lib modules**

## Performance

- **Duration:** 5 min
- **Started:** 2026-02-18T19:39:45Z
- **Completed:** 2026-02-18T19:44:35Z
- **Tasks:** 2/2 completed
- **Files modified:** 5 (3 created, 2 modified)

## Accomplishments

- Created lib/normalize.sh with normalize_vcf — special stderr/stdout capture bypasses run_cmd since bcftools norm emits warnings on stderr even on success
- Created lib/filter.sh with init_filter_pipeline (fill-tags init), apply_filter_stages (sequential BCF filter loop), write_filtered_output (format detection, PASS filter, auto-index)
- Created lib/stats.sh with generate_stats (run_cmd bcftools stats) and plot_stats_output (plot-vcfstats with output capture)
- Updated hardnormly.sh to source all 8 lib/ modules and delegate Steps 5-7 to module functions; 105 lines of inline pipeline code removed
- All 55 Phase 3 tests pass with tools available

## Task Commits

1. **Task 1: Create lib/normalize.sh, lib/filter.sh, lib/stats.sh** - `7784940` (feat)
2. **Task 2: Update hardnormly.sh to use normalize/filter/stats modules** - `b7ffb07` (feat)

**Plan metadata:** (pending final docs commit)

## Files Created/Modified

- `lib/normalize.sh` — normalize_vcf with manual stderr/stdout capture (bcftools norm special case)
- `lib/filter.sh` — init_filter_pipeline, apply_filter_stages, write_filtered_output
- `lib/stats.sh` — generate_stats, plot_stats_output
- `hardnormly.sh` — sources all 8 lib/ modules; Steps 5-7 replaced with module function calls
- `Makefile` — SH_FILES updated to include lib/normalize.sh, lib/filter.sh, lib/stats.sh

## Decisions Made

- normalize_vcf does NOT use run_cmd because bcftools norm emits warnings on stderr even on success — run_cmd would silently discard those warnings; manual mktemp capture preserves Warning/Lines logging
- apply_filter_stages takes stages as positional args (not nameref) — function already accepts "$@" after tmp_dir which is clean and testable
- Makefile SH_FILES updated to include all 8 lib/ modules for consistent lint coverage

## Deviations from Plan

None — plan executed exactly as written.

## Next Phase Readiness

- All 8 lib/ modules complete: logging.sh, cli.sh, genome.sh, bed.sh, annotate.sh, normalize.sh, filter.sh, stats.sh
- hardnormly.sh is now a slim orchestrator (227 lines with comments vs pre-refactoring ~320 lines)
- Ready for 04-04 (ERR trap removal / run_cmd integration cleanup)
- TFWK-05 regression tests pass — filter pipeline behavioral identity confirmed

---
*Phase: 04-refactoring*
*Completed: 2026-02-18*
