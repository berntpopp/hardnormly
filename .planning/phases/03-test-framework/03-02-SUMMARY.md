---
phase: 03-test-framework
plan: 02
subsystem: testing
tags: [bats, vcf, bcftools, gatk, freebayes, filters, unit-tests]

# Dependency graph
requires:
  - phase: 03-01
    provides: BATS infrastructure (test_helper/common.bash, _get_filter helper, smoke tests)
  - phase: 02-04
    provides: Expected output VCFs (tests/data/expected/), synthetic VCFs (tests/data/synthetic/)
provides:
  - GATK filter unit tests (18 tests, TFWK-02) — tests/02-gatk-filters.bats
  - Freebayes filter unit tests (18 tests, TFWK-03) — tests/03-freebayes-filters.bats
  - Per-variant FILTER tag assertions for all 7 GATK filters and all 9 Freebayes filters
affects: [04-refactor, CI, any future filter changes]

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "setup_file() runs pipeline once per bats file; @test blocks query output VCF"
    - "_get_filter <vcf> <chrom> <pos> returns FILTER field string for exact equality assertions"
    - "Inspect expected VCF with bcftools query before writing assertions (no placeholder values)"

key-files:
  created:
    - tests/02-gatk-filters.bats
    - tests/03-freebayes-filters.bats
  modified: []

key-decisions:
  - "Assertions use exact string equality ([  \"$filter\" = \"expected\" ]) not substring matching, catching ordering issues"
  - "Both positive and negative assertions included: positive for each filter trigger, negative on control variant"
  - "Multi-filter combined tag order verified exactly (e.g., NOT_IN_INCLUDE_REGION;VAFu02het;readPosBias)"
  - "setup_file() runs hardnormly.sh pipeline on synthetic VCF; pipeline output stored in BATS_FILE_TMPDIR"

patterns-established:
  - "Filter test pattern: inspect expected VCF first, record exact values, write assertions against those values"
  - "Each filter gets at minimum one positive test (variant IS tagged) and the control gets a negative test (NOT tagged)"
  - "Combined filter tags verified both as a whole (exact equality) and as components ([[ $filter == *tag* ]])"

# Metrics
duration: 3min
completed: 2026-02-18
---

# Phase 3 Plan 2: GATK and Freebayes Filter Unit Tests Summary

**36 BATS tests (18 GATK + 18 Freebayes) verifying per-variant FILTER tags for all 16 filter conditions using exact string equality assertions**

## Performance

- **Duration:** ~3 min
- **Started:** 2026-02-18T17:35:03Z
- **Completed:** 2026-02-18T17:37:45Z
- **Tasks:** 2
- **Files modified:** 2

## Accomplishments
- GATK filter unit tests: 18 tests covering all 7 filters (gatkSNPhard x6 conditions, gatkINDELhard x3 conditions, DPu10het, DPu5hom, VAFu02het, VAFo08het, VAFu095hom, region filters, negatives)
- Freebayes filter unit tests: 18 tests covering all 9 filters (lowQUAL, QUALperAO, strandBias x2, readPosBias, DPu10het, DPu5hom, VAFu02het, VAFo08het, VAFu095hom, region filters, negatives)
- Multi-filter combined tag verification: exact order checked for `IN_EXCLUDE_REGION;gatkINDELhard`, `NOT_IN_INCLUDE_REGION;VAFu02het;readPosBias`

## Task Commits

Each task was committed atomically:

1. **Task 1: Write 02-gatk-filters.bats** - `7c36c5d` (test)
2. **Task 2: Write 03-freebayes-filters.bats** - `a4df6d9` (test)

**Plan metadata:** (docs commit follows this summary)

## Files Created/Modified
- `tests/02-gatk-filters.bats` - 18 GATK filter unit tests using setup_file() pattern
- `tests/03-freebayes-filters.bats` - 18 Freebayes filter unit tests using setup_file() pattern

## Decisions Made
- Used exact string equality assertions (`[ "$filter" = "expected" ]`) rather than substring matching — catches filter ordering regressions
- Both positive and negative assertions: each filter has a positive test verifying the tag IS applied, plus negative tests on the control variant (PASS at POS 10101) confirming no false positives
- Multi-filter tags verified twice: once as exact combined string, once as component substring checks (for readability of what failed)
- Pipeline invoked via `"$HARDNORMLY"` in setup_file() writing to `$BATS_FILE_TMPDIR` — output persists for all @test blocks in file

## Deviations from Plan

None - plan executed exactly as written.

The critical first step (inspect actual expected output before writing assertions) was followed precisely: `bcftools query` on the expected VCFs produced exact FILTER strings that were directly embedded in test assertions.

## Issues Encountered
None.

## User Setup Required
None - no external service configuration required.

## Next Phase Readiness
- Filter unit tests (TFWK-02 and TFWK-03) complete and all passing
- 03-03 (remaining BATS tests: regression, real data, edge cases) can proceed
- Phase 4 (refactoring) has safety net: any filter behavior change will cause these tests to fail

---
*Phase: 03-test-framework*
*Completed: 2026-02-18*
