---
phase: 05-features-docs
plan: 02
subsystem: docs
tags: [bash, wget, bedtools, ENCODE, UCSC, exclusion-bed, genomics, documentation]

# Dependency graph
requires:
  - phase: 04-refactoring
    provides: hardnormly.sh modular architecture with lib/ modules and --exclude-bed option
provides:
  - scripts/generate_exclusion_bed.sh for one-command exclusion BED generation from 4 public sources
  - README Filter File Format section with three-column table and example
  - README Generating Exclusion BED Files section with source table and usage commands
affects:
  - 05-03 (integration tests may reference generate_exclusion_bed.sh)
  - 05-04 (changelog/release notes will include this feature)

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "Download-and-merge pattern: wget | gunzip | awk | bedtools sort > tmp; bedtools merge"
    - "Fallback logic with wc -l threshold for build-specific data availability (hg38 gap table fallback)"
    - "shfmt pipe redirect indentation: redirect > on continuation line gets extra indent level"

key-files:
  created:
    - scripts/generate_exclusion_bed.sh
    - .planning/phases/05-features-docs/05-02-SUMMARY.md
  modified:
    - README.md
    - Makefile
    - .github/workflows/ci.yml

key-decisions:
  - "generate_exclusion_bed.sh uses bedtools sort+merge (not multiinter) for union of all exclusion regions"
  - "hg38 centromere fallback: try gap.txt.gz first, use cytoBandIdeo.txt.gz if <10 rows"
  - "Individual source failures are warnings; only fail if ALL sources return 0 rows"
  - "Added to Makefile SH_FILES explicitly (consistent with other scripts/ entries)"
  - "Added to CI shellcheck/shfmt steps alongside existing explicit file list"
  - "Filter File Format section placed after Filter Expression Example (continues usage topic)"
  - "Old 'Exclusion BED Files Documentation' stub replaced with concrete generate script docs"

patterns-established:
  - "scripts/ additions: always add to both Makefile SH_FILES and CI lint steps"
  - "shfmt pipe-redirect: final redirect in multiline pipe gets one extra indent level vs pipe continuation"

# Metrics
duration: 5min
completed: 2026-02-18
---

# Phase 5 Plan 2: Exclusion BED Helper and README Documentation Summary

**Standalone exclusion BED generator (4 UCSC/ENCODE sources, hg19/hg38) and README sections for filter file format and BED source documentation**

## Performance

- **Duration:** ~5 min
- **Started:** 2026-02-18T20:27:40Z
- **Completed:** 2026-02-18T20:32:49Z
- **Tasks:** 2
- **Files modified:** 4

## Accomplishments
- Created `scripts/generate_exclusion_bed.sh` — downloads ENCODE Blacklist, UCSC SuperDups, RepeatMasker LCR, and UCSC gap/cytoBandIdeo into a merged exclusion BED with hg38 centromere fallback
- Added Filter File Format section to README with three-column table (`filter_name`, `action`, `bcftools_expression`) and example line
- Added Generating Exclusion BED Files section to README with source table, URLs, and usage examples for hg19/hg38
- Replaced outdated "Exclusion BED Files Documentation" stub with the new concrete generate script documentation

## Task Commits

Each task was committed atomically:

1. **Task 1: Create scripts/generate_exclusion_bed.sh helper** - `8436197` (feat)
2. **Task 2: Add filter file format and exclusion BED source documentation to README.md** - `3cab712` (docs)

**Plan metadata:** (this commit) (docs: complete plan)

## Files Created/Modified
- `scripts/generate_exclusion_bed.sh` - Standalone helper to download and merge 4 genomic exclusion sources; supports hg19/hg38 with hg38 fallback for cytoBandIdeo centromeres
- `README.md` - Added Filter File Format subsection; replaced old exclusion BED stub with Generating Exclusion BED Files section
- `Makefile` - Added `scripts/generate_exclusion_bed.sh` to SH_FILES list
- `.github/workflows/ci.yml` - Added `scripts/generate_exclusion_bed.sh` to shellcheck and shfmt steps

## Decisions Made
- `bedtools sort | bedtools merge` chosen over `multiinter` for the merge step — produces union of all exclusion regions as non-overlapping intervals, which is the correct semantics for exclusion regions
- hg38 centromere fallback threshold set at `<10 rows`: hg38 gap table legitimately has few centromere entries (UCSC uses cytoBandIdeo for hg38 reference acrocentric regions)
- Individual source download failures are warnings (not errors) to handle temporary network issues; only fail when ALL sources return empty
- Filter File Format placed as `###` subsection after "Filter Expression Example" (not a new `##` section) to maintain usage topic grouping

## Deviations from Plan

None - plan executed exactly as written.

## Issues Encountered
- Local WSL shellcheck v0.8.0 reports SC2310 on pre-existing `hardnormly.sh` code (`if ! plot_stats_output...`) due to `enable=check-set-e-suppressed` in `.shellcheckrc`. This is a local version artifact — CI Ubuntu runs a newer shellcheck that handles this correctly. The new `scripts/generate_exclusion_bed.sh` passes shellcheck cleanly in isolation.

## User Setup Required

None - no external service configuration required.

## Next Phase Readiness
- `scripts/generate_exclusion_bed.sh` is ready; 05-03 may write integration/smoke tests for it if needed
- README sections complete; 05-04 changelog can reference these as FEAT-05 and DOCS-02/DOCS-03

---
*Phase: 05-features-docs*
*Completed: 2026-02-18*
