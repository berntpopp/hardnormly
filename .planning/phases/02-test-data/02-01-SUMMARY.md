---
phase: 02-test-data
plan: 01
subsystem: testing
tags: [samtools, bedtools, fasta, vcf, bed, gitignore, test-data]

# Dependency graph
requires: []
provides:
  - .gitignore exceptions allowing tests/ directory and test data file extensions to be committed
  - tests/data/synthetic/mini_ref.fa - chr22 100kb FASTA extracted from hs37d5.fa, contig named 22
  - tests/data/synthetic/mini_ref.fa.fai - FASTA index for samtools/bcftools norm
  - tests/data/include_regions.bed - 3 include regions (10000-90000) for region annotation testing
  - tests/data/exclude_regions.bed - 2 exclude regions (10000-90000) for region annotation testing
  - tests/data/hg19_chr22.genome - chr22 length for bedtools slop operations
  - tests/data/{synthetic,real,expected}/ directory structure
affects:
  - 02-02 (synthetic VCFs need mini_ref.fa for bcftools norm)
  - 02-03 (expected outputs need BED files and genome file)
  - 02-04 (BATS tests need all files in tests/data/)
  - all subsequent Phase 2 plans

# Tech tracking
tech-stack:
  added: [samtools faidx (extraction + indexing)]
  patterns: [gitignore negation ordering - exceptions must come AFTER the pattern they override]

key-files:
  created:
    - tests/data/synthetic/mini_ref.fa
    - tests/data/synthetic/mini_ref.fa.fai
    - tests/data/include_regions.bed
    - tests/data/exclude_regions.bed
    - tests/data/hg19_chr22.genome
    - tests/data/synthetic/.gitkeep
    - tests/data/real/.gitkeep
    - tests/data/expected/.gitkeep
  modified:
    - .gitignore

key-decisions:
  - "Gitignore negations must appear AFTER *test* rule (not before) - later rules win in git"
  - "Extract chr22:1-100001 from hs37d5.fa (not mid-chromosome) so contig positions 1-100001 align with VCF positions"
  - "Used samtools via WSL since bioinformatics tools not available natively on Windows"

patterns-established:
  - "Test data in tests/data/{synthetic,real,expected}/ directory structure"
  - "Gitignore pattern: *test* rule first, then !/tests/ exceptions after"

# Metrics
duration: 25min
completed: 2026-02-18
---

# Phase 2 Plan 01: Test Data Foundation Summary

**Tab-delimited chr22 test infrastructure: mini_ref.fa (100001 bases from hs37d5.fa), 3-region include BED, 2-region exclude BED, and hg19 genome file - with gitignore fixed so tests/ can be committed**

## Performance

- **Duration:** ~25 min
- **Started:** 2026-02-18T (session start)
- **Completed:** 2026-02-18
- **Tasks:** 2
- **Files modified:** 9 (1 modified, 8 created)

## Accomplishments

- Fixed gitignore negation ordering so `!/tests/` exceptions override `*test*` rule
- Created `tests/data/{synthetic,real,expected}/` directory structure tracked by git
- Extracted 100kb chr22 FASTA from hs37d5.fa with proper `>22` header and samtools index
- Created tab-delimited BED and genome files for region annotation and bedtools slop testing

## Task Commits

Each task was committed atomically:

1. **Task 1: Fix gitignore and create directory structure** - `0f28cc4` (chore)
2. **Task 2: Create reference FASTA, BED files, and genome file** - `8e8ad4a` (feat)

**Plan metadata:** (see final commit below)

## Files Created/Modified

- `.gitignore` - Added `!/tests/` and extension exceptions after `*test*` rule
- `tests/data/synthetic/mini_ref.fa` - chr22 100001 bases from hs37d5.fa, header `>22`
- `tests/data/synthetic/mini_ref.fa.fai` - FASTA index (22, 100001, 4, 60, 61)
- `tests/data/include_regions.bed` - 3 regions: 22:10000-30000, 22:40000-60000, 22:70000-90000
- `tests/data/exclude_regions.bed` - 2 regions: 22:20000-25000, 22:80000-85000
- `tests/data/hg19_chr22.genome` - Single line: 22\t51304566
- `tests/data/synthetic/.gitkeep` - Directory placeholder
- `tests/data/real/.gitkeep` - Directory placeholder
- `tests/data/expected/.gitkeep` - Directory placeholder

## Decisions Made

- **Gitignore negation ordering:** Git processes `.gitignore` top-to-bottom; the LAST matching rule wins. The plan specified putting `!/tests/` BEFORE `*test*`, but that would be overridden. The correct ordering places `!/tests/` AFTER `*test*` so the negation wins. This is a fundamental git gitignore rule.

- **chr22:1-100001 extraction:** Used the first 100kb of chr22 (positions 1-100001) rather than a mid-chromosome region. After renaming the header to `>22`, bcftools sees contig 22 with positions 1-100001 - matching where synthetic VCFs will place variants (10001-90000).

- **WSL for bioinformatics tools:** samtools not available natively on Windows; used `wsl bash -c "samtools ..."` to run samtools from WSL while writing files to the Windows filesystem via `/mnt/c/...` paths.

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 1 - Bug] Fixed gitignore negation ordering (exceptions must come AFTER the pattern they override)**

- **Found during:** Task 1 (Fix gitignore)
- **Issue:** Plan specified placing `!/tests/` immediately BEFORE `*test*`, but git's last-match-wins rule means `*test*` (which comes after) would re-ignore files under tests/. Verified with `git check-ignore -v` which showed `*test*` still matching `tests/data/include_regions.bed` even with `!/tests/**` present before it.
- **Fix:** Moved `!/tests/` and all file-extension exceptions to AFTER the `*test*` rule (lines 21-29 instead of lines 18-25). This ensures the negations are the last matching rule for files under tests/.
- **Files modified:** .gitignore
- **Verification:** `git status` shows `tests/` as `??` (untracked, not ignored). `git add tests/` successfully stages all test files.
- **Committed in:** `0f28cc4` (Task 1 commit)

---

**Total deviations:** 1 auto-fixed (Rule 1 - Bug: gitignore ordering)
**Impact on plan:** Required fix - without it, no test files could be committed to the repo.

## Issues Encountered

- samtools not available natively on Windows; used WSL (`wsl bash -c "samtools ..."`) with `/mnt/c/development/hardnormly/` paths. chr22:1-100001 is mostly Ns (telomeric region) but provides real sequence for positions beyond ~2MB where actual gene sequence exists; for synthetic VCF testing the sequence content doesn't matter as long as the FASTA is valid and indexed.

## User Setup Required

None - no external service configuration required.

## Next Phase Readiness

- All shared foundation files are committed and gitignore-unblocked
- Plan 02-02 can proceed to create synthetic VCFs (needs mini_ref.fa for bcftools norm)
- Plan 02-03 can proceed to generate expected outputs using BED and genome files
- Plan 02-04 can proceed to write BATS tests using the full tests/data/ structure
- No blockers

---
*Phase: 02-test-data*
*Completed: 2026-02-18*
