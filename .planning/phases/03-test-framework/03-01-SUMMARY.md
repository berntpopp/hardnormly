---
phase: 03-test-framework
plan: 01
subsystem: testing
tags: [bats, bats-core, bats-assert, bats-support, bats-file, git-submodules, ci, shellcheck]

# Dependency graph
requires:
  - phase: 01-infrastructure
    provides: hardnormly.sh with --help, --version, and argument parsing
  - phase: 02-test-data
    provides: test data in tests/data/ (used by later test plans, not this one)
provides:
  - BATS test infrastructure (4 git submodules: bats-core, bats-assert, bats-support, bats-file)
  - Shared test helper (tests/test_helper/common.bash) with _common_setup, _require_tools, _get_filter
  - Smoke test suite (tests/01-smoke.bats) with 6 CLI validation tests
  - CI test job running BATS tests with bcftools/bedtools installed
  - Makefile test and test-debug targets
affects:
  - 03-02-PLAN (GATK filter tests — loads common.bash, runs via same bats binary)
  - 03-03-PLAN (Freebayes filter tests — same pattern)
  - 03-04-PLAN (Integration tests — same infrastructure)

# Tech tracking
tech-stack:
  added:
    - bats-core v1.13.0 (git submodule at tests/bats)
    - bats-assert v2.2.4 (git submodule at tests/test_helper/bats-assert)
    - bats-support v0.3.0 (git submodule at tests/test_helper/bats-support)
    - bats-file v0.2.0 (git submodule at tests/test_helper/bats-file)
  patterns:
    - Shared BATS helper pattern: common.bash sourced by every test file in setup()
    - Tool-availability skip pattern: _require_tools() skips if bcftools/bedtools missing
    - Smoke tests run without bioinformatics tools (no _require_tools call in 01-smoke.bats)
    - REPO_ROOT computed from BATS_TEST_FILENAME (portable, works from any directory)

key-files:
  created:
    - tests/test_helper/common.bash
    - tests/01-smoke.bats
    - .gitmodules
  modified:
    - .github/workflows/ci.yml
    - Makefile
    - .gitattributes

key-decisions:
  - "Git submodules for BATS libraries (not bats-action) — self-contained, uses load not bats_load_library"
  - "_get_filter takes 3 args (vcf, chrom, pos) for forward-compatibility with hg38 chr-prefixed VCFs"
  - ".gitattributes eol=lf for *.sh *.bash *.bats — fixes CRLF corruption on Windows dev with core.autocrlf=true"
  - "CI test job uses apt-get bcftools/bedtools (not conda) — sufficient for BATS tests, simpler CI"

patterns-established:
  - "Pattern: all test files load 'test_helper/common' and call _common_setup in setup()"
  - "Pattern: smoke tests omit _require_tools; integration tests call _require_tools to skip gracefully"
  - "Pattern: HARDNORMLY, SYNTH, REAL_DATA, EXPECTED, TEST_TMP exported by _common_setup for use in @test blocks"

# Metrics
duration: 25min
completed: 2026-02-18
---

# Phase 3 Plan 01: BATS Infrastructure and Smoke Tests Summary

**BATS testing framework installed via 4 git submodules; 6-test smoke suite (TFWK-01) passes in WSL and will pass in CI**

## Performance

- **Duration:** ~25 min
- **Started:** 2026-02-18T00:00:00Z
- **Completed:** 2026-02-18T00:25:00Z
- **Tasks:** 2
- **Files modified:** 7

## Accomplishments

- Installed bats-core, bats-assert, bats-support, bats-file as git submodules (tests/ directory)
- Created `tests/test_helper/common.bash` providing _common_setup, _require_tools, _get_filter helpers
- Wrote `tests/01-smoke.bats` with 6 passing tests for CLI argument validation (--help, --version, missing args, unknown flags)
- Updated CI test job: checkout with submodules:recursive, install bcftools/bedtools, run bats
- Added `make test` and `make test-debug` targets to Makefile

## Task Commits

Each task was committed atomically:

1. **Task 1: Install BATS submodules, create common.bash helper, write 01-smoke.bats** - `a96a1f0` (feat)
2. **Task 2: Update CI workflow and Makefile with BATS test targets** - `92e1a99` (feat)

**Plan metadata:** (docs commit below)

## Files Created/Modified

- `tests/bats/` - bats-core submodule (v1.13.0) — the BATS test runner
- `tests/test_helper/bats-support/` - bats-support submodule (v0.3.0) — diagnostic output for bats-assert
- `tests/test_helper/bats-assert/` - bats-assert submodule (v2.2.4) — assert_success, assert_failure, assert_output
- `tests/test_helper/bats-file/` - bats-file submodule (v0.2.0) — file existence assertions
- `tests/test_helper/common.bash` - Shared BATS helper with _common_setup, _require_tools, _get_filter
- `tests/01-smoke.bats` - 6 CLI smoke tests (TFWK-01)
- `.gitmodules` - Submodule registrations for 4 BATS libraries
- `.github/workflows/ci.yml` - Test job replaced: submodules checkout, bcftools/bedtools install, bats run
- `Makefile` - Added test and test-debug targets, updated .PHONY
- `.gitattributes` - Added eol=lf rules for *.sh, *.bash, *.bats

## Decisions Made

- **Git submodules for BATS** (not bats-action): Self-contained approach using `load` not `bats_load_library`. Requires `submodules: recursive` in CI checkout. Simpler than managing BATS_LIB_PATH.
- **_get_filter takes chrom + pos args**: Forward-compatible with hg38 VCFs using "chr22" prefix. Current synthetic VCFs use chrom="22".
- **CI uses apt-get bcftools/bedtools**: No conda setup needed in CI; apt packages sufficient for BATS filter/integration tests.

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 3 - Blocking] Added eol=lf to .gitattributes for shell/bats files**

- **Found during:** Task 1 (running smoke tests)
- **Issue:** core.autocrlf=true on Windows converted LF to CRLF in submodule bash files, causing `bash\r: No such file or directory` errors in WSL
- **Fix:** Added `*.sh`, `*.bash`, `*.bats` rules with `eol=lf` to `.gitattributes` to prevent CRLF conversion. Converted existing CRLF files in submodule working copies via dos2unix for local testing.
- **Files modified:** `.gitattributes`
- **Verification:** Smoke tests pass via WSL after fix; CI (Linux) will not have this issue since autocrlf defaults to false
- **Committed in:** a96a1f0 (Task 1 commit)

---

**Total deviations:** 1 auto-fixed (1 blocking)
**Impact on plan:** Essential fix for Windows development. CI on Linux is unaffected.

## Issues Encountered

- BATS submodule files had CRLF line endings on Windows due to `core.autocrlf=true`. Fixed with `.gitattributes` eol=lf rules and dos2unix conversion for the working copy. The in-repo submodule commit pointers are unaffected; on Linux CI the files will have correct LF endings.

## User Setup Required

None - no external service configuration required.

## Next Phase Readiness

- BATS infrastructure is installed and working; `tests/bats/bin/bats tests/01-smoke.bats` passes all 6 tests
- `tests/test_helper/common.bash` is ready for 03-02 and 03-03 to import via `load 'test_helper/common'`
- CI test job is ready and will run all `tests/*.bats` on push/PR
- 03-02 (GATK filter tests) can proceed immediately — no blockers

---
*Phase: 03-test-framework*
*Completed: 2026-02-18*
