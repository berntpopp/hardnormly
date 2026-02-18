---
phase: 03-test-framework
plan: 03
subsystem: testing
tags: [bats, integration-testing, regression-testing, 1000G, real-data, vcf, bcftools]

# Dependency graph
requires:
  - phase: 03-01
    provides: BATS infrastructure, _common_setup, _require_tools, _get_filter helpers
  - phase: 02-04
    provides: Expected golden files in tests/data/expected/ and real data in tests/data/real/

provides:
  - 13 BATS integration/regression/edge-case tests in tests/04-integration.bats
  - TFWK-04 coverage: full pipeline on real 1000G NA12878 data with hs37d5.fa
  - TFWK-05 coverage: regression diff against tests/data/expected/gatk_sample_A_filtered.vcf.gz
  - TFWK-06 coverage: -g/--genome flag verification closes issue #13
  - TFWK-07 coverage: empty VCF, no filters, no BED edge cases
affects:
  - 04-ci (CI configuration must know tests require ref/hs37d5.fa for TFWK-04)
  - 05-release (all 55 tests must pass before release)

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "Integration tests use setup() not setup_file() — each test is an independent pipeline run"
    - "Real-data tests skip via HAVE_FULL_REF guard when ref/hs37d5.fa absent (e.g., CI without ref/)"
    - "Regression failure reports diff count and first 20 changed lines for diagnosis"
    - "TFWK guard comments in test names enable traceability between tests and requirements"

key-files:
  created:
    - tests/04-integration.bats
  modified: []

key-decisions:
  - "Used ref/hs37d5.fa (not tests/data/real/chr22_16M_ref.fa) for real data — mini_ref only covers chr22:1-100001 but 1000G variants are at 16M+"
  - "1000G NA12878 run without filters (VCF has only GT FORMAT field; GATK filters require DP/VAF)"
  - "Regression test uses synthetic gatk_sample_A vs expected/ golden file (deterministic, tool-independent)"
  - "HAVE_FULL_REF guard allows TFWK-04 tests to skip in CI environments without ref/ directory"
  - "Regression failure message includes diff count + first 20 lines for immediate diagnosis"

patterns-established:
  - "Setup-per-test (setup()) pattern for tests with independent pipeline invocations"
  - "Resource availability guard: export HAVE_FULL_REF and skip in each test that needs it"
  - "Regression comparison via bcftools query CHROM/POS/FILTER (not binary diff — avoids header timestamp noise)"

# Metrics
duration: 25min
completed: 2026-02-18
---

# Phase 3 Plan 3: Integration, Regression, Genome Flag, and Edge Case Tests Summary

**13 BATS tests covering full pipeline on 1000G real data, regression diff against golden file, -g/--genome flag verification (issue #13), and 5 edge cases — all 55 tests in the suite pass**

## Performance

- **Duration:** ~25 min
- **Started:** 2026-02-18T18:30:00Z
- **Completed:** 2026-02-18T18:55:00Z
- **Tasks:** 2 (both complete)
- **Files modified:** 1 created

## Accomplishments

- 13 integration/regression/edge-case tests in `tests/04-integration.bats`
- Full 55-test suite passes (01-smoke: 6, 02-gatk-filters: 18, 03-freebayes-filters: 18, 04-integration: 13)
- All 4 TFWK requirements covered: TFWK-04 (integration), TFWK-05 (regression), TFWK-06 (genome flag), TFWK-07 (edge cases)
- Issue #13 verified resolved: both `-g` and `--genome` flags work and produce identical output

## Task Commits

Each task was committed atomically:

1. **Task 1: Write 04-integration.bats** - `9d4a773` (feat)
2. **Task 2: Run full test suite** - (verified during Task 1, no separate commit needed)

**Plan metadata:** (docs commit follows)

## Files Created/Modified

- `tests/04-integration.bats` - Integration (TFWK-04), regression (TFWK-05), genome flag (TFWK-06), and edge case (TFWK-07) tests (13 tests, 295 lines)

## Decisions Made

- **ref/hs37d5.fa for real data tests:** `tests/data/real/chr22_16M_ref.fa` only covers chr22:1-100001 but 1000G variants start at 16M+. The full `ref/hs37d5.fa` is required. Tests skip if file is absent (CI without large reference files).
- **No filters for 1000G data:** The 1000G VCF has only GT FORMAT field; GATK filters reference DP/VAF/etc and would fail. Integration tests verify pipeline mechanics (normalization, annotation, BED handling) not filter accuracy.
- **Synthetic data for regression:** Regression test uses `gatk_sample_A` (synthetic) because it has known, deterministic filter assignments. Real data filtering would depend on FORMAT fields absent from 1000G VCFs.
- **HAVE_FULL_REF guard:** Allows TFWK-04 tests to skip gracefully in CI environments without `ref/` directory while still passing.
- **Regression comparison strategy:** Extract CHROM/POS/FILTER with bcftools query instead of binary diff — avoids false failures from header timestamp differences between runs.

## Deviations from Plan

None - plan executed exactly as written. The BED coordinate mismatch noted in the plan (include_regions.bed targets synthetic chr22:10K-90K range, real data is at 16M+) was handled as the plan suggested: run without BED files for 1000G tests, use GIAB high-confidence BED for GIAB-specific test.

## Issues Encountered

- **`chr22_16M_ref.fa` only covers chr22:1-100001:** Discovered when testing pipeline on real 1000G data — bcftools norm reported "Reference allele mismatch at 22:16053659 .. REF_SEQ:'' vs VCF:'A'". Resolved by using `ref/hs37d5.fa` which covers full chr22. Added HAVE_FULL_REF guard so tests skip in CI without the large reference.
- **BATS tool PATH:** Tests skip when run without bioinformatics tools in PATH. Resolved by running BATS with `PATH=/home/bernt/miniconda3/envs/hardnormly/bin:$PATH` prefix — `_require_tools` correctly skips when tools are absent (CI behavior).

## User Setup Required

None - no external service configuration required.

## Next Phase Readiness

- Complete 55-test BATS suite (01-04) passes on dev environment
- Phase 03 test framework is complete (all 3 plans done: infrastructure, filter unit tests, integration tests)
- Phase 04 (CI integration) can use `tests/bats/bin/bats tests/*.bats` as the test command
- Real-data integration tests (TFWK-04) require `ref/hs37d5.fa` in CI — either bundle or skip with HAVE_FULL_REF guard
- All expected output files in `tests/data/expected/` are committed and versioned as golden files

---
*Phase: 03-test-framework*
*Completed: 2026-02-18*
