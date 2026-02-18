# Phase 3: Test Framework - Context

**Gathered:** 2026-02-18
**Status:** Ready for planning

<domain>
## Phase Boundary

BATS test suite that verifies hardnormly.sh works correctly before any refactoring. Covers smoke tests (CLI args, help), filter-level verification (GATK and Freebayes), full pipeline integration, and regression diffs against expected output. All tests run against the current pre-refactor script.

</domain>

<decisions>
## Implementation Decisions

### Filter test precision
- Exact FILTER tag verification on specific variants (not count-based)
- Tests grouped by caller: one GATK test checks all GATK filter tags, one Freebayes test checks all Freebayes filter tags
- Both positive and negative assertions: verify target variants GET the expected tag AND control variants are NOT tagged
- Include synthetic variants that trigger multiple filters simultaneously; verify combined FILTER field (e.g., FILTER=DPu10het;QDlt2)

### Regression comparison strategy
- Exact body match, ignore headers: strip VCF ##header lines, then exact diff on data lines
- On regression failure: show count of differences plus first 20 changed lines (summary + first N)
- Regression tests cover GATK pipeline only; Freebayes is covered by filter-level tests
- Expected output file strategy: Claude's discretion — research best practices for committed vs generated expected files and fit to this repo

### Test environment & dependencies
- Real bcftools/bedtools required for all tests (no mocking)
- Skip gracefully with message when tools are missing; smoke tests (--help, arg validation) still run without tools
- Per-test temp directories via BATS setup/teardown; full isolation, no shared state between tests
- Tests must work on both WSL and native Linux; no Windows-path assumptions

### Failure diagnostics
- On filter test failure: show expected vs actual FILTER tags side by side
- Preserve temp directories on failure only; clean up on pass
- Log the actual bcftools/bedtools commands on failure for manual reproduction outside BATS
- No custom verbose mode; rely on BATS built-in --verbose-run

### Claude's Discretion
- Expected output file strategy (committed to repo vs generated) — research best practices
- BATS file organization (one file per test category vs monolithic)
- Helper function structure and shared BATS libraries
- Exact BATS assertion patterns and helper utilities

</decisions>

<specifics>
## Specific Ideas

No specific requirements — open to standard approaches

</specifics>

<deferred>
## Deferred Ideas

None — discussion stayed within phase scope

</deferred>

---

*Phase: 03-test-framework*
*Context gathered: 2026-02-18*
