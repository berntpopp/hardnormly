# Phase 3: Test Framework - Research

**Researched:** 2026-02-18
**Domain:** BATS (Bash Automated Testing System) — bash script integration testing, filter verification, regression testing
**Confidence:** HIGH

## Summary

Phase 3 builds a BATS test suite that verifies hardnormly.sh behavior against the pre-refactor codebase. Four test categories are required: smoke tests (CLI args), filter unit tests (GATK + Freebayes), integration test (full pipeline on real data), and regression tests (diff against committed expected output). All test data was prepared in Phase 2 and is committed to `tests/data/`.

BATS 1.x (bats-core) is the established standard. The `bats-assert` and `bats-support` libraries provide readable assertion functions that also print diagnostics on failure — critical for the "show expected vs actual FILTER tags" requirement. Installation via git submodules is the recommended pattern for library management, with `bats-core/bats-action` as the GitHub Actions integration point.

The key architecture decision is **one file per test category** (4 files) with shared helpers extracted to `tests/test_helper/`. This follows the BATS tutorial's recommended structure and keeps each file independently runnable. For expected output files, the "golden file" pattern applies: commit expected outputs to the repo, update them intentionally when behavior changes, and use exact diff for regression detection.

**Primary recommendation:** Install bats-core and bats-assert/bats-support as git submodules in `tests/bats/` and `tests/test_helper/`. Use `bats-core/bats-action@4.0.0` in CI. Organize tests as four `.bats` files with a `tests/test_helper/common.bash` shared helper.

## Standard Stack

### Core

| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| bats-core | 1.x (latest) | BATS test runner | The maintained fork of the original bats; official project |
| bats-assert | 2.1.0 | Assertion functions (assert_output, assert_line, etc.) | Prints diagnostic context on failure; essential for filter tag comparisons |
| bats-support | 0.3.0 | Failure output formatting used by bats-assert | Required dependency of bats-assert |

### Supporting

| Library | Version | Purpose | When to Use |
|---------|---------|---------|-------------|
| bats-file | 0.4.0 | File existence assertions | Verifying output VCF exists after pipeline runs |
| bats-core/bats-action | 4.0.0 | GitHub Actions installer for bats + all libraries | CI setup only |

### Alternatives Considered

| Instead of | Could Use | Tradeoff |
|------------|-----------|----------|
| Git submodules for libs | npm install @bats-core/bats | npm only has bats-core, not bats-assert/bats-support; submodules get everything |
| bats-core/bats-action | Manual apt-get install bats | apt bats may be pre-1.0 (outdated) on ubuntu-latest; bats-action gives current 1.x |
| bats-assert | Raw `[ "$status" -eq 0 ]` | Raw assertions don't print diagnostics on failure — the "show expected vs actual" requirement mandates bats-assert |

**Installation (git submodules — recommended):**
```bash
git submodule add https://github.com/bats-core/bats-core.git tests/bats
git submodule add https://github.com/bats-core/bats-support.git tests/test_helper/bats-support
git submodule add https://github.com/bats-core/bats-assert.git tests/test_helper/bats-assert
git submodule add https://github.com/bats-core/bats-file.git tests/test_helper/bats-file
```

**GitHub Actions CI:**
```yaml
- name: Setup Bats and bats libs
  id: setup-bats
  uses: bats-core/bats-action@4.0.0

- name: Run BATS tests
  env:
    BATS_LIB_PATH: ${{ steps.setup-bats.outputs.lib-path }}
  run: tests/bats/bin/bats tests/
```

**Local run:**
```bash
# After cloning with submodules:
git submodule update --init --recursive

# Run all tests:
tests/bats/bin/bats tests/*.bats

# Run one file:
tests/bats/bin/bats tests/01-smoke.bats

# Run with verbose output (shows stdout for passing tests too):
tests/bats/bin/bats --verbose-run tests/
```

## Architecture Patterns

### Recommended Project Structure

```
tests/
├── bats/                          # bats-core submodule
├── test_helper/
│   ├── bats-support/              # submodule (required by bats-assert)
│   ├── bats-assert/               # submodule (assertion functions)
│   ├── bats-file/                 # submodule (file existence assertions)
│   └── common.bash                # shared setup: tool checks, path helpers
├── data/                          # all test data (committed, from Phase 2)
│   ├── synthetic/                 # gatk_sample_A/B, freebayes_sample_A/B, mini_ref.fa
│   ├── real/                      # GIAB, 1kg, freebayes_tiny, chr22_16M_ref.fa
│   ├── expected/                  # golden files (committed expected outputs)
│   │   └── BASELINE.txt           # commit hash + per-variant filter summary
│   ├── include_regions.bed
│   ├── exclude_regions.bed
│   └── hg19_chr22.genome
├── 01-smoke.bats                  # TFWK-01: help, version, missing args
├── 02-gatk-filters.bats           # TFWK-02: GATK filter tag verification
├── 03-freebayes-filters.bats      # TFWK-03: Freebayes filter tag verification
├── 04-integration.bats            # TFWK-04: full pipeline + TFWK-05: regression
│                                  #          + TFWK-06: genome flag + TFWK-07: edge cases
└── generate_test_data.sh          # public data download script (Phase 2)
```

**Note on grouping TFWK-04 through TFWK-07:** Requirements 04, 05, 06, and 07 all require the full pipeline to run. Grouping them in `04-integration.bats` avoids duplicating the pipeline invocation overhead and keeps tool-availability checks in one place.

### Pattern 1: Shared Helper (`tests/test_helper/common.bash`)

**What:** Sourced in every `setup()` call. Provides: library loading, path computation, tool-availability check helper, and temp-dir-on-failure preservation.

**When to use:** Always — every test file sources this in `setup()`.

```bash
# tests/test_helper/common.bash
# Source: BATS tutorial pattern (bats-core.readthedocs.io/en/stable/tutorial.html)

_common_setup() {
    load 'test_helper/bats-support/load'
    load 'test_helper/bats-assert/load'
    load 'test_helper/bats-file/load'

    # Compute repo root from BATS_TEST_FILENAME so tests run from any directory
    REPO_ROOT="$(cd "$(dirname "$BATS_TEST_FILENAME")/.." >/dev/null 2>&1 && pwd)"
    export REPO_ROOT

    # Canonical paths to key test assets
    HARDNORMLY="${REPO_ROOT}/hardnormly.sh"
    TEST_DATA="${REPO_ROOT}/tests/data"
    SYNTH="${TEST_DATA}/synthetic"
    EXPECTED="${TEST_DATA}/expected"
    export HARDNORMLY TEST_DATA SYNTH EXPECTED

    # Per-test temp directory (BATS provides $BATS_TEST_TMPDIR)
    # On failure, preserve; on pass, BATS cleans up automatically
    TEST_TMP="${BATS_TEST_TMPDIR}"
    export TEST_TMP
}

# Check that bioinformatics tools are available; call from setup() in non-smoke tests
_require_tools() {
    local missing=()
    for tool in bcftools bedtools bgzip tabix; do
        command -v "$tool" >/dev/null 2>&1 || missing+=("$tool")
    done
    if [[ ${#missing[@]} -gt 0 ]]; then
        skip "Required tools not found: ${missing[*]}"
    fi
}

# Extract FILTER field for a specific position from a VCF
# Usage: _get_filter <vcf> <pos>
# Prints: the FILTER value (e.g., "DPu10het" or "PASS" or ".")
_get_filter() {
    local vcf="$1"
    local pos="$2"
    bcftools query -r "22:${pos}-${pos}" -f '%POS\t%FILTER\n' "$vcf" \
        | awk -v p="$pos" '$1==p {print $2}'
}
```

### Pattern 2: Smoke Tests (`tests/01-smoke.bats`)

**What:** Tests that run without bioinformatics tools. --help exits 0, --version prints version, missing required args exits non-zero.

**When to use:** These run in all environments, including CI without bioinformatics tools.

```bash
# tests/01-smoke.bats
setup() {
    load 'test_helper/common'
    _common_setup
    # No _require_tools call — smoke tests run without bcftools/bedtools
}

@test "--help exits 0" {
    run "$HARDNORMLY" --help
    assert_success
}

@test "--version prints version string" {
    run "$HARDNORMLY" --version
    assert_success
    assert_output --partial "Version:"
}

@test "missing -v arg exits non-zero" {
    run "$HARDNORMLY" -f /dev/null
    assert_failure
}

@test "missing -f arg exits non-zero" {
    run "$HARDNORMLY" -v /dev/null
    assert_failure
}
```

### Pattern 3: Filter Unit Tests (`tests/02-gatk-filters.bats`)

**What:** Run hardnormly.sh on gatk_sample_A with GATK filters. Assert each expected variant POS has the correct FILTER tag. Assert control variants are NOT tagged.

**Key technique:** Extract FILTER field per-position with `bcftools query`, then use `assert_equal` for exact match. Run the pipeline once in `setup_file()` (expensive), then per-test assertions use the cached output VCF.

```bash
# tests/02-gatk-filters.bats
# Source: bats-core tutorial — setup_file() runs once per file

setup_file() {
    load 'test_helper/common'
    _common_setup
    _require_tools

    # Run the full pipeline once; all tests in this file use the output
    OUTPUT_VCF="${BATS_FILE_TMPDIR}/gatk_filtered.vcf.gz"
    export OUTPUT_VCF

    run "$HARDNORMLY" \
        -v "${SYNTH}/gatk_sample_A.vcf.gz" \
        -f "${SYNTH}/mini_ref.fa" \
        -b "${TEST_DATA}/include_regions.bed" \
        -e "${TEST_DATA}/exclude_regions.bed" \
        --filters-file "${REPO_ROOT}/defaults/gatk_filters.txt" \
        -g "${TEST_DATA}/hg19_chr22.genome" \
        -o "$OUTPUT_VCF" \
        --auto-index
    assert_success
}

setup() {
    load 'test_helper/common'
    _common_setup
}

@test "gatkSNPhard: AS_FS>60 variant is tagged" {
    # POS 10200: fail gatkSNPhard (AS_FS=65.0)
    run _get_filter "$OUTPUT_VCF" 10200
    assert_output "gatkSNPhard"
}

@test "gatkSNPhard: passing variant is NOT tagged" {
    # POS 10100: passes all filters
    run _get_filter "$OUTPUT_VCF" 10100
    assert_output "PASS"
}

@test "DPu10het;gatkSNPhard: multi-filter variant has combined tag" {
    # POS that triggers both filters simultaneously
    run _get_filter "$OUTPUT_VCF" <POS>
    assert_output "DPu10het;gatkSNPhard"
}
```

**Note on multi-filter combined tags:** bcftools filter with `-m+` accumulates tags with semicolons. The order in the FILTER field follows filter application order (as defined in the filter file). Verify the exact combined tag format by running the pipeline manually once before writing tests.

### Pattern 4: Regression Tests (`tests/04-integration.bats`)

**What:** Run the full GATK pipeline on gatk_sample_A, strip `##` header lines, diff data lines against committed expected output. On failure, print count of differences and first 20 changed lines.

**Golden file strategy:** Expected outputs are committed to `tests/data/expected/`. To update them intentionally, run `tests/generate_test_data.sh --regenerate-expected` (a flag to be added in plan 02-04) and commit the diff.

```bash
@test "GATK regression: output matches expected" {
    local actual="${TEST_TMP}/gatk_actual.vcf"
    local expected_vcf="${EXPECTED}/gatk_sample_A_filtered.vcf.gz"

    # Run pipeline
    run "$HARDNORMLY" \
        -v "${SYNTH}/gatk_sample_A.vcf.gz" \
        -f "${SYNTH}/mini_ref.fa" \
        ...other flags... \
        -o "${TEST_TMP}/out.vcf.gz" --auto-index
    assert_success

    # Strip headers; compare data lines
    bcftools view -H "${TEST_TMP}/out.vcf.gz" > "$actual"
    bcftools view -H "$expected_vcf" > "${TEST_TMP}/expected_data.vcf"

    # On failure: count + first 20 lines
    if ! diff -u "${TEST_TMP}/expected_data.vcf" "$actual" > "${TEST_TMP}/diff.txt"; then
        local count
        count=$(grep -c '^[+-]' "${TEST_TMP}/diff.txt" || true)
        echo "Regression failure: ${count} changed lines"
        echo "First 20 differences:"
        head -20 "${TEST_TMP}/diff.txt"
        fail "Output does not match expected (see diff above)"
    fi
}
```

### Pattern 5: Per-Test Temp Dir and Conditional Cleanup

**What:** BATS provides `$BATS_TEST_TMPDIR` as a unique per-test directory. On pass, BATS cleans it up. On failure, BATS also cleans it up by default — BUT the `--no-tempdir-cleanup` flag preserves it for debugging.

**Decision:** Do NOT implement custom conditional cleanup in teardown. Instead:
1. Use `$BATS_TEST_TMPDIR` for all per-test output files
2. In CI, run `bats` normally (auto-cleanup on pass/fail is fine in CI)
3. For local debugging of failures, run with `bats --no-tempdir-cleanup tests/`

This avoids relying on `$BATS_TEST_COMPLETED` (which is undocumented/unstable) and is simpler.

```bash
teardown() {
    # No conditional cleanup needed.
    # BATS cleans BATS_TEST_TMPDIR automatically.
    # For debugging failures, run with: bats --no-tempdir-cleanup tests/
    :
}
```

**Note:** Use `$BATS_FILE_TMPDIR` for outputs shared across tests within a file (e.g., the pipeline output shared by all filter assertions in 02-gatk-filters.bats). `$BATS_FILE_TMPDIR` persists for the file lifetime and is cleaned up after `teardown_file()`.

### Pattern 6: Failure Diagnostics for Filter Tests

**What:** On FILTER tag mismatch, show expected vs actual side-by-side.

**How:** `bats-assert`'s `assert_equal` automatically shows expected vs actual on failure. For the command-in-filter tests, extract the FILTER column via `bcftools query` and compare with `assert_equal`.

```bash
@test "DPu10het: variant at pos 30100 is tagged" {
    local actual_filter
    actual_filter=$(_get_filter "$OUTPUT_VCF" 30100)

    # On failure, assert_equal prints: "expected: DPu10het | actual: PASS"
    assert_equal "DPu10het" "$actual_filter"
}
```

For the "log bcftools commands on failure" requirement: the `--verbose-run` flag on BATS prints all `run` output. The pipeline command is part of `setup_file()` output, which BATS captures and can show on failure with `--print-output-on-failure`.

### Pattern 7: Tool-Availability Skip with BATS

```bash
# In setup() for tests that require bioinformatics tools:
_require_tools() {
    local missing=()
    for tool in bcftools bedtools bgzip tabix; do
        command -v "$tool" >/dev/null 2>&1 || missing+=("$tool")
    done
    if [[ ${#missing[@]} -gt 0 ]]; then
        skip "Required tools not found: ${missing[*]}"
    fi
}
```

When `skip` is called in `setup()`, BATS marks the test as skipped (not failed) and prints the reason. The smoke tests (01-smoke.bats) do NOT call `_require_tools`, so they run even without tools.

### Pattern 8: load command vs bats_load_library

- Use `load 'test_helper/bats-support/load'` for **submodule-installed libraries** (relative path to test file)
- Use `bats_load_library bats-assert` + `BATS_LIB_PATH` for **system-installed libraries** (GitHub Actions via bats-action)
- The CI workflow must set `BATS_LIB_PATH` from `bats-action`'s output when using `bats_load_library`
- For maximum portability (submodules work both locally and in CI), use `load` with relative paths, and configure CI to run `git submodule update --init --recursive` in the checkout step

**Recommended: use `load` for submodules** (simpler, no BATS_LIB_PATH needed):

```bash
setup() {
    load 'test_helper/bats-support/load'
    load 'test_helper/bats-assert/load'
    load 'test_helper/bats-file/load'
}
```

**Alternative: use bats_load_library for system installs (GitHub Actions with bats-action)**:

```bash
setup() {
    bats_load_library bats-support
    bats_load_library bats-assert
    bats_load_library bats-file
}
```

Pick one approach and use it consistently. The submodule approach is self-contained (no CI configuration needed beyond checkout). Mixing both requires guards.

### Anti-Patterns to Avoid

- **Raw `[ "$status" -eq 0 ]`:** Doesn't print diagnostics on failure. Use `assert_success` instead.
- **`[ "$output" = "..." ]`:** Same problem. Use `assert_output` or `assert_output --partial`.
- **Running the pipeline in every `@test`:** For filter tests, run once in `setup_file()` and share the output. Running hardnormly.sh 20 times per filter test file is slow.
- **Piped commands inside `run`:** BATS does not capture pipeline exit codes correctly. Instead of `run bash -c 'cmd1 | cmd2'`, assign to variable first: `result=$(cmd1 | cmd2)` then assert on `$result`.
- **Hard-coded absolute paths in test files:** Use `REPO_ROOT` computed from `BATS_TEST_FILENAME` so tests run from any directory.
- **Checking `$BATS_TEST_COMPLETED` in teardown:** This variable is undocumented and subject to change. Use `--no-tempdir-cleanup` flag for debugging instead of custom teardown logic.
- **Using apt-get bats:** Ubuntu apt packages BATS 0.x (pre-1.0, old fork). Always use bats-core 1.x.
- **Committing expected output as plain text VCF:** The compressed `.vcf.gz` + `.tbi` format matches what hardnormly.sh actually produces. Store expected outputs in the same format.

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| Failure output | Custom error message printing | `bats-assert` assertion functions | assert_equal/assert_output print expected vs actual automatically on failure |
| File existence checks | `[ -f "$file" ]` with custom message | `assert_file_exists` (bats-file) | Prints full path and reason on failure |
| Multi-filter tag ordering | Custom sort/normalize function | Run pipeline and read exact output format | bcftools determines tag order; test the actual output, don't normalize it |
| FILTER tag extraction | Complex awk/grep script | `bcftools query -r "22:POS-POS" -f '%FILTER\n'` | One-liner, handles compressed VCF natively |
| Running bats in CI | Custom bash wrapper | `bats-core/bats-action@4.0.0` | Installs bats + all libs, sets BATS_LIB_PATH in one step |

**Key insight:** bats-assert's diagnostic output is the main reason to use it over raw bash assertions. The "show expected vs actual FILTER tags side by side" requirement is handled automatically — no custom helper needed.

## Common Pitfalls

### Pitfall 1: setup_file() vs setup() for Pipeline Invocations

**What goes wrong:** Putting `hardnormly.sh` invocation in `setup()` runs it N times (once per test). For filter tests with 10+ `@test` blocks, this multiplies runtime by 10x.
**Why it happens:** `setup()` runs before each test; `setup_file()` runs once per file.
**How to avoid:** Run the pipeline in `setup_file()`, export `OUTPUT_VCF` to `$BATS_FILE_TMPDIR`. Per-test assertions read from this shared file.
**Warning signs:** BATS test suite takes 10+ minutes; `ps aux` shows 15 hardnormly.sh processes.

### Pitfall 2: Temp Directory Lifetime Mismatch

**What goes wrong:** Storing pipeline output in `$BATS_TEST_TMPDIR` inside `setup_file()`. `BATS_TEST_TMPDIR` is per-test and gets destroyed after each test — the output disappears before subsequent tests can read it.
**Why it happens:** `setup_file()` runs in a test context where BATS_TEST_TMPDIR is the FIRST test's temp dir, which is cleaned up after the first test.
**How to avoid:** Use `$BATS_FILE_TMPDIR` for output shared across tests within a file. `BATS_FILE_TMPDIR` persists for the entire file's lifetime.
**Warning signs:** Second test in the file fails with "file not found" on the pipeline output.

### Pitfall 3: load Path Resolution

**What goes wrong:** `load 'test_helper/bats-support/load'` fails with "no such file" when tests are run from a different working directory.
**Why it happens:** `load` resolves relative paths relative to the test file's directory (`$BATS_TEST_DIRNAME`), not the working directory.
**How to avoid:** BATS resolves `load` paths relative to the test file, so `load 'test_helper/bats-support/load'` works correctly as long as the submodule is at `tests/test_helper/bats-support/`. No changes needed for running from different directories.
**Warning signs:** Tests pass when run from `tests/` but fail when run from repo root as `bats tests/`.

### Pitfall 4: bcftools query FILTER Output Format

**What goes wrong:** `bcftools query -f '%FILTER\n'` returns `.` for variants with no filter applied (not "PASS"). A test asserting `FILTER == "PASS"` fails for unflagged variants.
**Why it happens:** VCF format distinguishes between "no filter applied" (`.`) and "explicitly passed" (`PASS`). hardnormly.sh may output either depending on its bcftools view invocation.
**How to avoid:** Verify the actual FILTER value for a known-passing variant by running `bcftools query -f '%POS\t%FILTER\n' expected/gatk_sample_A_filtered.vcf.gz` before writing assertions. Adjust assert to match actual output (either `.` or `PASS`).
**Warning signs:** All "control variant passes filter" assertions fail even though the pipeline ran correctly.

### Pitfall 5: Multi-Filter Tag Order

**What goes wrong:** A variant triggers both DPu10het and gatkSNPhard. Test asserts `FILTER == "DPu10het;gatkSNPhard"` but actual output is `FILTER == "gatkSNPhard;DPu10het"` (different order).
**Why it happens:** bcftools accumulates filter tags in application order. The filter file defines the order, but assertion writers may guess wrong.
**How to avoid:** Before writing multi-filter assertions, run the pipeline on the synthetic VCF and inspect the actual FILTER field: `bcftools query -f '%POS\t%FILTER\n' output.vcf.gz`. Use that exact string in assertions.
**Warning signs:** Multi-filter tests fail with "expected: DPu10het;gatkSNPhard | actual: gatkSNPhard;DPu10het".

### Pitfall 6: Expected Output Staleness

**What goes wrong:** Regression test compares current output against expected output generated at a different commit. Any pipeline change (including bug fixes) breaks regression tests.
**Why it happens:** Expected outputs are a snapshot of hardnormly.sh behavior at one point in time.
**How to avoid:** The "golden file" discipline: when pipeline behavior intentionally changes, regenerate expected outputs by running hardnormly.sh on the synthetic inputs and committing the new outputs. The BASELINE.txt records the commit hash for traceability. Never silently accept a regression failure — investigate it.
**Warning signs:** Regression tests fail after a pipeline commit; diff shows valid-looking output that just differs from the expected snapshot.

### Pitfall 7: WSL vs Native Linux Paths

**What goes wrong:** Tests use Windows-style paths (`C:/development/...`) or WSL paths (`/mnt/c/...`) that don't exist in native Linux CI.
**Why it happens:** Development on Windows via WSL introduces path assumptions.
**How to avoid:** Always compute paths from `BATS_TEST_FILENAME` using `$(dirname ...)` expansion. Never hard-code host paths. This is already enforced by the `REPO_ROOT` pattern in Pattern 1.
**Warning signs:** Tests pass locally (WSL) but fail in GitHub Actions with "no such file or directory".

### Pitfall 8: Expected Output Files Not Yet Generated (Phase 2 dependency)

**What goes wrong:** Regression tests reference `tests/data/expected/gatk_sample_A_filtered.vcf.gz` which doesn't exist because Phase 2 plan 04 hasn't completed.
**Why it happens:** Plans 02-03 and 02-04 are not committed yet (per current git status). Phase 3 planning must account for this dependency.
**How to avoid:** Phase 3 plan 03-04 (regression tests) is blocked on Phase 2 plan 02-04. The planner must sequence 02-04 as a prerequisite for 03-04.

## Code Examples

Verified patterns from official sources:

### Setup Function with Common Helper

```bash
# Source: bats-core tutorial (bats-core.readthedocs.io/en/stable/tutorial.html)

setup() {
    load 'test_helper/bats-support/load'
    load 'test_helper/bats-assert/load'
    load 'test_helper/bats-file/load'
    _common_setup
    _require_tools  # skip if bcftools/bedtools missing (omit in 01-smoke.bats)
}
```

### Running a Command and Checking Exit/Output

```bash
# Source: bats-core writing-tests docs

@test "--help exits 0" {
    run "$HARDNORMLY" --help
    assert_success
    assert_output --partial "Usage:"
}

@test "missing required arg exits non-zero" {
    run "$HARDNORMLY" -f /dev/null   # missing -v
    assert_failure
}
```

### FILTER Tag Assertion Pattern

```bash
# Source: derived from bats-assert docs + bcftools query usage

@test "DPu10het: DP<10 het variant is tagged" {
    local actual
    actual=$(bcftools query -r "22:30100-30100" \
             -f '%FILTER\n' "$OUTPUT_VCF")
    assert_equal "DPu10het" "$actual"
}

@test "DPu10het: DP>=10 variant is not tagged" {
    local actual
    actual=$(bcftools query -r "22:10100-10100" \
             -f '%FILTER\n' "$OUTPUT_VCF")
    # Assert it is either "." or "PASS" (not a filter tag)
    refute_line --partial "DPu10het"
    # Or if exact value is known:
    assert_equal "PASS" "$actual"
}
```

### Regression Diff with Diagnostics

```bash
# Source: golden file testing pattern, adapted for BATS

@test "GATK regression: data lines match expected" {
    local actual="${TEST_TMP}/actual_data.vcf"
    local expected="${TEST_TMP}/expected_data.vcf"
    local diff_out="${TEST_TMP}/diff.txt"

    bcftools view -H "${TEST_TMP}/output.vcf.gz" > "$actual"
    bcftools view -H "${EXPECTED}/gatk_sample_A_filtered.vcf.gz" > "$expected"

    if ! diff -u "$expected" "$actual" > "$diff_out" 2>&1; then
        local changed_lines
        changed_lines=$(grep -c '^[+-]' "$diff_out" || true)
        echo "Regression failure: ${changed_lines} changed lines"
        echo "=== First 20 differences ==="
        head -20 "$diff_out"
        fail "Output does not match expected output"
    fi
}
```

### GitHub Actions CI Integration

```yaml
# Source: bats-core/bats-action GitHub README

jobs:
  test:
    name: Test
    runs-on: ubuntu-latest
    steps:
      - name: Checkout (with submodules)
        uses: actions/checkout@v4
        with:
          submodules: recursive

      - name: Install bioinformatics tools
        run: sudo apt-get update && sudo apt-get install -y bcftools bedtools

      - name: Setup Bats and bats libs
        id: setup-bats
        uses: bats-core/bats-action@4.0.0

      - name: Run BATS tests
        env:
          BATS_LIB_PATH: ${{ steps.setup-bats.outputs.lib-path }}
        run: tests/bats/bin/bats --print-output-on-failure tests/
```

**Note:** If using git submodules for libraries (recommended), replace `bats_load_library` with `load` in test files and remove `BATS_LIB_PATH`. The submodule approach doesn't need bats-action's lib-path — it just needs `submodules: recursive` in the checkout step. In that case, replace `tests/bats/bin/bats` with `bash tests/bats/bin/bats` or add it to PATH.

**Simpler CI with submodules only:**
```yaml
- name: Checkout (with submodules)
  uses: actions/checkout@v4
  with:
    submodules: recursive

- name: Install bcftools + bedtools + bats
  run: sudo apt-get update && sudo apt-get install -y bcftools bedtools

- name: Run BATS tests
  run: tests/bats/bin/bats --print-output-on-failure tests/*.bats
```

### Makefile Target for Tests

```makefile
## Run BATS tests
test:
	tests/bats/bin/bats --print-output-on-failure tests/*.bats

## Run BATS tests preserving temp dirs (for debugging failures)
test-debug:
	tests/bats/bin/bats --no-tempdir-cleanup --verbose-run tests/*.bats
```

## Expected Output File Strategy (Claude's Discretion)

**Recommendation: Commit expected outputs as golden files.**

Rationale:
1. **Tests work immediately after clone** — no "generate expected outputs" step needed before running tests.
2. **Diffs are visible in git** — when expected outputs change, `git diff tests/data/expected/` shows exactly what changed. Reviewers can judge whether the change is intentional.
3. **Industry standard for this pattern** — "golden files" are the established practice for regression testing CLI tools that produce file output. Languages from Go to Rust have golden file libraries. The consensus is: commit golden files, update them intentionally.
4. **Matches Phase 2's design** — Phase 2 plan 04 already generates expected outputs and commits them with BASELINE.txt.

**Update workflow:** When hardnormly.sh behavior intentionally changes, run:
```bash
# Regenerate expected outputs (using committed mini_ref.fa and synthetic VCFs)
./hardnormly.sh -v tests/data/synthetic/gatk_sample_A.vcf.gz \
    -f tests/data/synthetic/mini_ref.fa \
    -b tests/data/include_regions.bed \
    -e tests/data/exclude_regions.bed \
    --filters-file defaults/gatk_filters.txt \
    -g tests/data/hg19_chr22.genome \
    -o tests/data/expected/gatk_sample_A_filtered.vcf.gz \
    --auto-index
# Update BASELINE.txt with new commit hash
echo "# Regenerated at commit: $(git rev-parse HEAD)" >> tests/data/expected/BASELINE.txt
# Review the diff, then commit
git diff tests/data/expected/
git add tests/data/expected/ && git commit -m "test: update expected outputs for pipeline change"
```

**What to compare:** Strip `##` header lines; diff data lines only. Header lines contain timestamps and tool versions that change on every run. Data lines (CHROM, POS, FILTER, FORMAT, sample columns) are deterministic.

## BATS File Organization Decision (Claude's Discretion)

**Recommendation: 4 separate files** (one per test category), NOT a monolithic file.

| File | Requirements | Rationale |
|------|-------------|-----------|
| `01-smoke.bats` | TFWK-01 | Runs without tools; distinct tool-availability behavior |
| `02-gatk-filters.bats` | TFWK-02 | Uses `setup_file()` pattern; one pipeline run for all assertions |
| `03-freebayes-filters.bats` | TFWK-03 | Same pattern as GATK; independent filter set |
| `04-integration.bats` | TFWK-04, 05, 06, 07 | Full pipeline + real data + regression + edge cases + genome flag |

Running `bats tests/` or `bats tests/*.bats` executes all four. Running `bats tests/01-smoke.bats` runs only smoke tests — useful when bioinformatics tools aren't available.

**Why not monolithic:** A single file cannot use `setup_file()` to share pipeline outputs between categories that need different pipeline invocations (GATK vs Freebayes). The filter tests would interfere with each other.

**Why not 7 files (one per requirement):** Over-splitting creates coordination overhead. TFWK-04 through TFWK-07 all require the same full-pipeline invocation pattern and tool availability — grouping them in one file avoids duplication.

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| Pre-1.0 bats (sstephenson) | bats-core 1.x | 2018 (fork) | `setup_file()`/`teardown_file()`, `BATS_TEST_TMPDIR`, better isolation |
| `[ "$status" -eq 0 ]` | `assert_success` (bats-assert) | With bats-assert library | Failure output shows actual vs expected automatically |
| apt install bats | git submodule or bats-action | With 1.x release | Ensures 1.x (not outdated 0.x) |
| `--no-tempdir-cleanup` was not available | Available since bats-core 1.x | bats-core 1.x | Debug temp dirs without custom teardown |
| No `setup_file()` | `setup_file()` / `teardown_file()` | bats-core 1.x | Expensive setup runs once per file, not per test |

**Deprecated/outdated:**
- `sstephenson/bats`: Archived, unmaintained. Use `bats-core/bats-core`.
- `apt install bats` on Ubuntu 22.04 (jammy): installs bats 1.2.1 (acceptable) but Ubuntu 20.04 installs pre-1.0. Use submodules or bats-action for guaranteed version.
- Manual `BATS_TMPDIR` management: replaced by `BATS_TEST_TMPDIR` (per-test) and `BATS_FILE_TMPDIR` (per-file) in 1.x.

## Open Questions

1. **Expected outputs not yet generated (Phase 2 plan 04 in progress)**
   - What we know: `tests/data/expected/` directory exists but is empty (`.gitkeep` only). Plan 02-04 generates expected outputs by running hardnormly.sh.
   - What's unclear: Whether plan 02-04 will complete before Phase 3 planning begins.
   - Recommendation: Plan 03-04 (regression tests) must list 02-04 as a dependency. If 02-04 hasn't run yet, the regression test plan must include a task to generate expected outputs first.

2. **Exact FILTER values for "passing" variants**
   - What we know: hardnormly.sh may output `PASS` or `.` for variants that pass all filters, depending on how bcftools view is invoked.
   - What's unclear: Which value hardnormly.sh currently produces for passing variants.
   - Recommendation: During plan 03-02 (GATK filter tests), the first task should be to inspect the actual expected output: `bcftools query -f '%POS\t%FILTER\n' tests/data/expected/gatk_sample_A_filtered.vcf.gz`. Use the actual value in assertions.

3. **Multi-filter combined tag order**
   - What we know: The synthetic VCFs include variants designed to trigger multiple filters (noted in CONTEXT.md and Phase 2 summaries). bcftools filter with `-m+` appends tags in application order.
   - What's unclear: The exact combined tag strings (e.g., is it `DPu10het;gatkSNPhard` or `gatkSNPhard;DPu10het`?). This depends on the order of filters in `defaults/gatk_filters.txt`.
   - Recommendation: `defaults/gatk_filters.txt` lists DPu10het before gatkSNPhard. So combined tags should be `DPu10het;gatkSNPhard`. Verify by inspecting the expected output VCF.

4. **bcftools + bedtools availability in GitHub Actions**
   - What we know: `sudo apt-get install -y bcftools bedtools` works on ubuntu-latest. Version available: bcftools 1.13+ (not 1.20 from conda, but functionally equivalent for these tests).
   - What's unclear: Whether any 1.20-specific bcftools features are required by the tests.
   - Recommendation: Use apt-get for CI (bcftools + bedtools). The hardnormly.sh tests only require standard bcftools features available in 1.13+. No conda setup needed in CI for the test job.

5. **Git submodules vs bats-action library loading method**
   - What we know: `load 'test_helper/bats-support/load'` works for submodule installs; `bats_load_library bats-support` works for bats-action installs. They cannot be mixed without guards.
   - What's unclear: The preferred approach for this specific repo.
   - Recommendation: Use git submodules + `load` for self-contained testing. Avoid bats-action's library system. The submodule approach requires `actions/checkout@v4` with `submodules: recursive` — simpler than managing `BATS_LIB_PATH`.

## Sources

### Primary (HIGH confidence)
- bats-core official docs: https://bats-core.readthedocs.io/en/stable/writing-tests.html — setup/teardown lifecycle, BATS_TEST_TMPDIR, BATS_FILE_TMPDIR, skip syntax, run command, load vs bats_load_library
- bats-core tutorial: https://bats-core.readthedocs.io/en/stable/tutorial.html — multi-file structure, shared helper pattern, submodule installation
- bats-core usage docs: https://bats-core.readthedocs.io/en/stable/usage.html — --verbose-run, --no-tempdir-cleanup, --print-output-on-failure flags
- bats-assert README: https://github.com/bats-core/bats-assert/blob/master/README.md — assert_output, assert_line, assert_equal, assert_success, assert_failure, refute_output, matching modes
- bats-file README: https://github.com/bats-core/bats-file/blob/master/README.md — assert_file_exists, assert_file_empty, assert_files_equal
- bats-core/bats-action GitHub: https://github.com/bats-core/bats-action — YAML usage, installed library versions, BATS_LIB_PATH output

### Secondary (MEDIUM confidence)
- BATS_TEST_COMPLETED behavior: https://github.com/bats-core/bats-core/issues/383 — BATS_TEST_COMPLETED is undocumented and unstable; --no-tempdir-cleanup is the recommended alternative
- Golden file testing: https://ro-che.info/articles/2017-12-04-golden-tests — rationale for committing expected outputs vs generating them; consensus: commit them
- Pytest golden file update workflow: https://johal.in/pytest-regressions-data-golden-file-updates-2025/ — update workflow pattern applicable across languages

### Tertiary (LOW confidence)
- WebSearch finding: apt bats on ubuntu-latest provides recent enough version (1.2+) for basic use, but not guaranteed to be 1.x current — prefer submodules for version control

## Metadata

**Confidence breakdown:**
- Standard stack (bats-core, bats-assert, bats-support): HIGH — official docs verified, versions confirmed from bats-action defaults
- Architecture patterns (file organization, setup_file, load): HIGH — from official tutorial and writing-tests docs
- Expected output strategy (golden files): HIGH — well-established industry pattern, corroborated by multiple sources
- BATS_TEST_COMPLETED behavior: MEDIUM — GitHub issue confirmed it's undocumented; --no-tempdir-cleanup alternative confirmed from usage docs
- Pitfall specifics (FILTER values, tag order): MEDIUM — derived from bcftools behavior and filter file order; must be verified against actual expected outputs

**Research date:** 2026-02-18
**Valid until:** 2026-08-18 (bats-core is stable; assertion API changes rarely)
