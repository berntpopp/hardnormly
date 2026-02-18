---
phase: 05-features-docs
verified: 2026-02-18T20:56:34Z
status: gaps_found
score: 3/5 must-haves verified
gaps:
  - truth: hardnormly.sh --help shows usage examples and filter file format description
    status: failed
    reason: >
      --help shows an Example: section (DOCS-01 satisfied) but does not describe the
      filter file format (DOCS-02 missing). Format is in README.md only. The pre-existing
      unit test show_help outputs options section (tests/06-lib-cli.bats:30) now fails
      because show_help was changed from Options: to Options (run-pipeline): in plan
      05-03 without updating the test.
    artifacts:
      - path: lib/cli.sh
        issue: >-
          show_help uses Options (run-pipeline): not Options: - breaks unit test at
          tests/06-lib-cli.bats line 30
      - path: lib/cli.sh
        issue: >-
          show_help has no filter file format description (DOCS-02 requires this in --help)
    missing:
      - Filter file format section in show_help: three-column format filter_name action bcftools_expression
      - Fix tests/06-lib-cli.bats line 30: update assertion from Options: to Options (run-pipeline):
  - truth: generate-inclusion-bed and generate-exclusion-bed produce valid merged BED files
    status: partial
    reason: >
      Both handlers fully implemented in hardnormly.sh (not stubs). They call lib/bed.sh
      correctly. No automated BATS tests verify actual BED file output - only dispatcher
      routing is smoke-tested.
    artifacts:
      - path: hardnormly.sh
        issue: >-
          cmd_generate_inclusion_bed and cmd_generate_exclusion_bed implemented but untested end-to-end
    missing:
      - Integration test: generate-inclusion-bed with real BED+genome produces non-empty sorted BED output
      - Integration test: generate-exclusion-bed with real BED produces non-empty BED output
---

# Phase 5: Features and Docs Verification Report

**Phase Goal:** New user-facing capabilities work and are discoverable via help text
**Verified:** 2026-02-18T20:56:34Z
**Status:** gaps_found
**Re-verification:** No - initial verification

## Goal Achievement

### Observable Truths

| # | Truth | Status | Evidence |
|---|-------|--------|----------|
| 1 | run-pipeline and legacy flag invocation produce identical output | VERIFIED | Dispatcher routes both to parse_args (hardnormly.sh lines 219-241) |
| 2 | generate-inclusion-bed and generate-exclusion-bed produce valid merged BED files | PARTIAL | Handlers fully implemented; no integration tests confirm actual BED output |
| 3 | --caller gatk uses gatk_filters.txt without --filters-file | VERIFIED | Implemented at hardnormly.sh lines 255-272 with conflict warning |
| 4 | plot-vcfstats failure is logged but pipeline continues | VERIFIED | Non-fatal pattern at hardnormly.sh lines 428-435 |
| 5 | --help shows usage examples and filter file format description | FAILED | Example: section present (DOCS-01 OK); filter format missing from --help (DOCS-02 gap); unit test regression at tests/06-lib-cli.bats:30 |

**Score:** 3/5 truths verified

### Required Artifacts

| Artifact | Expected | Status | Details |
|----------|----------|--------|---------|
| hardnormly.sh | Subcommand dispatcher | VERIFIED | 443 lines; dispatcher at lines 216-248 |
| hardnormly.sh | Non-fatal plot-vcfstats | VERIFIED | Lines 425-435: _plot_rc pattern with SC2310 disable |
| hardnormly.sh | --caller flag resolves to filter file | VERIFIED | Lines 254-272: gatk and freebayes resolution |
| hardnormly.sh | cmd_generate_inclusion_bed (full implementation) | VERIFIED | Lines 60-143: arg parser, validation, calls normalize_bed + merge_include_beds |
| hardnormly.sh | cmd_generate_exclusion_bed (full implementation) | VERIFIED | Lines 145-214: arg parser, validation, calls normalize_bed + merge_exclude_beds |
| lib/cli.sh | --caller and --strip-annotations in parse_args | VERIFIED | Lines 227-241: both flags with value validation |
| lib/cli.sh | Compact show_help with Example: section | PARTIAL | Example: present at line 46; filter format missing; Options (run-pipeline): label breaks unit test |
| lib/cli.sh | show_help_generate_inclusion_bed exits 0 | VERIFIED | Lines 62-78 |
| lib/cli.sh | show_help_generate_exclusion_bed exits 0 | VERIFIED | Lines 81-95 |
| lib/annotate.sh | strip_vcf_annotations function | VERIFIED | Lines 27-34: bcftools annotate -x strip_list |
| lib/stats.sh | plot_stats_output returns 1 on failure | VERIFIED | Lines 29-38: returns 1; hardnormly.sh handles non-fatally |
| scripts/generate_exclusion_bed.sh | Downloads 4 sources, produces merged BED | VERIFIED | 230 lines; ENCODE Blacklist, SuperDups, RepeatMasker LCR, UCSC gap/cytoBandIdeo |
| defaults/gatk_filters.txt | Exists for --caller resolution | VERIFIED | File exists |
| defaults/freebayes_filters.txt | Exists for --caller resolution | VERIFIED | File exists |

### Key Link Verification

| From | To | Via | Status | Details |
|------|----|-----|--------|---------|
| hardnormly.sh dispatcher | parse_args (run-pipeline) | case run-pipeline then parse_args | WIRED | Lines 219-221 |
| hardnormly.sh dispatcher | parse_args (legacy -* flags) | case -* then parse_args | WIRED | Lines 239-241 |
| hardnormly.sh --caller block | defaults/gatk_filters.txt | _SCRIPT_DIR/defaults/gatk_filters.txt | WIRED | Line 258 |
| hardnormly.sh --caller block | defaults/freebayes_filters.txt | _SCRIPT_DIR/defaults/freebayes_filters.txt | WIRED | Line 261 |
| hardnormly.sh step 4.5 | strip_vcf_annotations | if -n strip_annotations guard | WIRED | Lines 381-386 |
| hardnormly.sh plot block | plot_stats_output non-fatal | if ! plot_stats_output; then _plot_rc=1 | WIRED | Lines 429-435 |
| cmd_generate_inclusion_bed | normalize_bed + merge_include_beds | direct lib/bed.sh function calls | WIRED | Lines 129-137 |
| cmd_generate_exclusion_bed | normalize_bed + merge_exclude_beds | direct lib/bed.sh function calls | WIRED | Lines 200-208 |

### Requirements Coverage

| Requirement | Status | Blocking Issue |
|-------------|--------|----------------|
| FEAT-01: Non-fatal plot-vcfstats | SATISFIED | hardnormly.sh lines 425-435 |
| FEAT-02: Subcommand dispatcher | SATISFIED | hardnormly.sh lines 216-248 |
| FEAT-03: generate-inclusion-bed subcommand | SATISFIED | hardnormly.sh lines 60-143 |
| FEAT-04: generate-exclusion-bed subcommand | SATISFIED | hardnormly.sh lines 145-214 |
| FEAT-05: Helper script for exclusion BED from public sources | SATISFIED | scripts/generate_exclusion_bed.sh 230 lines |
| FEAT-06: --caller auto-selects filter file | SATISFIED | hardnormly.sh lines 254-272 |
| FEAT-07: No-subcommand defaults to run-pipeline | SATISFIED | case -* routes to parse_args |
| FEAT-08: --strip-annotations flag | SATISFIED | parse_args + step 4.5 in hardnormly.sh |
| DOCS-01: --help includes usage examples | SATISFIED | Example: section in show_help (lib/cli.sh line 46) |
| DOCS-02: --help describes filter file format | BLOCKED | Filter format in README.md only, not in --help output |
| DOCS-03: Exclusion BED sources documented | SATISFIED | README.md lines 161-177; scripts/generate_exclusion_bed.sh header |

### Anti-Patterns Found

| File | Line | Pattern | Severity | Impact |
|------|------|---------|----------|--------|
| tests/06-lib-cli.bats | 30 | assert_output --partial Options: fails against Options (run-pipeline): | BLOCKER | Test show_help outputs options section FAILS |
| lib/cli.sh | 21 | Options (run-pipeline): introduced in 05-03 without updating the phase-4 unit test | BLOCKER | 1 of 118 BATS tests now fails |

### Human Verification Required

#### 1. generate-inclusion-bed produces valid BED output

**Test:** Run hardnormly.sh generate-inclusion-bed -b tests/data/ref/target_region.bed -g defaults/hg19.genome -o /tmp/test_inc.bed with bioinformatics tools active
**Expected:** /tmp/test_inc.bed is a non-empty valid 3-column BED file with sorted coordinates
**Why human:** Test environment lacks bcftools/bedtools; all tool-dependent BATS tests are skipped

#### 2. generate-exclusion-bed produces valid BED output

**Test:** Run hardnormly.sh generate-exclusion-bed -e tests/data/ref/exclusion_region.bed -o /tmp/test_exc.bed
**Expected:** /tmp/test_exc.bed is a non-empty valid 3-column BED file
**Why human:** Test environment lacks bioinformatics tools

#### 3. --caller gatk wires end-to-end

**Test:** Run full pipeline with --caller gatk on a test VCF, without --filters-file. Check FILTER tags in output.
**Expected:** Output VCF has FILTER tags matching gatk_filters.txt entries (gatkSNPhard, gatkINDELhard, etc.)
**Why human:** Full pipeline requires bioinformatics tools

#### 4. plot-vcfstats failure is non-fatal in real execution

**Test:** Run pipeline with --generate-stats --plot-stats --plot-output-dir pointing to a write-protected directory so plot-vcfstats fails. Verify exit 0 and warning in log.
**Expected:** Pipeline exits 0; log contains Warning: plot-vcfstats failed; pipeline continues.
**Why human:** Reaching the plot step requires bioinformatics tools to complete all prior pipeline steps

### Gaps Summary

Two gaps block full goal achievement:

**Gap 1 - DOCS-02 missing from --help (Critical):** Success criterion 5 requires --help to show the filter file format description. The three-column format (filter_name action bcftools_expression) is documented in README.md lines 85-103 but NOT in --help output. The help mentions --filters and --filters-file with one-line descriptions but does not explain the file format. Fix: add a brief filter file format section to show_help in lib/cli.sh.

**Gap 2 - Unit test regression from show_help rewrite (Critical):** Plan 05-03 changed the label from Options: to Options (run-pipeline): in show_help. The existing unit test at tests/06-lib-cli.bats line 30 asserts --partial Options: which does not match Options (run-pipeline): as a substring (no colon adjacent to Options in the new label). This is a test regression. Fix: update the test assertion to --partial Options (run-pipeline): or to --partial run-pipeline.

Both gaps are targeted fixes. All core features - subcommand dispatcher, --caller, --strip-annotations, non-fatal plot, generate-*-bed handlers, and scripts/generate_exclusion_bed.sh - are correctly implemented and wired.

---

_Verified: 2026-02-18T20:56:34Z_
_Verifier: Claude (gsd-verifier)_
