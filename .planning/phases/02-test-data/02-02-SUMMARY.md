---
phase: 02-test-data
plan: 02
subsystem: testing
tags: [bcftools, bgzip, tabix, vcf, synthetic-data, gatk, freebayes, fill-tags, vcf-filters]

# Dependency graph
requires:
  - phase: 02-01
    provides: mini_ref.fa (all-N contig 22, 100001 bases) and indexed FASTA for bcftools norm
provides:
  - tests/data/synthetic/gatk_sample_A.vcf.gz - 16-variant GATK VCF with all 7 filter conditions covered
  - tests/data/synthetic/gatk_sample_B.vcf.gz - 8-variant GATK VCF for batch testing
  - tests/data/synthetic/freebayes_sample_A.vcf.gz - 12-variant Freebayes VCF with all 9 filter conditions covered
  - tests/data/synthetic/freebayes_sample_B.vcf.gz - 7-variant Freebayes VCF for batch testing
  - tabix .tbi indexes for all 4 VCF files
affects:
  - 02-03 (expected outputs generated from these VCFs using filter files)
  - 02-04 (BATS tests validate filters against these VCFs)
  - Phase 3 filter unit tests (core test data)

# Tech tracking
tech-stack:
  added: [bgzip (BGZF compression), tabix (VCF indexing), bcftools +fill-tags (VAF computation)]
  patterns:
    - Write .vcf plaintext then bgzip (never write .vcf.gz directly)
    - REF=N for all variants when reference contig is all-N (avoids REF_MISMATCH in norm --check-ref w)
    - FORMAT/AD with Number=R enables fill-tags VAF computation without declaring FORMAT/VAF in header
    - Freebayes strand/read-position fields as INFO (site-level) not FORMAT

key-files:
  created:
    - tests/data/synthetic/gatk_sample_A.vcf.gz
    - tests/data/synthetic/gatk_sample_A.vcf.gz.tbi
    - tests/data/synthetic/gatk_sample_B.vcf.gz
    - tests/data/synthetic/gatk_sample_B.vcf.gz.tbi
    - tests/data/synthetic/freebayes_sample_A.vcf.gz
    - tests/data/synthetic/freebayes_sample_A.vcf.gz.tbi
    - tests/data/synthetic/freebayes_sample_B.vcf.gz
    - tests/data/synthetic/freebayes_sample_B.vcf.gz.tbi
  modified: []

key-decisions:
  - "REF=N for all variants - mini_ref.fa is entirely N (chr22:1-100001 is telomeric); bcftools TYPE() correctly classifies N>A as SNP and N>NA as INDEL"
  - "Freebayes AO/SAF/SAR/RPR/RPL as INFO fields (Type=Integer) matching filter expression syntax INFO/SAF==0"
  - "FORMAT/VAF not declared in header - bcftools +fill-tags computes it from FORMAT/AD at test time"
  - "WSL with /home/bernt/miniforge3/envs/hardnormly/bin in PATH provides bcftools 1.20 on Windows"

patterns-established:
  - "Synthetic VCF pattern: write .vcf with tabs, bgzip -f, bcftools index -t, bcftools view validate"
  - "N-reference pattern: use REF=N when reference is all-N to avoid REF_MISMATCH warnings"
  - "Filter coverage pattern: one pass + one fail variant per filter condition (not just per filter)"

# Metrics
duration: 45min
completed: 2026-02-18
---

# Phase 2 Plan 02: Synthetic VCF Creation Summary

**4 bgzip-indexed VCFs (43 total variants) covering all 7 GATK and 9 Freebayes filter conditions with REF=N to match all-N mini_ref.fa contig**

## Performance

- **Duration:** ~45 min
- **Started:** 2026-02-18T (session start)
- **Completed:** 2026-02-18
- **Tasks:** 2
- **Files modified:** 8 created (4 .vcf.gz + 4 .tbi)

## Accomplishments

- Created gatk_sample_A.vcf.gz with 16 variants covering all 7 GATK filter conditions (each sub-condition of gatkSNPhard and gatkINDELhard has its own triggering variant)
- Created freebayes_sample_A.vcf.gz with 12 variants covering all 9 Freebayes filter conditions including SAF=0 and SAR=0 as separate cases for strandBias
- All 4 VCFs: zero REF mismatches against mini_ref.fa, FORMAT/VAF computable for all variants via bcftools +fill-tags

## Task Commits

Each task was committed atomically:

1. **Task 1: Create synthetic GATK VCFs (samples A and B)** - `eba15ec` (feat)
2. **Task 2: Create synthetic Freebayes VCFs (samples A and B)** - `63af0f8` (feat)

**Plan metadata:** (see final commit below)

## Files Created/Modified

- `tests/data/synthetic/gatk_sample_A.vcf.gz` - 16 GATK variants: 7 SNP filter-fail scenarios + 4 INDEL fail scenarios + 5 DP/VAF fail scenarios; REF=N throughout
- `tests/data/synthetic/gatk_sample_A.vcf.gz.tbi` - tabix index
- `tests/data/synthetic/gatk_sample_B.vcf.gz` - 8 GATK variants with non-overlapping positions (11001-51001)
- `tests/data/synthetic/gatk_sample_B.vcf.gz.tbi` - tabix index
- `tests/data/synthetic/freebayes_sample_A.vcf.gz` - 12 Freebayes variants: lowQUAL, QUALperAO, SAF=0, SAR=0, readPosBias, DPu10het, DPu5hom, VAFu02het, VAFo08het, VAFu095hom triggers
- `tests/data/synthetic/freebayes_sample_A.vcf.gz.tbi` - tabix index
- `tests/data/synthetic/freebayes_sample_B.vcf.gz` - 7 Freebayes variants with non-overlapping positions (11001-61001)
- `tests/data/synthetic/freebayes_sample_B.vcf.gz.tbi` - tabix index

## Decisions Made

**REF=N for all variants:** The mini_ref.fa (chr22:1-100001 from hs37d5) is entirely composed of N bases (telomeric region). Using `REF=N` in all VCF records avoids `REF_MISMATCH` warnings from `bcftools norm --check-ref w`. bcftools correctly classifies `N>A` as SNP and `N>NA` as INDEL, so TYPE-based filters (gatkSNPhard, gatkINDELhard) work as expected.

**Freebayes INFO field placement:** AO, SAF, SAR, RPR, RPL are declared as `##INFO` (site-level) not `##FORMAT`. The Freebayes filter expressions use `INFO/SAF==0` syntax. Declaring them as FORMAT would break all strand bias and read position bias filters.

**FORMAT/VAF not in header:** The plan specifies that `bcftools +fill-tags` computes VAF at test time from `FORMAT/AD`. Declaring VAF in the header was explicitly avoided to test the fill-tags pipeline.

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 1 - Bug] Used REF=N instead of REF=A throughout all VCFs**

- **Found during:** Task 1 verification
- **Issue:** Plan said to "verify REF bases against mini_ref.fa" and use the actual base. Ran verification and found all positions in mini_ref.fa are N (telomeric region). Initial VCFs used REF=A which caused 16 REF_MISMATCH warnings per file, failing verification criterion 6.
- **Fix:** Rewrote all VCF records with REF=N. For SNPs: REF=N, ALT=A. For insertions: REF=N, ALT=NA. Confirmed bcftools TYPE() correctly classifies these as SNP/INDEL.
- **Files modified:** gatk_sample_A.vcf, gatk_sample_B.vcf, then freebayes_sample_A.vcf, freebayes_sample_B.vcf (all before bgzip compression)
- **Verification:** `bcftools norm --check-ref w ... 2>&1 | grep -c 'REF_MISMATCH'` returns 0 for all 4 files
- **Committed in:** eba15ec (Task 1), 63af0f8 (Task 2)

---

**Total deviations:** 1 auto-fixed (Rule 1 - Bug: REF base mismatch with all-N reference)
**Impact on plan:** Required fix - without it, bcftools norm --check-ref verification would fail. REF=N is semantically correct and functionally equivalent for filter testing purposes.

## Issues Encountered

**Verification command V5 ambiguity:** The plan's verification command `grep -c '\.'` on VAF output values is ambiguous - `.` in regex matches any character, so it would match all numeric values (e.g., `0.666667`). The correct command to detect missing VAF markers is `grep -c '^\.$'`. All VCFs produce 0 missing VAF values when using the correct regex. The plan's intent (no missing VAF) is met; the verification expression was corrected in practice.

**bioinformatics tools on Windows:** Tools not available natively; accessed via WSL (`wsl bash -c '...'`) with `/home/bernt/miniforge3/envs/hardnormly/bin` in PATH. Same approach as Plan 02-01.

## User Setup Required

None - no external service configuration required.

## Next Phase Readiness

- All 4 synthetic VCFs committed and ready for filter testing
- Plan 02-03 can generate expected output VCFs by running hardnormly.sh with each sample and filter file
- Plan 02-04 can write BATS tests that run filters against these VCFs and compare to expected outputs
- No blockers

---
*Phase: 02-test-data*
*Completed: 2026-02-18*
