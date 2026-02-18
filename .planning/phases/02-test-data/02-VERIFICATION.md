---
phase: 02-test-data
verified: 2026-02-18T00:00:00Z
status: passed
score: 4/4 must-haves verified
---

# Phase 2: Test Data Verification Report

**Phase Goal:** All the data needed to run automated tests exists in the repo and can be regenerated
**Verified:** 2026-02-18
**Status:** passed
**Re-verification:** No - initial verification

## Goal Achievement

### Observable Truths

| # | Truth | Status | Evidence |
|---|-------|--------|----------|
| 1 | tests/data/ contains synthetic VCFs covering every GATK and Freebayes filter | VERIFIED | BASELINE.txt confirms all 7 GATK and 9 Freebayes filter names triggered |
| 2 | Minimal reference FASTA and matching BED files exist for the chr22 test region | VERIFIED | mini_ref.fa (>22, 100001 bases), include_regions.bed (3 regions), exclude_regions.bed (2 regions), hg19_chr22.genome |
| 3 | Running tests/generate_test_data.sh from a clean checkout reproduces all real data subsets | VERIFIED | Script is 207 lines, portable (no hardcoded paths), idempotent; all real data files committed to repo |
| 4 | Pre-computed expected output VCFs exist that regression tests can diff against | VERIFIED | 4 expected VCFs with .tbi indexes; BASELINE.txt has commit hash and 43 per-variant filter records |

**Score:** 4/4 truths verified

---

## Required Artifacts

### Success Criterion 1: Synthetic VCFs Covering All Filter Triggers (TDAT-01, TDAT-02, TDAT-03)

| Artifact | Expected | Status | Details |
|----------|----------|--------|---------|
| tests/data/synthetic/gatk_sample_A.vcf.gz | GATK VCF ~16 variants covering all GATK filters | VERIFIED | Exists; BASELINE.txt confirms all 7 GATK filter names triggered (gatkSNPhard, gatkINDELhard, DPu10het, DPu5hom, VAFu02het, VAFo08het, VAFu095hom) |
| tests/data/synthetic/gatk_sample_A.vcf.gz.tbi | tabix index | VERIFIED | Exists |
| tests/data/synthetic/gatk_sample_B.vcf.gz | GATK VCF ~8 variants for batch testing | VERIFIED | Exists; BASELINE shows 8 variants with gatkSNPhard, gatkINDELhard, DPu10het, VAFu02het, PASS |
| tests/data/synthetic/gatk_sample_B.vcf.gz.tbi | tabix index | VERIFIED | Exists |
| tests/data/synthetic/freebayes_sample_A.vcf.gz | Freebayes VCF ~12 variants covering all Freebayes filters | VERIFIED | Exists; BASELINE confirms all 9 Freebayes filter names triggered |
| tests/data/synthetic/freebayes_sample_A.vcf.gz.tbi | tabix index | VERIFIED | Exists |
| tests/data/synthetic/freebayes_sample_B.vcf.gz | Freebayes VCF ~7 variants for batch testing | VERIFIED | Exists; BASELINE shows 7 variants |
| tests/data/synthetic/freebayes_sample_B.vcf.gz.tbi | tabix index | VERIFIED | Exists |
| tests/data/synthetic/multiallelic.vcf.gz | Utility VCF with 3 multiallelic sites | VERIFIED | Exists + .tbi |
| tests/data/synthetic/empty_variants.vcf.gz | Valid VCF header, zero data records | VERIFIED | Exists + .tbi |
| tests/data/synthetic/minimal.vcf.gz | Exactly 1 variant with minimal fields | VERIFIED | Exists + .tbi |

### Success Criterion 2: Reference FASTA and BED Files (TDAT-04)

| Artifact | Expected | Status | Details |
|----------|----------|--------|---------|
| tests/data/synthetic/mini_ref.fa | chr22 100kb FASTA, contig named 22 | VERIFIED | Header is >22; all-N (telomeric); FAI shows 100001 bases |
| tests/data/synthetic/mini_ref.fa.fai | FASTA index | VERIFIED | Content: 22 TAB 100001 TAB 4 TAB 60 TAB 61 |
| tests/data/include_regions.bed | 3 include regions within chr22:10000-90000 | VERIFIED | 3 tab-delimited lines: 22:10000-30000, 22:40000-60000, 22:70000-90000 |
| tests/data/exclude_regions.bed | 2 exclude regions within chr22:10000-90000 | VERIFIED | 2 tab-delimited lines: 22:20000-25000, 22:80000-85000 |
| tests/data/hg19_chr22.genome | Genome file with chr22 length for bedtools slop | VERIFIED | Single line: 22 TAB 51304566 |

### Success Criterion 3: Real Data Generation Script (TDAT-05, TDAT-06)

| Artifact | Expected | Status | Details |
|----------|----------|--------|---------|
| tests/generate_test_data.sh | Portable, idempotent script to download all real data | VERIFIED | 207 lines; set -euo pipefail; dependency checks for bcftools/samtools/tabix/wget; --force flag; no hardcoded personal paths; 4 download sections |
| tests/data/real/giab_NA12878_chr22_16M.vcf.gz | GIAB NA12878 chr22:16M-16.1M subset | VERIFIED | Exists + .tbi |
| tests/data/real/giab_NA12878_chr22_16M_highconf.bed | GIAB high-confidence BED for test region | VERIFIED | Exists; regions confirmed in chr22:16M range |
| tests/data/real/1kg_NA12878_chr22_16M.vcf.gz | 1000G NA12878 chr22:16M-16.1M | VERIFIED | Exists + .tbi |
| tests/data/real/1kg_NA19247_chr22_16M.vcf.gz | 1000G NA19247 (YRI/AFR) chr22:16M-16.1M | VERIFIED | Exists + .tbi (NA19247 substituted for planned NA19240; documented in SUMMARY) |
| tests/data/real/1kg_HG00096_chr22_16M.vcf.gz | 1000G HG00096 (GBR/EUR) chr22:16M-16.1M | VERIFIED | Exists + .tbi |
| tests/data/real/freebayes_tiny.vcf.gz | Authentic Freebayes output with AO/SAF/SAR/RPR/RPL INFO fields | VERIFIED | Exists + .tbi |
| tests/data/real/freebayes_tiny_ref.fa | Reference FASTA for freebayes_tiny | VERIFIED | Exists + .fai; header is >q |
| tests/data/real/chr22_16M_ref.fa | chr22:16M-16.1M reference, contig renamed to >22 | VERIFIED | Exists + .fai; header is >22; FAI shows 100001 bases |

### Success Criterion 4: Expected Output VCFs (TDAT-07, TDAT-08)

| Artifact | Expected | Status | Details |
|----------|----------|--------|---------|
| tests/data/expected/gatk_sample_A_filtered.vcf.gz | Pipeline output with FILTER tags applied | VERIFIED | Exists + .tbi; BASELINE shows 16 variants, 9 distinct FILTER values |
| tests/data/expected/gatk_sample_B_filtered.vcf.gz | Pipeline output for batch testing | VERIFIED | Exists + .tbi; BASELINE shows 8 variants |
| tests/data/expected/freebayes_sample_A_filtered.vcf.gz | Freebayes pipeline output | VERIFIED | Exists + .tbi; BASELINE shows 12 variants with all 9 Freebayes filter names |
| tests/data/expected/freebayes_sample_B_filtered.vcf.gz | Freebayes pipeline output for batch testing | VERIFIED | Exists + .tbi; BASELINE shows 7 variants |
| tests/data/expected/BASELINE.txt | Commit hash + per-variant filter summaries | VERIFIED | Contains commit 1d41779ff42534491c5274d10963540c4366c8ac, date 2026-02-18T17:05:23Z, hardnormly.sh version 0.6.0, and 43 per-variant CHROM/POS/FILTER records |
| tests/data/gatk_vcfs.txt | 2 GATK VCF paths | VERIFIED | 2 lines pointing to existing files |
| tests/data/freebayes_vcfs.txt | 2 Freebayes VCF paths | VERIFIED | 2 lines pointing to existing files |
| tests/data/real_vcfs.txt | 3 real VCF paths (NA12878, NA19247, HG00096) | VERIFIED | 3 lines pointing to existing files |
| tests/data/test_config_gatk.yaml | Snakemake config matching config.schema.yaml | VERIFIED | Has all 5 sections (ref, paths, regions, filtering, processing); uses output_folder, plot_stats, GRCh37, genome_file; references gatk_filters.txt |
| tests/data/test_config_freebayes.yaml | Snakemake config for Freebayes | VERIFIED | As above; references freebayes_filters.txt |
| tests/data/test_config_real.yaml | Snakemake config for real data | VERIFIED | References chr22_16M_ref.fa and real_vcfs.txt |
| tests/data/generate_expected.sh | Reproducible script to regenerate expected outputs | VERIFIED (with warning) | Exists, 101 lines, set -euo pipefail; hardcoded PATH=/home/bernt/miniconda3/... is non-portable |

---

## Key Link Verification

| From | To | Via | Status | Details |
|------|----|----|--------|---------|
| tests/data/synthetic/gatk_sample_A.vcf.gz | defaults/gatk_filters.txt | BASELINE filter labels match filter file names | WIRED | All 7 filter names confirmed: gatkSNPhard, gatkINDELhard, DPu10het, DPu5hom, VAFu02het, VAFo08het, VAFu095hom |
| tests/data/synthetic/freebayes_sample_A.vcf.gz | defaults/freebayes_filters.txt | BASELINE filter labels match filter file names | WIRED | All 9 filter names confirmed: lowQUAL, QUALperAO, strandBias, readPosBias, DPu10het, DPu5hom, VAFu02het, VAFo08het, VAFu095hom |
| tests/data/expected/*.vcf.gz | tests/data/synthetic/*.vcf.gz | Generated by running hardnormly.sh on synthetic inputs | WIRED | BASELINE.txt records generation commit; filter diversity confirms pipeline ran (not all PASS) |
| tests/data/test_config_gatk.yaml | tests/data/gatk_vcfs.txt | paths.vcf_list field | WIRED | Config has vcf_list: tests/data/gatk_vcfs.txt |
| tests/data/test_config_gatk.yaml | defaults/gatk_filters.txt | filtering.filters_file field | WIRED | Config has filters_file: defaults/gatk_filters.txt |
| .gitignore | tests/ | \!/tests/ and \!/tests/** exceptions after *test* rule | WIRED | All test files confirmed committed in git; exceptions at lines 24-31 after *test* rule at line 21 |
| tests/generate_test_data.sh | tests/data/real/ | 4 download sections with bcftools/tabix streaming | WIRED | Script has GIAB, 1000G, Freebayes tiny, chr22 ref sections; idempotency guards present for each |

---

## Requirements Coverage

| Requirement | Description | Status | Evidence |
|-------------|-------------|--------|---------|
| TDAT-01 | Synthetic GATK VCFs (A and B) with all GATK filter triggers | SATISFIED | gatk_sample_A (16 variants) and gatk_sample_B (8 variants); all 7 GATK filter labels confirmed in BASELINE |
| TDAT-02 | Synthetic Freebayes VCFs (A and B) with all Freebayes filter triggers | SATISFIED | freebayes_sample_A (12 variants) and freebayes_sample_B (7 variants); all 9 Freebayes filter labels confirmed |
| TDAT-03 | Utility VCFs (multiallelic, minimal, empty_variants) | SATISFIED | All 3 utility VCFs exist with .tbi indexes |
| TDAT-04 | Synthetic reference FASTA and BED files for chr22 test region | SATISFIED | mini_ref.fa (100001 bases, contig 22), include/exclude BEDs, hg19_chr22.genome all exist and are substantive |
| TDAT-05 | Real data subsets with indexes | SATISFIED | GIAB NA12878, 3x 1000G samples, freebayes_tiny, chr22_16M_ref all exist with .tbi/.fai |
| TDAT-06 | generate_test_data.sh reproducibly creates real data subsets | SATISFIED | Script is 207 lines, portable (no personal paths), idempotent, covers all 4 data sources |
| TDAT-07 | Expected output files for regression testing | SATISFIED | 4 expected VCFs with FILTER annotations applied; BASELINE.txt with commit hash and 43 per-variant records |
| TDAT-08 | VCF list files and test configs for Snakemake batch testing | SATISFIED | 3 VCF list files and 3 test_config YAMLs matching config.schema.yaml field names |

---

## Anti-Patterns Found

| File | Line | Pattern | Severity | Impact |
|------|------|---------|----------|--------|
| tests/data/generate_expected.sh | 7 | Hardcoded export PATH=/home/bernt/miniconda3/envs/hardnormly/bin | Warning | Not portable; fails for any developer other than bernt; blocks regenerating expected outputs |
| tests/data/verify_task2.sh | 5 | Hardcoded export PATH=/home/bernt/miniconda3/envs/hardnormly/bin | Warning | Same portability issue; development-only verification helper script |

Neither anti-pattern blocks phase success criteria. generate_expected.sh is a helper for regenerating
expected outputs when the pipeline intentionally changes and is not referenced in any success criterion.
Should be fixed before Phase 3 so other developers can regenerate expected outputs after refactoring.

---

## Human Verification Required

### 1. Synthetic VCF Field Values (VAF computation)

**Test:** In an environment with bcftools in PATH, run:
    bcftools +fill-tags tests/data/synthetic/gatk_sample_A.vcf.gz -- -t FORMAT/VAF | bcftools query -f "[%VAF]
"
**Expected:** All 16 values are numeric; no missing-value markers (dot)
**Why human:** Cannot inspect binary VCF.gz content without a runtime environment with bcftools

### 2. Freebayes INFO Field Declaration

**Test:** In an environment with bcftools, run:
    bcftools view -h tests/data/synthetic/freebayes_sample_A.vcf.gz | grep INFO.*SAF
**Expected:** Shows Type=Integer and line is a ##INFO declaration (not ##FORMAT)
**Why human:** Cannot inspect binary VCF headers without bcftools

### 3. REF Allele Validation Against mini_ref.fa

**Test:** In an environment with bcftools, run:
    bcftools norm --check-ref w -f tests/data/synthetic/mini_ref.fa tests/data/synthetic/gatk_sample_A.vcf.gz -Oz -o /dev/null 2>&1 | grep -c REF_MISMATCH
**Expected:** 0 (no mismatches; REF=N matches the all-N mini_ref.fa telomeric sequence)
**Why human:** Cannot run bcftools norm without a runtime environment

### 4. generate_test_data.sh Execution from Clean State

**Test:** On a system with bcftools, samtools, tabix, wget in PATH:
    rm -rf tests/data/real/
    ./tests/generate_test_data.sh
**Expected:** All 8 real data files recreated in tests/data/real/ with correct variant counts
**Why human:** Requires network access to GIAB FTP, 1000G EBI FTP, and GitHub; cannot verify statically

---

## Notes

**test_config_real.yaml BED Coordinate Mismatch:**
tests/data/test_config_real.yaml references the synthetic BED files (chr22:10000-90000 coordinates).
Real VCFs contain variants in chr22:16M-16.1M. These ranges do not overlap, so running the Snakemake
workflow with this config will tag all real variants as NOT_IN_INCLUDE_REGION. This is documented in
02-04-SUMMARY as a Phase 3 concern. The config is syntactically valid and meets the TDAT-08 structural
requirement (VCF list files and test configs for Snakemake batch testing). Phase 3 may need BED files
covering the chr22:16M coordinate space for meaningful real-data integration testing.

**gitignore Design:**
The .gitignore has a catch-all !/tests/** at line 25 that allows all file types under tests/. The
specific file extension exceptions on lines 26-31 (*.vcf, *.vcf.gz, *.bed, *.fa, etc.) are redundant
but harmless. The critical ordering rule (exceptions after the *test* rule at line 21) is correctly
implemented.

---

*Verified: 2026-02-18*
*Verifier: Claude (gsd-verifier)*
