---
phase: 02-test-data
plan: 03
subsystem: testing
tags: [vcf, bgzip, tabix, bcftools, samtools, giab, 1000genomes, freebayes, test-data]

# Dependency graph
requires:
  - phase: 02-01
    provides: mini_ref.fa (chr22 100kb reference), tests/data/ directory structure, .gitignore exceptions
provides:
  - tests/data/synthetic/multiallelic.vcf.gz - 3 multiallelic sites that split into 6 biallelic records via bcftools norm -m-any
  - tests/data/synthetic/empty_variants.vcf.gz - valid VCF header with zero data records
  - tests/data/synthetic/minimal.vcf.gz - exactly 1 variant with GT+DP+AD fields only
  - tests/generate_test_data.sh - reproducible download script for all real data subsets (206 lines, idempotent, --force support)
  - tests/data/real/giab_NA12878_chr22_16M.vcf.gz - GIAB HG001 v4.2.1 benchmark, 7 variants in chr22:16M-16.1M
  - tests/data/real/giab_NA12878_chr22_16M_highconf.bed - GIAB high-confidence regions for chr22:16M-16.1M
  - tests/data/real/1kg_NA12878_chr22_16M.vcf.gz - 1000G Phase 3 NA12878, 20 variants
  - tests/data/real/1kg_NA19247_chr22_16M.vcf.gz - 1000G Phase 3 NA19247 (YRI), 30 variants
  - tests/data/real/1kg_HG00096_chr22_16M.vcf.gz - 1000G Phase 3 HG00096 (GBR), 17 variants
  - tests/data/real/freebayes_tiny.vcf.gz - authentic Freebayes output with AO/SAF/SAR/RPR/RPL INFO fields
  - tests/data/real/freebayes_tiny_ref.fa - matching reference FASTA for freebayes_tiny.vcf.gz
  - tests/data/real/chr22_16M_ref.fa - chr22:16M-16.1M from hs37d5, contig renamed to >22
affects:
  - 02-04 (BATS tests need all VCFs; multiallelic, empty, minimal test edge cases)
  - 03-xx (integration tests need real data subsets and chr22_16M_ref.fa)
  - all future test phases

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "WSL bioinformatics tool invocation: wsl -d Ubuntu -u bernt -- sh -c with MSYS_NO_PATHCONV considerations"
    - "Idempotent download scripts with [[ ! -f ]] existence checks and --force flag"
    - "Remote 1000G tabix streaming: tabix -h URL REGION | bcftools view -s SAMPLE --min-ac 1"
    - "VCF creation workflow: write plain .vcf -> bgzip -f -> bcftools index -t"
    - "MSYS path rewriting: avoid shell variable expansion via env -i and absolute paths in wsl invocations"

key-files:
  created:
    - tests/data/synthetic/multiallelic.vcf.gz
    - tests/data/synthetic/multiallelic.vcf.gz.tbi
    - tests/data/synthetic/empty_variants.vcf.gz
    - tests/data/synthetic/empty_variants.vcf.gz.tbi
    - tests/data/synthetic/minimal.vcf.gz
    - tests/data/synthetic/minimal.vcf.gz.tbi
    - tests/generate_test_data.sh
    - tests/data/real/giab_NA12878_chr22_16M.vcf.gz
    - tests/data/real/giab_NA12878_chr22_16M.vcf.gz.tbi
    - tests/data/real/giab_NA12878_chr22_16M_highconf.bed
    - tests/data/real/1kg_NA12878_chr22_16M.vcf.gz
    - tests/data/real/1kg_NA12878_chr22_16M.vcf.gz.tbi
    - tests/data/real/1kg_NA19247_chr22_16M.vcf.gz
    - tests/data/real/1kg_NA19247_chr22_16M.vcf.gz.tbi
    - tests/data/real/1kg_HG00096_chr22_16M.vcf.gz
    - tests/data/real/1kg_HG00096_chr22_16M.vcf.gz.tbi
    - tests/data/real/freebayes_tiny.vcf.gz
    - tests/data/real/freebayes_tiny.vcf.gz.tbi
    - tests/data/real/freebayes_tiny_ref.fa
    - tests/data/real/freebayes_tiny_ref.fa.fai
    - tests/data/real/chr22_16M_ref.fa
    - tests/data/real/chr22_16M_ref.fa.fai
  modified:
    - tests/generate_test_data.sh (NA19240 -> NA19247 fix during Task 3)

key-decisions:
  - "REF bases in synthetic VCFs are N because mini_ref.fa (chr22:1-100001 from hs37d5) is all-N telomeric sequence"
  - "NA19247 substituted for NA19240 — NA19240 is not present in 1000G Phase 3 chr22 dataset; NA19247 is also YRI/AFR"
  - "chr22_16M_ref.fa uses sed to rename contig from >22:16000000-16100000 to >22 for bcftools norm compatibility"
  - "Remote tabix streaming used for 1000G to avoid downloading 196MB chr22 file"
  - "WSL invocations use env -i with explicit LD_LIBRARY_PATH to avoid MSYS path rewriting issues"

patterns-established:
  - "WSL bioinformatics tool pattern: env -i HOME=/home/bernt LD_LIBRARY_PATH=$BCENV/lib $BCENV/bin/tool args"
  - "Idempotent download: [[ FORCE == true ]] || [[ ! -f OUTPUT ]] guards every download block"
  - "Multiallelic VCF splitting test: bcftools norm -m-any reduces 3 sites to 6 biallelic records"

# Metrics
duration: ~45min
completed: 2026-02-18
---

# Phase 2 Plan 3: Utility VCFs and Real Data Subsets Summary

**3 utility VCFs (multiallelic/empty/minimal) + reproducible generate_test_data.sh script + 14 real data files from GIAB, 1000G, and Freebayes — all committed and indexed**

## Performance

- **Duration:** ~45 min
- **Started:** 2026-02-18T17:20:00Z
- **Completed:** 2026-02-18T18:00:00Z
- **Tasks:** 3
- **Files modified:** 22 (3 utility VCFs + their indexes, generate_test_data.sh, 8 real VCFs + indexes, 1 BED, 2 FASTAs + indexes)

## Accomplishments

- Created 3 utility VCFs for edge case testing: multiallelic (3 sites -> 6 biallelic via norm -m-any), empty_variants (valid header, 0 records), minimal (1 variant with GT/DP/AD only)
- Wrote 206-line generate_test_data.sh script covering all 4 data sources with idempotency, dependency checks, and --force support
- Downloaded and committed all real data subsets: GIAB NA12878 (7 variants), 3 x 1000G samples (17-30 variants each), Freebayes tiny (16 records with authentic INFO fields), chr22 reference FASTA

## Task Commits

Each task was committed atomically:

1. **Task 1: Create utility VCFs** - `0ecece8` (feat)
2. **Task 2: Write generate_test_data.sh** - `3b59382` (feat)
3. **Task 3: Run script and download real data** - `194b056` (feat)

**Plan metadata:** (this commit)

## Files Created/Modified

- `tests/data/synthetic/multiallelic.vcf.gz` - 3 multiallelic sites, splits to 6 biallelic via bcftools norm -m-any; REF=N throughout (mini_ref.fa is all-N telomeric sequence)
- `tests/data/synthetic/empty_variants.vcf.gz` - valid VCFv4.2 header with zero data records
- `tests/data/synthetic/minimal.vcf.gz` - exactly 1 variant (22:10001 N->G, het 0/1, DP=20, AD=10,10)
- `tests/generate_test_data.sh` - reproducible download script, 206 lines, idempotent
- `tests/data/real/giab_NA12878_chr22_16M.vcf.gz` - GIAB HG001 v4.2.1 benchmark, 7 variants
- `tests/data/real/giab_NA12878_chr22_16M_highconf.bed` - high-confidence regions
- `tests/data/real/1kg_NA12878_chr22_16M.vcf.gz` - 20 variants
- `tests/data/real/1kg_NA19247_chr22_16M.vcf.gz` - 30 variants (YRI/AFR sample)
- `tests/data/real/1kg_HG00096_chr22_16M.vcf.gz` - 17 variants (GBR/EUR sample)
- `tests/data/real/freebayes_tiny.vcf.gz` - 16 records, AO/SAF/SAR/RPR/RPL INFO fields confirmed
- `tests/data/real/freebayes_tiny_ref.fa` - Freebayes test tiny reference
- `tests/data/real/chr22_16M_ref.fa` - chr22:16M-16.1M, contig header fixed to >22

## Decisions Made

- REF bases in utility VCFs are N because mini_ref.fa (chr22:1-100001 from hs37d5) is all-N telomeric sequence; bcftools accepts N as REF without complaint and norm -m-any splits correctly
- NA19247 used instead of planned NA19240 (NA19240 not in 1000G Phase 3 chr22 dataset; NA19247 is same population YRI/AFR)
- chr22_16M_ref.fa contig renamed via `sed '1s/.*/\>22/'` to match VCF CHROM=22 for bcftools norm compatibility

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 1 - Bug] NA19240 substituted with NA19247**
- **Found during:** Task 3 (run generate_test_data.sh)
- **Issue:** Plan specified 1000G sample NA19240 (labeled AFR/YRI), but NA19240 is not present in the 1000G Phase 3 chr22 dataset. `bcftools view -s NA19240` returned "subset called for sample that does not exist in header"
- **Fix:** Verified available NA19xxx samples; NA19247 (also YRI/AFR) is present and has 30 variants in the region. Updated generate_test_data.sh and file naming throughout
- **Files modified:** tests/generate_test_data.sh (NA19240 -> NA19247), output file is 1kg_NA19247_chr22_16M.vcf.gz
- **Verification:** 1kg_NA19247_chr22_16M.vcf.gz created with 30 variants; bcftools query confirms valid content
- **Committed in:** `194b056` (Task 3 commit)

**2. [Rule 3 - Blocking] MSYS path rewriting prevents shell variable use in WSL invocations**
- **Found during:** Task 1 (compress VCFs)
- **Issue:** Git Bash on Windows rewrites /bin/, /tmp/, /usr/ paths inside shell arguments even when quoted. Variables like `BCENV=/home/bernt/...` are corrupted to empty strings because `/home` is rewritten. This broke all compound WSL commands using shell variables.
- **Fix:** Used `env -i HOME=... LD_LIBRARY_PATH=... /full/path/to/binary args` pattern consistently. All tool invocations use fully explicit paths with no variable expansion in critical segments.
- **Files modified:** No files; execution approach adapted
- **Verification:** All bcftools/bgzip/tabix operations completed successfully with env -i pattern
- **Committed in:** N/A (execution approach, no file change)

---

**Total deviations:** 2 (1 auto-fixed bug, 1 blocking issue resolved)
**Impact on plan:** NA19240->NA19247 substitution is equivalent coverage (both YRI/AFR, same dataset); MSYS workaround is execution-only and doesn't affect output files.

## Issues Encountered

- Initial bgzip compression of multiallelic.vcf prompted Vim pager in background tasks; solved by using explicit output redirection (`-Oz -o /dev/null`) and `bcftools query` (no pager) for verification
- First generate_test_data.sh run timed out at 120s during GIAB download; subsequent run with 600s timeout completed successfully; GIAB + NA12878 1000G files already partially downloaded from first attempt

## User Setup Required

None - no external service configuration required. All data is downloaded automatically by generate_test_data.sh.

## Next Phase Readiness

- All test data in place for Phase 2 Plan 4 (BATS test framework)
- Utility VCFs tested: multiallelic splits correctly, empty handles gracefully, minimal is minimal
- Real data covers diverse variant callers (GATK via GIAB, 1000G, Freebayes) and populations (EUR, AFR)
- generate_test_data.sh is idempotent and reproducible from a clean checkout
- Note for Plan 04: `1kg_NA19240` is now `1kg_NA19247` — update any BATS test references accordingly

---
*Phase: 02-test-data*
*Completed: 2026-02-18*
