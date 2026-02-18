# Milestones: hardnormly

## v0.7.0 Code Quality, Testing & Features (Shipped: 2026-02-18)

**Delivered:** Transformed hardnormly from a 463-line monolithic bash script with no tests into a well-structured, linted, tested, modular tool with new subcommand capabilities.

**Phases completed:** 1-5 (19 plans total)

**Key accomplishments:**

- Linting infrastructure with ShellCheck + shfmt, pre-commit hooks, and GitHub Actions CI enforcing quality automatically on every commit and push
- 130-test BATS suite covering smoke, GATK/Freebayes filter units, integration on real 1000G data, regression against golden files, and lib/ module unit tests
- Synthetic + real test data: VCFs covering all filter triggers for both callers, plus 1000 Genomes NA12878/NA19247/HG00096 real data subsets
- Modularized into 8 lib/ modules (logging, cli, genome, bed, annotate, normalize, filter, stats); main orchestrator reduced from 463 to ~219 lines
- Array-based filter pipeline replacing `eval` on dynamically built strings with safe varargs construction
- New subcommands (`generate-inclusion-bed`, `generate-exclusion-bed`, `run-pipeline`), `--caller` and `--strip-annotations` flags, enhanced help text

**Stats:**

- 169 files created/modified
- ~14,043 lines of bash + bats
- 5 phases, 19 plans, 44 requirements (100% shipped)
- 1 day (2026-02-18 → 2026-02-18) — intensive single-day delivery
- 92 commits, 22,666 insertions

**Git range:** `507cef9` → `cc56540`

**Archive:** `.planning/milestones/v0.7.0-ROADMAP.md`

**What's next:** Next milestone to be determined via `/gsd:new-milestone`

---

## v0.6.0 (Pre-GSD — Existing)

**Status:** Shipped (before GSD tracking)

**What shipped:**
- Core VCF normalization and hard filtering pipeline
- Region-based annotation (INCLUDE_REGION/EXCLUDE_REGION)
- Caller-specific filter files (GATK HaplotypeCaller, Freebayes)
- BED file normalization and merging
- Genome file creation from UCSC MySQL
- Stats generation and plotting (plot-vcfstats)
- Compressed/indexed output support
- Snakemake workflow for batch processing
- Config validation, profile-based resources
- Cluster profiles (default, charite SLURM)

**Last phase:** 0 (pre-GSD)

---
