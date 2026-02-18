# Phase 5: Features & Docs - Context

**Gathered:** 2026-02-18
**Status:** Ready for planning

<domain>
## Phase Boundary

Add user-facing CLI capabilities to hardnormly: subcommand dispatcher, convenience flags (--caller, --strip-annotations), BED generation utilities, non-fatal plot errors, and compact help text. This phase builds on the modular lib/ structure from Phase 4.

</domain>

<decisions>
## Implementation Decisions

### Subcommand UX
- No-args behavior: show full --help output (like git without args)
- Each subcommand (run-pipeline, generate-inclusion-bed, generate-exclusion-bed) has its own --help
- generate-*-bed subcommands are silent by default; use -v/--verbose for progress
- Subcommands are the primary/promoted interface going forward
- Old style (`hardnormly.sh -v input.vcf ...`) still works but isn't featured in docs — backward compat only

### --caller flag
- Accepts `gatk` and `freebayes` only; unknown values produce an error
- If both --caller and --filters-file are provided, --filters-file wins (with a warning)
- --caller and --filters-file shown at equal prominence in help text
- --caller maps to defaults/gatk_filters.txt or defaults/freebayes_filters.txt

### --strip-annotations flag
- Comma-separated syntax: `--strip-annotations INFO/CSQ,INFO/ANN`
- Matches bcftools annotate -x syntax
- Applied in the pipeline before filtering step

### Plot-vcfstats error handling
- Plot-vcfstats failures are logged but non-fatal — pipeline continues and exits 0

### Exclusion BED generation
- Helper script downloads and merges from 4 public sources: ENCODE blacklist, segmental duplications (UCSC Superdups), low-complexity regions (RepeatMasker/SDUST), centromeres/telomeres (UCSC gaps)
- Supports both hg19 and hg38 genome builds
- Download script model: user runs when needed, no pre-built BEDs in repo
- Produces one merged exclusion BED per build (hg19_exclusion.bed, hg38_exclusion.bed)

### Help text
- Compact --help: flag list + one-liner descriptions, fits one terminal screen
- One usage example: minimal invocation showing only required args
- Filter file format documented in README only (not in --help) — deviates from original DOCS-02 wording
- Exclusion BED source documentation (URLs, builds, regions) in README section, not in script --help

### Claude's Discretion
- Subcommand dispatcher implementation pattern
- Per-subcommand --help layout and wording
- Warning message wording for --caller + --filters-file conflict
- Error message for unknown --caller value
- Helper script internal structure and download mechanism

</decisions>

<specifics>
## Specific Ideas

- Subcommand style should feel modern — subcommands primary, flags-only as legacy compat
- --strip-annotations follows bcftools annotate -x comma-separated convention for familiarity
- Exclusion BED helper is a standalone download-and-merge script, not integrated into the main pipeline

</specifics>

<deferred>
## Deferred Ideas

None — discussion stayed within phase scope

</deferred>

---

*Phase: 05-features-docs*
*Context gathered: 2026-02-18*
