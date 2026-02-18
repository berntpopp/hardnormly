# Phase 5: Features & Docs - Research

**Researched:** 2026-02-18
**Domain:** Bash CLI subcommand dispatcher, bcftools annotate -x, exclusion BED download scripts, non-fatal error handling in set -Eeuo pipefail
**Confidence:** HIGH (verified against actual codebase and bcftools 1.20 installation; patterns cross-referenced with established project conventions)

## Summary

Phase 5 adds user-facing CLI capabilities to hardnormly.sh on top of the modular lib/ structure from Phase 4. The work divides into five independent technical areas: (1) subcommand dispatcher in hardnormly.sh and cli.sh, (2) --caller flag mapping to default filter files, (3) --strip-annotations flag via bcftools annotate -x, (4) non-fatal plot-vcfstats error handling, and (5) a standalone exclusion BED download helper script.

The standard approach for all five areas is verified: the subcommand dispatcher uses a case statement routing to shell functions, with the first argument either being a subcommand name or going straight to parse_args for backward compatibility. The --strip-annotations flag maps directly to bcftools annotate -x (confirmed syntax: comma-separated `INFO/CSQ,INFO/ANN`). Non-fatal plot-vcfstats errors are achieved by catching the return code before it propagates under set -e. The exclusion BED script downloads from four public UCSC/ENCODE sources via wget/curl and merges with bedtools.

The biggest implementation constraint is the SC2310/SC2311 pattern already established in this project: **do not use `||` with function calls under `set -Eeuo pipefail`**. Non-fatal handling requires the `local exit_code=0; fn ... || exit_code=$?` idiom or an `if fn; then ... fi` block.

**Primary recommendation:** Implement the subcommand dispatcher as a thin routing function in hardnormly.sh (not in cli.sh) that inspects `$1` before calling `parse_args`. Keep cli.sh focused on argument parsing. All new flags (--caller, --strip-annotations) are parsed in cli.sh and handled in hardnormly.sh.

---

## Standard Stack

No new library dependencies. All tooling was established in prior phases.

### Core Tools (unchanged)

| Tool | Version | Purpose | Notes |
|------|---------|---------|-------|
| bash | 4+ | Script runtime | Arrays, `[[ ]]`, namerefs all available |
| bcftools | 1.20 | VCF annotation, filtering | `annotate -x` confirmed syntax |
| bedtools | 2.31 | BED merging, sorting | Used in exclusion BED helper |
| bgzip + tabix | htslib 1.20 | BED compression/indexing | For merged output BEDs |
| wget/curl | system | Downloading exclusion source files | Either works; wget simpler for UCSC |
| shfmt | current | Code formatting | Flags: `-i 0 -bn -ci` (locked) |
| shellcheck | current | Static analysis | SC2310/SC2311 pattern must be followed |
| bats-core | 1.x | Testing | Already installed (Phase 3) |

### No New Dependencies

Phase 5 is purely about wiring new logic through existing tools. No new bash libraries, Python packages, or external tools are introduced.

**Installation:** No new packages needed.

---

## Architecture Patterns

### Recommended Project Structure (additions)

```
hardnormly/
├── hardnormly.sh                    # Gains: subcommand dispatcher, --caller, --strip-annotations
├── lib/
│   ├── cli.sh                       # Gains: --caller, --strip-annotations parsing; subcommand help functions
│   ├── stats.sh                     # plot_stats_output stays fatal internally; caller handles non-fatal
│   └── [unchanged modules]
└── scripts/
    └── generate_exclusion_bed.sh    # NEW: standalone download-and-merge helper for exclusion BEDs
```

### Pattern 1: Subcommand Dispatcher

**What:** Inspect `$1` in hardnormly.sh before passing to `parse_args`. If it matches a known subcommand, route to a subcommand handler function. If not, fall through to the legacy `parse_args "$@"` path.

**When to use:** Entry point only (hardnormly.sh). cli.sh handles parsing within each subcommand's scope.

**Key constraint:** The dispatcher must NOT be inside `parse_args`. It runs before `parse_args` so that backward-compatible invocations (`hardnormly.sh -v input.vcf ...`) still work by routing to `run-pipeline` implicitly.

```bash
# hardnormly.sh — dispatcher (runs BEFORE parse_args)
_dispatch() {
    local subcommand="${1:-}"
    case "$subcommand" in
        run-pipeline)
            shift
            # fall through to existing parse_args path
            parse_args "$@"
            ;;
        generate-inclusion-bed)
            shift
            cmd_generate_inclusion_bed "$@"
            exit 0
            ;;
        generate-exclusion-bed)
            shift
            cmd_generate_exclusion_bed "$@"
            exit 0
            ;;
        "" | -h | --help)
            show_help
            ;;
        -*)
            # Starts with a flag — legacy mode, route to run-pipeline
            parse_args "$@"
            ;;
        *)
            echo "Unknown subcommand: $subcommand" >&2
            show_help
            ;;
    esac
}
```

**Note on exit:** Subcommand handlers call `exit 0` after completing. `show_help` exits with 1 (existing behavior). The `run-pipeline` and legacy-flag paths do NOT exit here — execution continues to the pipeline.

**Note on no-args behavior:** `""` case routes to `show_help`, which exits 1 (matching existing smoke test: "no arguments exits non-zero"). This is consistent with both the user decision ("show full --help output like git without args") and the existing smoke test.

### Pattern 2: Per-Subcommand Help Functions

**What:** Each subcommand has its own help function (`show_help_generate_inclusion_bed`, `show_help_generate_exclusion_bed`) in cli.sh. The subcommand handler inspects its own args for `-h`/`--help` and calls its help function.

**When to use:** Inside each subcommand handler when `$1` is `--help` or `-h`.

```bash
# lib/cli.sh
show_help_generate_exclusion_bed() {
    echo "Usage: hardnormly.sh generate-exclusion-bed -o <output.bed> [options]"
    echo ""
    echo "Options:"
    echo "  -o, --output    Output BED file path (required)"
    echo "  --genome-build  Genome build: hg19 or hg38 (default: hg19)"
    echo "  -v, --verbose   Show download progress"
    echo "  -h, --help      Show this help"
    exit 0
}
```

**Note:** generate-*-bed subcommands use `exit 0` for their help (not `exit 1` like the main `show_help`). The decision was: subcommand help is informational, not an error.

### Pattern 3: --caller Flag Mapping

**What:** `--caller` is parsed in cli.sh's `parse_args`. A new `caller` variable is set. In hardnormly.sh, after `parse_args` and before `validate_args`, translate `caller` to `filters_file` if not already set.

**When to use:** After `parse_args` returns, before `validate_args` runs.

**Key constraint (from CONTEXT.md):** If both `--caller` and `--filters-file` are provided, `--filters-file` wins with a warning. This means translation happens only when `filters_file` is empty.

```bash
# hardnormly.sh — after parse_args, before validate_args
_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

if [[ -n "$caller" ]]; then
    case "$caller" in
        gatk)
            _caller_file="${_SCRIPT_DIR}/defaults/gatk_filters.txt"
            ;;
        freebayes)
            _caller_file="${_SCRIPT_DIR}/defaults/freebayes_filters.txt"
            ;;
        *)
            echo "Error: Unknown --caller value '$caller'. Valid values: gatk, freebayes" >&2
            exit 1
            ;;
    esac
    if [[ -n "$filters_file" ]]; then
        log_msg "Warning: Both --caller and --filters-file provided; --filters-file takes precedence."
    else
        filters_file="$_caller_file"
    fi
fi
```

**Note:** `caller` variable must be initialized to `""` in the global defaults section of hardnormly.sh (alongside `filters_file`). `_caller_file` is a local-scope variable here (not exported).

### Pattern 4: --strip-annotations Flag

**What:** `--strip-annotations INFO/CSQ,INFO/ANN` parses the comma-separated list in cli.sh. In hardnormly.sh, a `bcftools annotate -x "$strip_annotations"` step runs against the input VCF before the existing normalization step.

**Placement in pipeline:** Step 4.5 — after exclusion annotation (Step 4), before normalization (Step 5). This ensures INFO fields like CSQ/ANN are stripped from the working VCF before norm and filter steps.

**Verified syntax** (bcftools 1.20 confirmed):
```bash
# lib/annotate.sh or hardnormly.sh
strip_vcf_annotations() {
    local vcf_file="$1"
    local strip_list="$2"   # comma-separated: "INFO/CSQ,INFO/ANN"
    local output_vcf="$3"
    run_cmd bcftools annotate -x "$strip_list" "$vcf_file" -Oz -o "$output_vcf"
}
```

**Note on function placement:** `strip_vcf_annotations` belongs in `lib/annotate.sh` (alongside `annotate_vcf_with_regions`) since both are annotation-stage operations. This keeps filter.sh focused on the filter pipeline.

**Note on parse_args:** Add `strip_annotations=""` to the global defaults. Parse `--strip-annotations` just like `--filters-file` (takes a value argument). No validation needed in `validate_args` — an empty string means no stripping; bcftools is not called.

### Pattern 5: Non-Fatal Plot-vcfstats (FEAT-01)

**What:** `plot_stats_output` currently returns 1 on failure. Under `set -Eeuo pipefail`, this would exit the script. The fix is in hardnormly.sh: catch the return code before it propagates.

**Critical constraint (SC2310/SC2311 project pattern):** Do NOT use `plot_stats_output ... || true` or `plot_stats_output ... || log_msg "..."` — ShellCheck flags `||` on function calls as SC2310/SC2311. Instead, use the `if` block pattern:

```bash
# hardnormly.sh — non-fatal plot call (correct pattern for this project)
if [[ "$plot_stats" == "true" ]]; then
    log_msg "Plotting stats to $plot_output_dir"
    if ! plot_stats_output "$stats_output" "$plot_output_dir" "$tmp_dir"; then
        log_msg "Warning: plot-vcfstats failed; pipeline continues."
    fi
fi
```

**Why `if !` works but `||` does not:** `if ! fn` is a conditional — bash does not apply the `set -e` abort-on-failure rule inside `if` conditions. ShellCheck does not flag this as SC2310/SC2311. The `||` form triggers those warnings because it modifies function exit codes in a way that conflicts with errexit propagation detection.

**Note:** `plot_stats_output` itself does NOT change — it still returns 1 on failure. The non-fatal behavior is enforced at the call site in hardnormly.sh, not inside the function. This preserves testability of `plot_stats_output` (unit tests can still assert failure).

### Pattern 6: generate-inclusion-bed Subcommand (FEAT-03)

**What:** Wraps the existing `merge_include_beds` + `compress_index_bed` logic from bed.sh into a subcommand that accepts BED files as positional/flag args and writes a merged output.

**When to use:** User calls `hardnormly.sh generate-inclusion-bed -b file1.bed -b file2.bed -o merged.bed`.

**Key insight:** All the merging logic already exists in `lib/bed.sh`. This subcommand is just a CLI wrapper that calls `normalize_bed`, `merge_include_beds`, and optionally `compress_index_bed`. It requires a genome file for slop (same as the pipeline). The subcommand handler lives in hardnormly.sh (not cli.sh) since it orchestrates lib/ calls.

**Silent by default:** generate-*-bed subcommands log nothing by default. Use `-v/--verbose` to see progress. Implement by passing a verbosity flag through to logging — or simply suppress `log_msg` output by redirecting stdout to /dev/null unless verbose mode is on. The simpler approach: run the subcommand's own local arg parser that sets a `_verbose` variable, then only call `log_msg` if `_verbose=true`.

### Pattern 7: generate-exclusion-bed Subcommand (FEAT-04) and Helper Script (FEAT-05)

**What:** The subcommand calls a standalone helper script (`scripts/generate_exclusion_bed.sh`) that downloads from 4 public sources and merges into a single exclusion BED per genome build.

**Why a separate script:** The download logic is large (4 sources, awk transforms, wget calls) and user-run when needed — not part of the main pipeline. Keeping it in `scripts/` matches the existing `scripts/setup-hooks.sh` and `scripts/run_snakemake.sh` pattern.

**The 4 sources and their UCSC download URLs:**

| Source | hg19 URL | hg38 URL | Notes |
|--------|----------|----------|-------|
| ENCODE blacklist | `https://github.com/Boyle-Lab/Blacklist/raw/master/lists/hg19-blacklist.v2.bed.gz` | `https://github.com/Boyle-Lab/Blacklist/raw/master/lists/hg38-blacklist.v2.bed.gz` | Direct GitHub raw; confirmed in Boyle-Lab/Blacklist `lists/` folder |
| Segmental dups (UCSC SuperDups) | `https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/genomicSuperDups.txt.gz` | `https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/genomicSuperDups.txt.gz` | Tab-delimited; convert with `cut -f 2-4` to BED |
| Low complexity (RepeatMasker) | `https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/rmsk.txt.gz` | `https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/rmsk.txt.gz` | Filter for `Low_complexity` class; `cut -f 6-8` for BED |
| Centromeres/telomeres (UCSC gap) | `https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/gap.txt.gz` | hg38 gap table may be empty (T2T assembly) — use `cytoBandIdeo.txt.gz` for centromeres | `cut -f 2-4` for BED; filter type `centromere|telomere` from column 8 |

**MEDIUM confidence warning on hg38 gap/centromere URL:** hg38 moved to T2T centromere sequences; the gap table may have fewer/no centromere entries. Planner should note this as a verification step: download and inspect row count for hg38 gap.

**Download and transform pattern for helper script:**

```bash
# Example: SuperDups → BED
wget -q -O - "https://hgdownload.soe.ucsc.edu/goldenPath/${build}/database/genomicSuperDups.txt.gz" \
    | gunzip -c \
    | awk '{OFS="\t"; print $2, $3, $4}' \
    > "${tmp_dir}/superdups.bed"

# Example: RepeatMasker Low_complexity → BED
wget -q -O - "https://hgdownload.soe.ucsc.edu/goldenPath/${build}/database/rmsk.txt.gz" \
    | gunzip -c \
    | awk '$12 == "Low_complexity" {OFS="\t"; print $6, $7, $8}' \
    > "${tmp_dir}/lcr.bed"
```

**Final merge:**
```bash
bedtools multiinter -i "${tmp_dir}/blacklist.bed" "${tmp_dir}/superdups.bed" \
    "${tmp_dir}/lcr.bed" "${tmp_dir}/gaps.bed" \
    | bedtools sort -i - \
    | awk '{OFS="\t"; print $1, $2, $3}' \
    > "${output_dir}/${build}_exclusion.bed"
```

### Pattern 8: Compact --help Text (DOCS-01)

**What:** Rewrite `show_help()` in cli.sh to be compact (one terminal screen), list only subcommands prominently at top, then flags with one-liner descriptions, and one usage example showing only required args.

**Current problem:** The existing help is verbose (23 lines of option descriptions, each 2-3 sentences). The new help should be ~30 lines maximum.

**Structure:**
```
Usage: hardnormly.sh <subcommand> [options]
       hardnormly.sh run-pipeline -v input.vcf.gz -f ref.fasta [options]

Subcommands:
  run-pipeline             Normalize and filter a VCF file (default)
  generate-inclusion-bed   Merge BED files into a combined inclusion region
  generate-exclusion-bed   Merge BED files into a combined exclusion region

Options (run-pipeline):
  -v, --vcf FILE           Input VCF file (required)
  -f, --fasta FILE         Reference FASTA (required)
  ...
  --caller CALLER          Auto-select filter file: gatk, freebayes
  --filters-file FILE      Filter expression file (overrides --caller)
  --strip-annotations LIST Remove INFO fields before filtering (e.g. INFO/CSQ,INFO/ANN)
  ...

Example:
  hardnormly.sh run-pipeline -v input.vcf.gz -f ref.fasta -o output.vcf.gz

Run 'hardnormly.sh <subcommand> --help' for subcommand-specific options.
```

**Exit code:** `show_help` exits 1 (existing behavior, confirmed by smoke tests).

### Anti-Patterns to Avoid

- **Don't put dispatcher logic in parse_args:** parse_args handles flag parsing, not routing. The dispatcher is a separate function called before parse_args.
- **Don't use `||` on function calls for non-fatal handling:** SC2310/SC2311 is enforced by shellcheck in this project. Use `if ! fn; then ... fi` instead.
- **Don't call `exit` inside parse_args for unknown subcommands:** The dispatcher handles unknown subcommands before parse_args is reached.
- **Don't hardcode defaults/ paths as relative paths:** Use `_SCRIPT_DIR` for caller-file lookups, same pattern as library sourcing.
- **Don't make generate-*-bed subcommands verbose by default:** Decision is silent by default with `-v/--verbose` opt-in.
- **Don't put strip_annotations step after normalization:** Strip first (Step 4.5), normalize second (Step 5). INFO fields should be clean before bcftools norm runs.

---

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| Argument parsing per-subcommand | Full getopts re-implementation | Simple local while/case loop same as parse_args pattern | Consistent style, shellcheck-clean, no external deps |
| BED merging logic | Custom awk merge | `bedtools multiinter` + `bedtools sort` | Already used in bed.sh; handles overlaps correctly |
| Low-complexity region extraction | Custom RepeatMasker parser | `awk '$12 == "Low_complexity"'` filter on rmsk.txt | rmsk.txt column layout is stable; simple awk is sufficient |
| Progress display for downloads | Custom progress bar | `wget -q` (silent) with log_msg milestones | Project uses log_msg pattern; wget progress bar is not needed for large files in pipelines |
| VCF INFO field removal | Custom VCF parser | `bcftools annotate -x` | This is exactly what bcftools annotate is for |
| Caller→filter file resolution | Registry file/config | Inline case statement in hardnormly.sh | Only 2 callers (gatk, freebayes); a case statement is the simplest correct solution |

**Key insight:** Every "new" capability in Phase 5 is wiring through existing tools (bcftools, bedtools, wget). The bash logic is routing and argument parsing only.

---

## Common Pitfalls

### Pitfall 1: SC2310/SC2311 on Non-Fatal Function Calls

**What goes wrong:** Writing `plot_stats_output "$a" "$b" "$c" || log_msg "Warning: failed"` triggers ShellCheck SC2310/SC2311 because `||` on a function call prevents errexit from propagating correctly in `set -Eeuo pipefail` scripts.

**Why it happens:** Under `set -e`, bash treats `|| expr` as "this failure is handled" — but ShellCheck knows this can mask errors in nested function calls. The project has established this as a known pitfall (see STATE.md: SC2310/SC2311 pattern).

**How to avoid:** Use `if ! fn; then ... fi` pattern:
```bash
if ! plot_stats_output "$stats_output" "$plot_output_dir" "$tmp_dir"; then
    log_msg "Warning: plot-vcfstats failed; pipeline continues."
fi
```

**Warning signs:** ShellCheck warnings about SC2310 or SC2311 in CI.

### Pitfall 2: --caller Flag Checked Before _SCRIPT_DIR Is Set

**What goes wrong:** The `--caller` → `filters_file` translation uses `_SCRIPT_DIR` to build the defaults path. But `_SCRIPT_DIR` is set at the top of hardnormly.sh, before `_dispatch` is called. If translation happens inside parse_args (in cli.sh), `_SCRIPT_DIR` is not available.

**How to avoid:** Translation happens in hardnormly.sh (not cli.sh), after `parse_args` returns. cli.sh only sets the `caller` variable. hardnormly.sh does the path lookup.

**Warning signs:** `filters_file` resolves to a relative path that breaks when the script is called from a different directory.

### Pitfall 3: Subcommand Exit Code Propagation

**What goes wrong:** `cmd_generate_inclusion_bed "$@"` returns 0, but the handler calls `exit 0` before cleanup_handler can run. If an error occurs inside the handler, cleanup_handler (EXIT trap) still fires — that's correct. But `exit 0` after success bypasses the rest of hardnormly.sh's pipeline, which is the intended behavior.

**How to avoid:** Always use explicit `exit $?` or `exit 0` in dispatcher routing for subcommands that are complete pipelines. Do NOT let execution fall through to the run-pipeline steps after a successful `generate-*-bed` subcommand.

### Pitfall 4: hg38 Gap Table May Be Empty for Centromeres

**What goes wrong:** The hg38 assembly introduced T2T centromere sequences — the UCSC gap table for hg38 may have zero or few centromere rows. If the download step produces an empty BED, `bedtools multiinter` will still succeed but the centromere exclusion will be missing.

**How to avoid:** In `generate_exclusion_bed.sh`, after downloading the gap table, check row count. Log a warning if zero centromere/telomere rows found. Consider using `cytoBandIdeo.txt.gz` as a fallback for hg38 centromere positions.

**Warning signs:** Empty gaps.bed for hg38 during merge step. Confirmed by inspection of row count after awk filter.

### Pitfall 5: BED Output Missing Header Strip for UCSC Tables

**What goes wrong:** UCSC database text files (genomicSuperDups.txt, rmsk.txt, gap.txt) have a first line starting with `#` or column headers. Passing these directly to `bedtools sort` causes failures.

**How to avoid:** Always pipe through `grep -v '^#'` or `awk 'NR > 1'` before feeding to bedtools. Alternatively, use `tail -n +2` to skip the first line if the format is confirmed header-first.

**Warning signs:** bedtools error: "ERROR: illegal BED3 record" or similar.

### Pitfall 6: generate-*-bed Subcommand Requires Genome File for Slop

**What goes wrong:** `merge_include_beds` requires a genome file (chromosome sizes) for `bedtools slop`. If the user doesn't provide one via `-g`, the subcommand needs to either fetch it or skip slop.

**How to avoid:** Make `-g/--genome` required for `generate-inclusion-bed` (since slop requires it). For `generate-exclusion-bed`, slop is not needed (exclusion regions don't need padding) — so no genome file is required for that subcommand.

---

## Code Examples

Verified patterns from official sources and codebase analysis:

### Subcommand Dispatcher (hardnormly.sh)

```bash
# Source: established pattern from gist.github.com/waylan/4080362 + project SC2310/SC2311 constraints
# Runs BEFORE parse_args in hardnormly.sh
_subcommand="${1:-}"
case "$_subcommand" in
    run-pipeline)
        shift
        parse_args "$@"
        ;;
    generate-inclusion-bed)
        shift
        cmd_generate_inclusion_bed "$@"
        exit 0
        ;;
    generate-exclusion-bed)
        shift
        cmd_generate_exclusion_bed "$@"
        exit 0
        ;;
    "" | -h | --help)
        show_help
        ;;
    -*)
        # Legacy mode: first arg is a flag → run-pipeline implied
        parse_args "$@"
        ;;
    *)
        echo "Unknown subcommand: ${_subcommand}" >&2
        show_help
        ;;
esac
```

### bcftools annotate -x Strip (verified bcftools 1.20 syntax)

```bash
# Source: bcftools 1.20 --help, samtools.github.io/bcftools/bcftools.html#annotate
# Comma-separated list, no spaces. Matches --strip-annotations CLI decision.
bcftools annotate -x "INFO/CSQ,INFO/ANN" input.vcf.gz -Oz -o output.vcf.gz

# In strip_vcf_annotations (lib/annotate.sh):
strip_vcf_annotations() {
    local vcf_file="$1"
    local strip_list="$2"
    local output_vcf="$3"
    run_cmd bcftools annotate -x "$strip_list" "$vcf_file" -Oz -o "$output_vcf"
}
```

### Non-Fatal Plot Error (SC2310/SC2311 safe pattern)

```bash
# Source: established project SC2310/SC2311 pattern (STATE.md)
# In hardnormly.sh, replaces the current direct call to plot_stats_output
if [[ "$plot_stats" == "true" ]]; then
    log_msg "Plotting stats to $plot_output_dir"
    if ! plot_stats_output "$stats_output" "$plot_output_dir" "$tmp_dir"; then
        log_msg "Warning: plot-vcfstats failed; pipeline continues."
    fi
fi
```

### --caller Resolution (hardnormly.sh, after parse_args)

```bash
# Source: codebase analysis + CONTEXT.md decisions
if [[ -n "$caller" ]]; then
    case "$caller" in
        gatk)      _caller_file="${_SCRIPT_DIR}/defaults/gatk_filters.txt" ;;
        freebayes) _caller_file="${_SCRIPT_DIR}/defaults/freebayes_filters.txt" ;;
        *)
            echo "Error: Unknown --caller '$caller'. Valid values: gatk, freebayes" >&2
            exit 1
            ;;
    esac
    if [[ -n "$filters_file" ]]; then
        log_msg "Warning: Both --caller and --filters-file provided; --filters-file takes precedence."
    else
        filters_file="$_caller_file"
    fi
fi
```

### UCSC Table Download for Exclusion BED (helper script pattern)

```bash
# Source: UCSC goldenPath file convention, verified URLs
# For genomicSuperDups (hg19):
wget -q -O - "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/genomicSuperDups.txt.gz" \
    | gunzip -c \
    | grep -v '^#' \
    | awk '{OFS="\t"; print $2, $3, $4}' \
    | bedtools sort -i - \
    > "${tmp_dir}/superdups.bed"

# For rmsk Low_complexity (hg19):
wget -q -O - "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/rmsk.txt.gz" \
    | gunzip -c \
    | grep -v '^#' \
    | awk '$12 == "Low_complexity" {OFS="\t"; print $6, $7, $8}' \
    | bedtools sort -i - \
    > "${tmp_dir}/lcr.bed"

# For gap table centromeres/telomeres (hg19):
wget -q -O - "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/gap.txt.gz" \
    | gunzip -c \
    | grep -v '^#' \
    | awk '$8 ~ /centromere|telomere/ {OFS="\t"; print $2, $3, $4}' \
    | bedtools sort -i - \
    > "${tmp_dir}/gaps.bed"

# For ENCODE blacklist (hg19):
wget -q -O - "https://github.com/Boyle-Lab/Blacklist/raw/master/lists/hg19-blacklist.v2.bed.gz" \
    | gunzip -c \
    | awk '{OFS="\t"; print $1, $2, $3}' \
    | bedtools sort -i - \
    > "${tmp_dir}/blacklist.bed"

# Merge all sources:
bedtools multiinter -i "${tmp_dir}/blacklist.bed" "${tmp_dir}/superdups.bed" \
    "${tmp_dir}/lcr.bed" "${tmp_dir}/gaps.bed" \
    | awk '{OFS="\t"; print $1, $2, $3}' \
    > "${output_dir}/${build}_exclusion.bed"
```

---

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| `hardnormly.sh -v input.vcf ...` (flags only) | `hardnormly.sh run-pipeline -v input.vcf ...` (subcommands) | Phase 5 | Subcommands become primary; old flags-only still works as backward compat |
| No caller shortcut (must provide --filters-file) | `--caller gatk` auto-selects filter file | Phase 5 | Reduces CLI verbosity for standard workflows |
| plot-vcfstats failure exits pipeline | plot-vcfstats failure logged, pipeline continues | Phase 5 | Prevents PDF dependency from breaking VCF output |
| No exclusion BED helper | `scripts/generate_exclusion_bed.sh` downloads from 4 public sources | Phase 5 | Users get curated exclusion BEDs without manual assembly |

**Deprecated/outdated in this phase:**
- `show_help()` current implementation: replaced with compact version including subcommand listing and --caller/--strip-annotations flags.

---

## Open Questions

1. **hg38 centromere gap table coverage**
   - What we know: UCSC gap.txt.gz for hg19 contains centromere/telomere rows (confirmed pattern); hg38 transitioned to T2T assembly which may have fewer gap entries
   - What's unclear: Whether `https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/gap.txt.gz` has centromere rows for hg38
   - Recommendation: In `generate_exclusion_bed.sh`, after downloading hg38 gap table, log the row count. If zero centromere rows, log a warning and use `cytoBandIdeo.txt.gz` col 4 filtering `acen` stain for centromeres as fallback. Plan for this branch in the helper script.

2. **rmsk.txt.gz column index for Low_complexity**
   - What we know: hg19 rmsk.txt column 12 (0-indexed: 11) is the repeat class; "Low_complexity" is the class name; columns 6-8 are chrom/start/end
   - What's unclear: Whether hg38 rmsk.txt has the same column layout
   - Recommendation: Treat as HIGH confidence for hg19 (verified pattern); add a row-count sanity check after awk filter for both builds. If fewer than 1000 rows, emit a warning.

3. **Test coverage for generate-*-bed subcommands in CI**
   - What we know: CI uses apt-get bcftools/bedtools; wget is available; but UCSC/GitHub downloads require internet access
   - What's unclear: Whether download tests should be guarded with a skip (like HAVE_FULL_REF guard)
   - Recommendation: BATS tests for generate-*-bed subcommands should either (a) mock the download step with local test BED files, or (b) skip if no network. Pattern (a) is cleaner — pass local BED files as inputs to test the merge logic, not the download.

4. **DOCS-02 scope (filter file format)**
   - What we know: CONTEXT.md decision is "README only, not --help" — deviates from original DOCS-02 wording
   - What's unclear: Whether the plan should mark DOCS-02 as "implemented in README" or rephrase it
   - Recommendation: Plan should note DOCS-02 is satisfied by README section, not --help. The success criterion in the phase spec ("--help shows filter file format description") conflicts with the CONTEXT.md decision. README is the correct location per user decision.

---

## Sources

### Primary (HIGH confidence)

- bcftools 1.20 installed locally — `bcftools annotate --help` output verified; `-x` flag syntax confirmed
- `lib/cli.sh`, `lib/stats.sh`, `lib/filter.sh`, `lib/annotate.sh`, `lib/bed.sh`, `hardnormly.sh` — direct codebase inspection for all integration points
- `.planning/STATE.md` — SC2310/SC2311 pattern, `_SCRIPT_DIR` convention, all established project decisions
- `.planning/phases/04-refactoring/04-RESEARCH.md` — lib/ module architecture, include guard pattern, run_cmd pattern

### Secondary (MEDIUM confidence)

- [Boyle-Lab/Blacklist GitHub](https://github.com/Boyle-Lab/Blacklist) — ENCODE blacklist confirmed in `lists/` folder; hg19/hg38 URLs from `lists/` directory confirmed
- [UCSC goldenPath hg19/hg38 database](https://hgdownload.soe.ucsc.edu/goldenPath/) — genomicSuperDups.txt.gz, rmsk.txt.gz, gap.txt.gz URL pattern confirmed via WebSearch cross-reference
- [samtools.github.io/bcftools/bcftools.html#annotate](https://samtools.github.io/bcftools/howtos/annotate.html) — bcftools annotate -x examples

### Tertiary (LOW confidence)

- [gist.github.com/waylan/4080362](https://gist.github.com/waylan/4080362) — subcommand dispatcher pattern (cross-verified with project constraints; the exact `sub_${subcommand}` indirect function call approach is NOT used here due to shellcheck SC2086 concerns)
- WebSearch results on bash subcommand patterns — informed dispatcher shape but not used directly

---

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH — no new dependencies; all tools confirmed in existing project
- Architecture: HIGH — patterns derived directly from codebase + established project decisions
- Subcommand dispatcher: HIGH — simple case statement, verified against shellcheck constraints
- bcftools -x syntax: HIGH — confirmed against bcftools 1.20 local installation
- UCSC download URLs: MEDIUM — URL patterns confirmed via WebSearch; actual file column layouts verified for hg19, less certain for hg38
- hg38 gap/centromere coverage: LOW — possible T2T transition issue; flagged as open question

**Research date:** 2026-02-18
**Valid until:** 2026-03-18 (30 days — stable bash/bcftools domain; UCSC URLs rarely change)
