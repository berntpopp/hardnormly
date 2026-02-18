# hardnormly — Issue Review, Prioritization & Roadmap

> Comprehensive review of all [open GitHub issues](https://github.com/berntpopp/hardnormly/issues), augmented with best-practice recommendations for linting, formatting, testing, and modularization.

---

## Executive Summary

hardnormly is a 463-line monolithic bash script with no tests, no linting, no CI/CD, and no modularization. The open issues cluster into four themes:

1. **Code quality** — modularization, linting, formatting (#8, #18)
2. **Testing** — test data and automated tests (#6, #10, #15)
3. **Features** — subcommands, BED generation, plot-vcfstats fix (#5, #11, #12, #13, #14)
4. **Documentation** — help text, exclusion BED docs, filter docs (#16, #17)

The recommended approach applies DRY, KISS, SOLID (adapted for bash), and modularization principles in a dependency-ordered sequence: **infrastructure first, then refactor, then features**.

---

## Current State Assessment

### ShellCheck Results (5 findings)

| Code | Severity | Location | Issue |
|------|----------|----------|-------|
| SC2155 | warning | L74 | `local timestamp=$(...)` — declare and assign separately |
| SC2181 | style | L275, L288, L404, L414 | `$?` checks — use direct exit code testing |

The script is remarkably clean for ShellCheck but has no `.shellcheckrc`, no `shfmt` formatting, and no pre-commit hooks.

### Architecture Concerns

- **Single monolithic file** — all 463 lines in one script with only 3 functions (`cleanup_tmp_dir`, `show_help`, `normalize_bed`); the remaining ~300 lines are top-level procedural code
- **`eval` for pipeline construction** (L401) — works but fragile and hard to test
- **No `set -euo pipefail`** — missing bash strict mode
- **Duplicate cleanup logic** — trap on L28 *and* manual cleanup on L448-453
- **Snakemake workflow** uses `yaml.safe_load(open(...))` without context manager (resource leak)
- **No CI/CD** — no GitHub Actions, no automated checks

---

## Prioritized Issue Roadmap

### Priority 0 — Infrastructure (New, not yet tracked)

These are prerequisites that unblock everything else. Do these first.

#### P0.1: Add ShellCheck + shfmt linting

- **Rationale**: Catches bugs before they reach production; enforces consistent style. [Microsoft recommends ShellCheck for all bash projects](https://microsoft.github.io/code-with-engineering-playbook/code-reviews/recipes/bash/). [GitLab requires it for all shell scripts](https://docs.gitlab.com/development/shell_scripting_guide/).
- **Implementation**:
  1. Add `.shellcheckrc` with project-level config (shell=bash, severity=warning)
  2. Add `.editorconfig` for consistent whitespace
  3. Fix the 5 existing ShellCheck findings (SC2155, SC2181)
  4. Run `shfmt -i 4 -ci -w hardnormly.sh` for consistent formatting
- **Effort**: Small (1-2 hours)

#### P0.2: Add `set -euo pipefail` (bash strict mode)

- **Rationale**: Without this, failed commands silently continue. Critical for a pipeline that chains bcftools/bedtools commands. The `eval` pipeline on L401 is especially risky without `pipefail`.
- **Implementation**: Add after the shebang, with targeted exception handling where needed (e.g., `|| true` for expected failures)
- **Effort**: Small but needs careful testing

#### P0.3: Add GitHub Actions CI

- **Rationale**: No point in having linting if it doesn't run automatically.
- **Implementation**:
  1. `.github/workflows/lint.yml` — runs ShellCheck + shfmt check
  2. `.github/workflows/test.yml` — runs BATS tests (once created)
  3. Use [Snakemake's recommended CI actions](https://snakemake.readthedocs.io/en/stable/snakefiles/best_practices.html) for the workflow: snakemake lint + snakefmt check
- **Effort**: Medium (half day)

#### P0.4: Add pre-commit hooks

- **Rationale**: Catch issues before commit, not in CI.
- **Implementation**: `.pre-commit-config.yaml` with shellcheck, shfmt, snakefmt, and trailing-whitespace hooks
- **Effort**: Small

---

### Priority 1 — Testing Foundation

#### Issue #15 / #6: Add test dataset
**Status**: OPEN (duplicate pair)
**Priority**: **Critical** — blocks all other testing
**Recommendation**:
- Create `tests/data/` with a **minimal synthetic VCF** (10-20 variants covering SNPs, indels, multiallelic sites, and edge cases)
- Include a **tiny reference FASTA** (1-2 chromosomes, small regions) and matching **BED files**
- Include a **pre-computed expected output** VCF for regression testing
- Use [bcftools to generate test VCFs](https://samtools.github.io/bcftools/bcftools.html) or hand-craft them
- Keep data tiny (<100KB total) so CI runs fast
- **Effort**: Medium (half day)
- **Closes**: #6 and #15 together

#### Issue #10: Add BATS tests
**Status**: OPEN
**Priority**: **Critical** — enables safe refactoring
**Recommendation**:
- Use [bats-core](https://github.com/bats-core/bats-core) with bats-support and bats-assert as git submodules
- Structure:
  ```
  tests/
    test_helper/       # bats-support, bats-assert, bats-file submodules
    data/              # test dataset (from #15/#6)
    test_cli.bats      # argument parsing, help, version, error handling
    test_bed.bats      # BED normalization, merge, slop
    test_filter.bats   # filter file parsing, inline filters, PASS filtering
    test_pipeline.bats # end-to-end: input VCF → output VCF matches expected
  ```
- Test categories (ordered by importance):
  1. **Smoke tests** — `--help` exits 0, `--version` prints version, missing args exits non-zero
  2. **Unit tests** — individual functions after modularization (normalize_bed, filter pipeline building)
  3. **Integration tests** — full pipeline run on test dataset, compare output to expected VCF (`bcftools isec` or `diff`)
  4. **Edge cases** — empty VCF, no filters, no BED files, uncompressed output
- **Effort**: Medium-Large (1-2 days)
- **Blocked by**: #15/#6 (test dataset), partially by #8/#18 (easier to test modular functions)

---

### Priority 2 — Refactoring & Modularization

#### Issue #8 / #18: Refactor into modular functions
**Status**: OPEN (duplicate pair — #8 from Aug 2024, #18 from Jan 2025)
**Priority**: **High** — enables testability, readability, and future features
**Recommendation**:

Modularize into a `lib/` directory with sourced function libraries. Each module handles one concern ([DRY](https://en.wikipedia.org/wiki/Don%27t_repeat_yourself), Single Responsibility):

```
lib/
  logging.sh      # log_msg, debug_msg, error handling
  cli.sh          # argument parsing, validation, show_help
  genome.sh       # genome file creation/validation
  bed.sh          # normalize_bed, merge, intersect, compress, index
  annotate.sh     # VCF annotation with BED regions
  normalize.sh    # bcftools norm wrapper
  filter.sh       # filter pipeline construction and execution
  stats.sh        # stats generation and plotting
hardnormly.sh     # main entry point: source libs, orchestrate steps
```

**Key refactoring patterns**:

1. **Extract functions** — every `# Step N:` block becomes a function
2. **Use `local` variables** — avoid global state leaks
3. **Return exit codes** — functions signal success/failure, caller decides what to do
4. **Eliminate duplicate cleanup** — the trap on L28 and the manual cleanup on L448 do the same thing; keep only the trap
5. **Replace `eval` pipeline** — build pipeline as an array or use named pipes; `eval` is error-prone and untestable
6. **Consistent error handling** — every external command should be checked; use a `run_cmd()` wrapper:
   ```bash
   run_cmd() {
       "$@" || { log_msg "Error: $1 failed"; exit 1; }
   }
   ```

**Effort**: Large (2-3 days)
**Closes**: #8 and #18 together

---

### Priority 3 — Bug Fixes

#### Issue #12: Error in plot-vcfstats (tectonic PDF failure)
**Status**: OPEN (bug)
**Priority**: **Medium** — affects stats output but not core filtering
**Root cause**: `tectonic` fails on some LaTeX generated by `plot-vcfstats`. This is a known issue with plot-vcfstats + tectonic compatibility.
**Recommendation**:
1. Make plotting **non-fatal** — log the error but don't `exit 1` (currently L434-438 kills the pipeline for a cosmetic feature)
2. Add `--plot-backend` option to choose between tectonic and pdflatex
3. Consider replacing `plot-vcfstats` with a direct Python/matplotlib solution for better control
4. Document the tectonic requirement and workaround in README
- **Effort**: Small-Medium

---

### Priority 4 — Feature Enhancements

#### Issue #5 / #14: Add subcommands for specific steps
**Status**: OPEN (duplicate pair)
**Priority**: **Medium** — nice-to-have, depends on modularization
**Recommendation**:
- After modularization (#8/#18), add a subcommand dispatcher:
  ```bash
  case "${1:-}" in
      generate-inclusion-bed)  shift; generate_inclusion_bed "$@" ;;
      generate-exclusion-bed)  shift; generate_exclusion_bed "$@" ;;
      run-pipeline|"")         shift; run_pipeline "$@" ;;
      *)                       show_help ;;
  esac
  ```
- Each subcommand sources only the libs it needs (KISS — don't load everything)
- Maintain backward compatibility: running without subcommand = `run-pipeline`
- **Effort**: Medium (after modularization is done)
- **Blocked by**: #8/#18

#### Issue #13: Add option to provide genome file and skip generation
**Status**: OPEN
**Priority**: **Low** — **already implemented!**
**Recommendation**: Close this issue. The `-g/--genome` flag already exists (L104-106) and the script skips generation when a genome file is provided (L187-194). Just verify with a test and close.

#### Issue #11: Generate combined exclusion BED files for hg19 and hg38
**Status**: OPEN
**Priority**: **Low** — useful but independent of core pipeline
**Recommendation**:
- Create a helper script `scripts/generate_exclusion_beds.sh` that:
  1. Downloads exclusion regions from [excluderanges](https://github.com/dozmorovlab/excluderanges) and [ENCODE](https://hgdownload.cse.ucsc.edu/goldenpath/hg19/encodeDCC/wgEncodeMapability)
  2. Merges and sorts them with bedtools
  3. Outputs `defaults/hg19_exclusions.bed` and `defaults/hg38_exclusions.bed`
- Include pre-built files in the repo for convenience
- **Effort**: Medium

---

### Priority 5 — Documentation

#### Issue #17: Enhance show_help()
**Status**: OPEN
**Priority**: **Low** — help text is already decent after previous enhancement
**Recommendation**:
- Add examples section to `--help` output
- Add `--filters-file` format description in help text
- Consider a `--help-filters` subcommand that prints the filter file format
- **Effort**: Small
- **Note**: May already be partially addressed; current help (L44-70) is reasonably detailed

#### Issue #16: Document exclusion BED file sources
**Status**: OPEN
**Priority**: **Low**
**Recommendation**:
- Already partially addressed in README (L93-124)
- Add a `defaults/exclusion_beds.md` documenting each source BED file, its origin URL, genome build, and what regions it covers
- **Effort**: Small

---

## New Issues to Create

Beyond the existing open issues, the following should be tracked:

### N1: Add `set -euo pipefail` and strict error handling
- **Why**: Pipeline can silently produce corrupt output on intermediate failures
- **Maps to**: DRY principle (single error handling pattern), defensive programming

### N2: Add ShellCheck + shfmt + pre-commit + CI
- **Why**: No automated quality checks exist; all linting/formatting is manual
- **Maps to**: Industry standard for bash projects per [Microsoft](https://microsoft.github.io/code-with-engineering-playbook/code-reviews/recipes/bash/) and [GitLab](https://docs.gitlab.com/development/shell_scripting_guide/) engineering playbooks

### N3: Add snakefmt + snakemake --lint to CI
- **Why**: Snakemake workflow has no formatting/linting either
- **Maps to**: [Snakemake best practices](https://snakemake.readthedocs.io/en/stable/snakefiles/best_practices.html) recommend lint and format checks

### N4: Replace `eval` pipeline construction with safer alternative
- **Why**: `eval` on L401 executes a dynamically built string — hard to debug, test, or secure
- **Alternative**: Use a temporary script file, or build pipeline with process substitution

### N5: Fix Snakemake workflow to use context managers
- **Why**: `yaml.safe_load(open('config.yaml'))` and `open(config['vcf_files'])` leak file handles
- **Fix**: Use `with open(...) as f:` pattern

### N6: Add `--caller` flag to auto-select filter file
- **Why**: Users must manually match filter files to callers; error-prone
- **Implementation**: `--caller gatk` → uses `defaults/gatk_filters.txt`; `--caller freebayes` → uses `defaults/freebayes_filters.txt`

---

## Recommended Implementation Order

```
Phase 1: Infrastructure (P0)          ~1 day
  P0.1  ShellCheck + shfmt
  P0.2  set -euo pipefail
  P0.3  GitHub Actions CI
  P0.4  pre-commit hooks
  N5    Fix Snakemake file handle leaks
         │
Phase 2: Test Foundation (P1)          ~1-2 days
  #15/#6  Create test dataset
  #10     Add BATS tests (smoke + integration)
  #13     Close as already implemented (verify with test)
         │
Phase 3: Refactor (P2)                ~2-3 days
  #8/#18  Modularize into lib/ functions
  N1      Strict error handling
  N4      Replace eval pipeline
          Add unit tests for each module
         │
Phase 4: Features (P3-P4)             ~2-3 days
  #12     Fix plot-vcfstats (make non-fatal)
  #5/#14  Add subcommands
  #11     Generate combined exclusion BEDs
  N6      Add --caller flag
         │
Phase 5: Documentation (P5)           ~0.5 day
  #17     Enhance show_help()
  #16     Document exclusion BED sources
  N3      Snakemake lint/format in CI
```

---

## Summary Table

| # | Title | State | Priority | Effort | Blocked By | Recommendation |
|---|-------|-------|----------|--------|------------|----------------|
| — | ShellCheck + shfmt + CI | NEW | P0 | S | — | Add linting infra first |
| — | `set -euo pipefail` | NEW | P0 | S | — | Essential safety net |
| — | GitHub Actions CI | NEW | P0 | M | — | Automate all checks |
| — | pre-commit hooks | NEW | P0 | S | — | Catch issues locally |
| #15/#6 | Add test dataset | OPEN | P1 | M | — | Synthetic minimal VCF + BED |
| #10 | Add BATS tests | OPEN | P1 | L | #15 | Smoke, unit, integration tests |
| #13 | Genome file option | OPEN | P1 | — | — | **Close: already implemented** |
| #8/#18 | Modularize | OPEN | P2 | L | #10 | lib/ directory, source pattern |
| — | Replace `eval` | NEW | P2 | M | #8 | Safer pipeline construction |
| #12 | plot-vcfstats bug | OPEN | P3 | M | — | Make non-fatal |
| #5/#14 | Subcommands | OPEN | P4 | M | #8 | Dispatcher pattern |
| #11 | Exclusion BED gen | OPEN | P4 | M | — | Helper script + pre-built files |
| — | `--caller` flag | NEW | P4 | S | — | Auto-select filter file |
| #17 | Enhance help | OPEN | P5 | S | — | Add examples, filter format |
| #16 | Document excl. BEDs | OPEN | P5 | S | — | Expand existing README section |

**Legend**: S = Small (<2h), M = Medium (half day), L = Large (1-3 days)

---

## Key Principles Applied

- **DRY**: Deduplicate cleanup logic; shared functions in lib/ sourced by both main script and subcommands; single error handling pattern
- **KISS**: Don't over-abstract — bash functions, not a framework; simple subcommand dispatcher, not a CLI framework; test with BATS, not a custom harness
- **SOLID** (adapted for bash): Single Responsibility per module (bed.sh does BED ops only); Open/Closed via filter files (add filters without modifying code); Dependency Inversion via sourced libraries (swap implementations)
- **Modularization**: lib/ directory with focused modules; main script becomes orchestrator; each module independently testable

---

## References

- [ShellCheck — static analysis for shell scripts](https://github.com/koalaman/shellcheck)
- [shfmt — shell formatter](https://github.com/mvdan/sh)
- [Microsoft Bash Code Review Playbook](https://microsoft.github.io/code-with-engineering-playbook/code-reviews/recipes/bash/)
- [GitLab Shell Scripting Standards](https://docs.gitlab.com/development/shell_scripting_guide/)
- [BATS — Bash Automated Testing System](https://github.com/bats-core/bats-core)
- [Snakemake Best Practices](https://snakemake.readthedocs.io/en/stable/snakefiles/best_practices.html)
- [Designing Modular Bash: Functions, Namespaces, and Library Patterns](https://www.lost-in-it.com/posts/designing-modular-bash-functions-namespaces-library-patterns/)
- [bcftools Documentation](https://samtools.github.io/bcftools/bcftools.html)
- [AMP/CAP Bioinformatics Pipeline Validation Guidelines](https://www.sciencedirect.com/science/article/pii/S1525157817303732)
