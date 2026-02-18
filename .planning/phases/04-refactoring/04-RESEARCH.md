# Phase 4: Refactoring - Research

**Researched:** 2026-02-18
**Domain:** Bash script modularization — lib/ modules, include guards, error handling wrappers, cross-module patterns
**Confidence:** MEDIUM (bash has no authoritative spec like Context7 libraries; findings are verified across multiple credible sources)

## Summary

Phase 4 splits hardnormly.sh (544 lines) into 8 focused lib/ modules and replaces the current `trap 'err_handler ${LINENO}' ERR` pattern with an explicit `run_cmd` wrapper function. The research domain is bash script architecture — there is no single authoritative source like a library's official docs. Findings were cross-verified across credible sources (betterdev.blog, gabrielstaples.com, wafaicloud.com, sap1ens.com, bats-core docs) plus direct codebase analysis.

The standard approach for bash modularization is: explicit `source` calls using `$(dirname "${BASH_SOURCE[0]}")` for path resolution, one include guard per module, logging as a global-variable module (the one exception to "avoid globals"), all other state passed as function arguments, and a `run_cmd` wrapper that captures the command, its exit code, and stderr before printing a rich error and propagating the failure.

Filter parsing ownership belongs in `cli.sh` (CLI validates and tokenizes) while `filter.sh` owns the execution mechanism. This cleanly separates "what to filter" from "how to apply filters." The PASS filter and output format step belongs in `filter.sh` because it is the last filter stage in the pipeline, not a separate orchestration concern.

**Primary recommendation:** Explicit source lines in hardnormly.sh using `SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)`, include guards in every lib/ module, logging as global-variable module (special case), all other modules use function arguments only, `run_cmd` wrapper replaces the ERR trap for external commands.

---

## Standard Stack

This phase has no new library dependencies. All tooling was established in prior phases.

### Tools (unchanged from prior phases)

| Tool | Version | Purpose | Notes |
|------|---------|---------|-------|
| bash | 4+ | Script runtime | Requires arrays, `[[ ]]`, `BASH_SOURCE` |
| shfmt | current | Code formatting | Flags: `-i 0 -bn -ci` (locked in Phase 1) |
| shellcheck | current | Static analysis | Already in use |
| bats-core | 1.x | Unit tests for lib/ modules | Already installed (Phase 3) |
| bats-assert | 2.1.0 | Assertions in unit tests | Already installed (Phase 3) |

### No New Dependencies

The refactoring is structural — no new libraries are introduced. bcftools, bedtools, bgzip, tabix remain the external tools called through `run_cmd`.

---

## Architecture Patterns

### Recommended Project Structure

```
hardnormly/
├── hardnormly.sh            # Orchestrator: sources lib/, parses args, calls pipeline functions
└── lib/
    ├── logging.sh           # log_msg, debug_msg, error_msg, run_cmd; global _LOG_FILE, _DEBUG
    ├── cli.sh               # parse_args, validate_args, show_help, show_version
    ├── genome.sh            # create_genome_file, validate_genome_file
    ├── bed.sh               # normalize_bed, merge_beds, compress_index_bed, create_header_file
    ├── annotate.sh          # annotate_vcf_with_regions
    ├── normalize.sh         # normalize_vcf
    ├── filter.sh            # build_filter_stages, apply_filter_stages, write_output
    └── stats.sh             # generate_stats, plot_stats
```

**Note on bed/annotate split:** `bed.sh` handles all BED file operations including creating the `.hdr` header files (because headers are a direct output of BED preparation, not a VCF concern). `annotate.sh` receives ready-compressed BED files and header paths as arguments and runs `bcftools annotate`. This split is clean: bed.sh knows about BED regions, annotate.sh knows about bcftools VCF annotation.

### Pattern 1: Include Guard (every module)

**What:** Each `lib/*.sh` module checks a sentinel variable on source; returns immediately if already loaded.

**When to use:** Every module, unconditionally.

**Source:** Verified against wafaicloud.com (multiple sources agree on this exact pattern).

```bash
# lib/logging.sh — include guard at top of every module
[[ -n "${_LIB_LOGGING_LOADED:-}" ]] && return 0
readonly _LIB_LOGGING_LOADED=1
```

**Key detail:** Use `${_VAR:-}` (not `${_VAR}`) in the guard check so that `set -u` doesn't abort when the variable is unset on first source. Use `readonly` after setting to prevent accidental overwrite.

**Naming convention:** `_LIB_<MODULENAME>_LOADED` — underscore prefix signals internal-only, module name in uppercase, avoids collision across modules.

### Pattern 2: Module Loading in hardnormly.sh

**What:** Main script resolves its own directory via `BASH_SOURCE[0]`, then sources each lib module with an explicit relative path.

**When to use:** The orchestrator script only. Modules do NOT source each other; the orchestrator loads all modules.

**Source:** gabrielstaples.com (verified with multiple sources — explicit over auto-load is the consensus for maintainability).

```bash
# hardnormly.sh — module loading section (near top, after shebang and set -Eeuo pipefail)
_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${_SCRIPT_DIR}/lib/logging.sh"
source "${_SCRIPT_DIR}/lib/cli.sh"
source "${_SCRIPT_DIR}/lib/genome.sh"
source "${_SCRIPT_DIR}/lib/bed.sh"
source "${_SCRIPT_DIR}/lib/annotate.sh"
source "${_SCRIPT_DIR}/lib/normalize.sh"
source "${_SCRIPT_DIR}/lib/filter.sh"
source "${_SCRIPT_DIR}/lib/stats.sh"
```

**Why explicit over auto-loader:** An auto-loader (scanning lib/ and sourcing everything) adds complexity for no gain in a fixed 8-module system. Explicit source lines are self-documenting, grep-able, and immediately obvious to any reader. Load order is visible at a glance.

**Why `BASH_SOURCE[0]` not `$0`:** `$0` is the basename only when bash is invoked certain ways. `BASH_SOURCE[0]` is always the path to the currently executing script. Use `realpath` or the `cd && pwd` idiom (both work; the `cd && pwd` approach avoids a dependency on `realpath` being available).

### Pattern 3: Logging as Global-Variable Module (justified exception)

**What:** `logging.sh` sets two global variables that all modules read directly: `_LOG_FILE` and `_DEBUG`. All other modules use function arguments only.

**Why this exception is justified:**

Logging is universally accepted as the correct use case for global state in bash. The alternatives are worse:
- Passing `log_file` and `debug` as arguments to every function in every module creates enormous noise (every function signature changes) with no architectural benefit
- A return-value or output-capture approach for logging is circular (you'd need to log errors in the logging function itself)

In production bash tooling, the pattern is: globals for logger configuration (log destination, verbosity), arguments for everything else. This is analogous to Python's `logging.basicConfig()` — set once, used everywhere.

```bash
# lib/logging.sh — globals for logger config (the ONE module that uses globals)
_LOG_FILE=""   # set by cli.sh after arg parsing via: set_log_file "$path"
_DEBUG=false   # set by cli.sh after arg parsing via: set_debug true

set_log_file() { _LOG_FILE="$1"; }
set_debug()    { _DEBUG="$1"; }

log_msg() {
    local timestamp
    timestamp=$(date +"%Y-%m-%d %H:%M:%S")
    local msg="[$timestamp] $1"
    if [[ -n "$_LOG_FILE" ]]; then
        printf '%s\n' "$msg" >> "$_LOG_FILE"
    else
        printf '%s\n' "$msg"
    fi
}

debug_msg() {
    [[ "$_DEBUG" == "true" ]] || return 0
    log_msg "[DEBUG] $1"
}

error_msg() {
    log_msg "ERROR: $1" >&2
}
```

### Pattern 4: run_cmd Wrapper (replaces ERR trap)

**What:** All external tool calls go through `run_cmd`, which captures stderr to a temp file, checks exit code, and emits a rich error message before failing.

**When to use:** Every call to bcftools, bedtools, bgzip, tabix, mysql, plot-vcfstats. NOT for bash built-ins or simple test commands.

**Source:** Synthesized from betterdev.blog error handling template, multiple bash scripting guides on stderr capture patterns; the specific design decisions below are informed by the project's locked requirements.

```bash
# lib/logging.sh — run_cmd lives here because it uses log_msg/error_msg
# Usage: run_cmd [description] -- cmd [args...]
# Or:    run_cmd cmd [args...]
run_cmd() {
    local description=""
    local -a cmd=()

    # Optional description before "--"
    if [[ "$1" == "--" ]]; then
        shift
        cmd=("$@")
    elif [[ "${2:-}" == "--" ]]; then
        description="$1"
        shift 2
        cmd=("$@")
    else
        cmd=("$@")
    fi

    local stderr_file
    stderr_file=$(mktemp)

    local exit_code=0
    "${cmd[@]}" 2>"$stderr_file" || exit_code=$?

    if [[ "$exit_code" -ne 0 ]]; then
        local stderr_content
        stderr_content=$(cat "$stderr_file")
        rm -f "$stderr_file"
        local label="${description:-${cmd[0]}}"
        error_msg "'${label}' failed (exit ${exit_code})"
        error_msg "Command: ${cmd[*]}"
        if [[ -n "$stderr_content" ]]; then
            error_msg "Stderr: ${stderr_content}"
        fi
        return "$exit_code"
    fi

    rm -f "$stderr_file"
    return 0
}
```

**Why remove the ERR trap:** The locked decision is to replace `trap 'err_handler ${LINENO}' ERR` with `run_cmd`. The ERR trap fires on any command exit code != 0 anywhere in the script including in test conditions, `||` chains, and `if` blocks — which creates false positives and requires careful scoping. With `run_cmd`, only intentional external commands trigger the rich error message. The EXIT trap (cleanup_handler) is kept unchanged — it is responsible for temp dir cleanup, not error reporting.

**Retry for network operations (UCSC MySQL only):**

```bash
# lib/genome.sh — retry wrapper for MySQL only
run_cmd_with_retry() {
    local max_attempts="${1}"
    shift
    local attempt=1
    while [[ "$attempt" -le "$max_attempts" ]]; do
        run_cmd "$@" && return 0
        log_msg "Attempt ${attempt}/${max_attempts} failed. Retrying..."
        ((attempt++))
        sleep 2
    done
    error_msg "All ${max_attempts} attempts failed."
    return 1
}
```

### Pattern 5: Function Arguments (all non-logging modules)

**What:** Module functions receive all state they need as explicit arguments. No global state reads (except `_LOG_FILE`/`_DEBUG` from logging.sh).

**When to use:** genome.sh, bed.sh, annotate.sh, normalize.sh, filter.sh, stats.sh functions.

**Source:** wafaicloud.com, sap1ens.com, bertvv.github.io/cheat-sheets/Bash.html — consistent consensus across sources.

```bash
# lib/bed.sh example — explicit argument passing
normalize_bed() {
    local bed_file="$1"
    local annotation="$2"
    local output_file="$3"
    debug_msg "Normalizing BED file: ${bed_file} with annotation: ${annotation}"
    awk -v annot="$annotation" '{OFS="\t"; print $1, $2, $3, annot}' "$bed_file" \
        | bedtools sort -i - > "$output_file" \
        || { error_msg "normalize_bed failed for ${bed_file}"; return 1; }
    debug_msg "Normalized BED file written to: ${output_file}"
}
```

**Why not globals:** Globals make functions untestable in isolation (BATS can't control the global state without side effects on other tests). Function arguments make BATS unit tests trivial: call the function with known inputs, assert the output file contents. REFR-13 requires unit tests per module.

### Pattern 6: Init Function — Ready-on-Source (no explicit init needed)

**What:** Modules define functions and set include guard — no `init()` function required. logging.sh initializes its globals to empty/false defaults at source time.

**When to use:** All modules.

**Rationale:** Modules in this project don't need initialization work beyond defining functions. They have no connections to open, no files to pre-create, no services to start. An `init()` pattern adds boilerplate with no payoff. The only initialization needed (setting log file path, enabling debug) happens via the `set_log_file`/`set_debug` setter functions in logging.sh, called once by hardnormly.sh after argument parsing.

### Pattern 7: Filter Parsing Ownership (cli.sh → filter_stages array passed to filter.sh)

**What:** `cli.sh` parses `--filters` / `--filters-file` arguments and builds the `filter_stages` array (the same `name|action|expression` format already in use). This array is passed to `filter.sh`'s apply function.

**Why cli.sh owns parsing:** CLI parsing is the responsibility of the CLI module. Filters are just another CLI argument type. The existing format (`name action expression` → `name|action|expression`) is a parsing/tokenization concern, not a filter-execution concern.

**Why filter.sh does NOT re-parse:** filter.sh owns the bcftools filter execution loop. It receives ready-to-execute `filter_stages` array entries. This is clean separation: cli.sh knows about user input format, filter.sh knows about bcftools mechanics.

```bash
# cli.sh: parse filter arguments → populate global array
# Note: filter_stages is RETURNED via nameref or populated in caller's scope
# Pattern: pass filter_stages as nameref argument (bash 4.3+)
parse_filter_args() {
    local -n _stages_ref="$1"   # nameref to caller's array
    shift
    local filters_file="$1"
    shift
    local -a inline_filters=("$@")

    # Region-based filters come from BED processing — not parsed here
    # Inline filters from --filters args
    for filter in "${inline_filters[@]}"; do
        local name action expr
        IFS=" " read -r name action expr <<< "$filter"
        _stages_ref+=("${name}|${action}|${expr}")
    done

    # File-based filters
    if [[ -n "$filters_file" ]]; then
        while IFS=" " read -r name action expr; do
            expr=$(tr -d '\r\n' <<< "$expr")
            _stages_ref+=("${name}|${action}|${expr}")
        done < "$filters_file"
    fi
}
```

**Alternative considered:** cli.sh populates an exported global `FILTER_STAGES` array, filter.sh reads it. Rejected: global arrays make testing harder and break the "explicit arguments" principle.

### Pattern 8: PASS Filter and Output Format (owned by filter.sh)

**What:** The final `bcftools view -f PASS` step and output format detection (`-Oz` / `-Ov`) belong in `filter.sh`'s `write_filtered_output` function.

**Rationale:** PASS filtering IS a filter stage — the last one. Output format is a concern of the step that writes the final VCF. Both belong together in filter.sh. Putting them in hardnormly.sh main makes filter.sh's interface awkward (filter.sh produces a BCF temp file, main has to know about it).

**Interface:** `write_filtered_output <current_bcf> <output_vcf> <only_pass> <auto_index>`

### Anti-Patterns to Avoid

- **Global variables for pipeline state:** Don't make `tmp_dir`, `vcf_file`, `fasta_file`, etc. global. Pass them as arguments. Globals survive between function calls in BATS tests and cause test pollution.
- **Modules sourcing other modules:** Don't have filter.sh source logging.sh. All modules are sourced by hardnormly.sh. Cross-module source creates ordering dependencies and double-source risks.
- **Sourcing with relative paths from CWD:** `source lib/logging.sh` breaks when the script is called from a different directory. Always use `${_SCRIPT_DIR}/lib/module.sh`.
- **Shadowing `$BASH_SOURCE` in modules:** Don't use `BASH_SOURCE[0]` inside modules for self-location — it resolves to the module file, not the main script. Only hardnormly.sh uses `BASH_SOURCE[0]` for `_SCRIPT_DIR`.
- **ERR trap on subshells:** The existing `err_handler` fires in subshells and command substitutions. The replacement `run_cmd` avoids this entirely by checking exit code inline.
- **Skipping the include guard:** Without include guards, sourcing modules in test setup AND having the module source itself again creates duplicate function definitions and can corrupt state.

---

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| Bash argument parsing | Custom getopts wrapper | The existing while/case loop in cli.sh | The existing pattern works, is tested, and shfmt-compliant; no external tool needed |
| Retry logic | Complex retry framework | Simple `while attempt <= max` loop | Network retry for MySQL is 2-3 attempts max; a loop is sufficient |
| Module dependency tracking | Import graph / topological sort | Explicit source order in hardnormly.sh | 8 fixed modules with no inter-module dependencies; a graph is over-engineering |
| Stderr capture | Named pipes / process substitution | mktemp temp file in run_cmd | Temp file is simple, reliable, works with set -e; process substitution has FIFO ordering issues |
| Log rotation | logrotate integration | Simple file append | hardnormly.sh is a one-shot tool, not a daemon; log rotation is not needed |

**Key insight:** In bash, the "don't hand-roll" principle mainly applies to avoiding re-implementation of bcftools/bedtools capabilities in pure bash. For the structural refactoring itself, simple patterns beat frameworks.

---

## Common Pitfalls

### Pitfall 1: `BASH_SOURCE` Path Resolution Breaks When Symlinked

**What goes wrong:** `dirname "${BASH_SOURCE[0]}"` returns the symlink's directory, not the real script directory. If lib/ is relative to the real script, sourcing fails.

**Why it happens:** `BASH_SOURCE[0]` gives the path as invoked, which may be through a symlink.

**How to avoid:** Use `cd "$(dirname "${BASH_SOURCE[0]}")" && pwd` — this resolves symlinks because `pwd` returns the physical path after `cd`. Or use `realpath "${BASH_SOURCE[0]}"` if available.

**Warning signs:** Tests pass when run from project root but fail when hardnormly.sh is called via a symlink in PATH.

### Pitfall 2: Include Guard Variable Collides Across Modules

**What goes wrong:** Two modules both define `_LIB_LOADED=1`. Second module never loads.

**Why it happens:** Generic guard name without module namespace.

**How to avoid:** Use `_LIB_<MODULENAME>_LOADED` naming. Example: `_LIB_BED_LOADED`, `_LIB_FILTER_LOADED`.

### Pitfall 3: `set -u` Fires on Unset Include Guard Variable

**What goes wrong:** First time a module is sourced, `[[ -n "$_LIB_BED_LOADED" ]]` causes `unbound variable` error under `set -u`.

**Why it happens:** The variable doesn't exist yet, and `set -u` treats unset variables as errors.

**How to avoid:** Always use the `:-` default expansion: `[[ -n "${_LIB_BED_LOADED:-}" ]]`. The `:-` provides an empty default, satisfying `set -u`.

### Pitfall 4: `run_cmd` Temp File Not Cleaned on Signal Interruption

**What goes wrong:** `run_cmd` creates a temp file for stderr. If a SIGINT arrives mid-command, the temp file leaks.

**Why it happens:** The mktemp temp file inside run_cmd isn't registered with the EXIT trap.

**How to avoid:** The existing `cleanup_handler` already removes `$tmp_dir`. Register `run_cmd`'s stderr temp file in `$tmp_dir` rather than in `/tmp` directly:
```bash
stderr_file=$(mktemp "${tmp_dir}/runcmd-stderr-XXXXXX")
```
This way, `cleanup_handler`'s `rm -rf "$tmp_dir"` covers it.

### Pitfall 5: filter_stages Array Scoping in BATS Tests

**What goes wrong:** BATS test sources filter.sh and calls `apply_filter_stages`. But the array built in cli.sh's `parse_filter_args` (using nameref) isn't visible because tests call functions directly.

**Why it happens:** Bash namerefs only work when the referenced variable exists in the caller's scope.

**How to avoid:** In BATS unit tests for filter.sh, construct a `filter_stages` array directly in the test and pass it to the function. Document that `parse_filter_args` is tested separately in cli.sh unit tests.

### Pitfall 6: `bcftools norm` Stderr Contains Warnings Not Errors

**What goes wrong:** `run_cmd` captures all stderr as an "error" message. bcftools norm emits `Warning:` lines on stderr even for successful runs.

**Why it happens:** bcftools uses stderr for both warnings and errors; exit code 0 means success regardless of warnings.

**How to avoid:** In `normalize.sh`, do NOT use bare `run_cmd`. Instead, use the existing pattern: capture stderr separately, check exit code, log warnings on success, log errors on failure. This is the one case where the current inline stderr-capture approach is the right design. The `run_cmd` wrapper is for commands where any stderr output + non-zero exit = error.

**Specific to this codebase:** The current hardnormly.sh already handles this correctly (lines 386-412). Preserve that pattern in normalize.sh rather than routing through run_cmd.

### Pitfall 7: Global Array `filter_stages` Persists Between BATS Tests

**What goes wrong:** BATS runs setup() before each test. If `filter_stages` is a global array populated by parse_filter_args, tests that don't reset it inherit entries from previous tests.

**Why it happens:** Global arrays aren't cleared between BATS tests (only local variables go away).

**How to avoid:** Use the nameref pattern in cli.sh so the caller controls the array lifecycle. In tests, declare a fresh `local filter_stages=()` per test.

### Pitfall 8: Duplicate Cleanup Logic (REFR-12)

**What goes wrong:** cleanup_handler in hardnormly.sh plus ad-hoc `rm -f` calls at error sites both try to remove the same files.

**Why it happens:** Pre-refactor code has `rm -f "$norm_output" "$norm_stdout"` inline (lines 395, 410) AND cleanup_handler removes `$tmp_dir`. After refactoring, inline rm calls should be removed — cleanup_handler handles everything.

**How to avoid:** After refactoring, all temp files must live in `$tmp_dir`. cleanup_handler removes `$tmp_dir` on EXIT. No individual `rm -f` calls at error sites.

---

## Code Examples

### Include Guard (canonical form)

```bash
# Every lib/*.sh file starts with this — no exceptions
[[ -n "${_LIB_<MODULENAME>_LOADED:-}" ]] && return 0
readonly _LIB_<MODULENAME>_LOADED=1
```

### Module Loading in hardnormly.sh

```bash
#!/bin/bash
# hardnormly.sh

set -Eeuo pipefail

_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

source "${_SCRIPT_DIR}/lib/logging.sh"
source "${_SCRIPT_DIR}/lib/cli.sh"
source "${_SCRIPT_DIR}/lib/genome.sh"
source "${_SCRIPT_DIR}/lib/bed.sh"
source "${_SCRIPT_DIR}/lib/annotate.sh"
source "${_SCRIPT_DIR}/lib/normalize.sh"
source "${_SCRIPT_DIR}/lib/filter.sh"
source "${_SCRIPT_DIR}/lib/stats.sh"
```

### run_cmd Pattern (core)

```bash
# Simplified production version in lib/logging.sh
run_cmd() {
    local stderr_file
    stderr_file=$(mktemp "${tmp_dir}/runcmd-XXXXXX")
    local exit_code=0
    "$@" 2>"$stderr_file" || exit_code=$?
    if [[ "$exit_code" -ne 0 ]]; then
        local stderr_content
        stderr_content=$(< "$stderr_file")
        rm -f "$stderr_file"
        error_msg "'$1' failed (exit ${exit_code})"
        error_msg "Command: $*"
        [[ -n "$stderr_content" ]] && error_msg "Stderr: ${stderr_content}"
        return "$exit_code"
    fi
    rm -f "$stderr_file"
}

# Usage in lib/annotate.sh:
run_cmd bcftools annotate \
    -a "$include_bed_gz" \
    -h "$include_hdr" \
    -c CHROM,FROM,TO,INCLUDE_REGION \
    "$vcf_file" -Oz -o "$output_vcf"
```

### BATS Unit Test for a lib/ Module Function

```bash
# tests/05-lib-bed.bats (example for REFR-13)
setup() {
    load 'test_helper/common'
    _common_setup
    # Source just the module under test
    source "${REPO_ROOT}/lib/logging.sh"
    source "${REPO_ROOT}/lib/bed.sh"
}

@test "normalize_bed adds annotation column and sorts" {
    local input="${BATS_TEST_TMPDIR}/input.bed"
    local output="${BATS_TEST_TMPDIR}/output.bed"
    printf '22\t100\t200\n22\t50\t75\n' > "$input"

    normalize_bed "$input" "INCLUDE" "$output"

    assert_success
    run grep -c "INCLUDE" "$output"
    assert_output "2"
}
```

### Nameref Pattern for filter_stages (bash 4.3+)

```bash
# lib/cli.sh
parse_inline_filters() {
    local -n _ref="$1"   # nameref to caller's array
    shift
    local filters_file="$1"
    shift

    for filter in "$@"; do
        local name action expr
        IFS=" " read -r name action expr <<< "$filter"
        _ref+=("${name}|${action}|${expr}")
    done

    if [[ -n "$filters_file" ]]; then
        while IFS=" " read -r name action expr; do
            expr=$(tr -d '\r\n' <<< "$expr")
            _ref+=("${name}|${action}|${expr}")
        done < "$filters_file"
    fi
}

# hardnormly.sh caller:
filter_stages=()
parse_inline_filters filter_stages "$filters_file" "${filters[@]}"
# Region-based entries prepended separately (after BED processing)
```

---

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|---|---|---|---|
| Single monolithic script | lib/ modules sourced by orchestrator | Phase 4 | Testability, maintainability |
| `trap 'err_handler ${LINENO}' ERR` | `run_cmd` wrapper per external call | Phase 4 | Rich error messages, no false positives |
| Inline `rm -f` at error sites + cleanup_handler | cleanup_handler only (all temp in $tmp_dir) | Phase 4 | Single cleanup path (REFR-12) |
| Global variable soup (all state in globals) | Arguments for pipeline state, globals for logging config only | Phase 4 | Module testability |
| eval-based filter pipeline | Array-based filter_stages with | delimiter | Phase 1 (already done) | N/A for Phase 4 |

**Already resolved from roadmap:**
- "Replace eval" — eval was removed in Phase 1 (REFR-10 criterion already satisfied; no Phase 4 work needed here)

---

## Open Questions

1. **bash version constraint for `local -n` (nameref)**
   - What we know: `local -n` (namerefs) requires bash 4.3+. Ubuntu 20.04+ ships bash 5.x. macOS ships bash 3.x by default.
   - What's unclear: Does the project need to support macOS default bash?
   - Recommendation: If macOS support needed, replace nameref pattern with a global `_FILTER_STAGES` array (set by cli.sh, consumed by filter.sh). If macOS support not needed (bioinformatics tool, likely Linux-only), use nameref pattern. Check CLAUDE.md — it mentions Ubuntu/WSL2 environments. The project appears Linux-targeted; nameref should be safe.

2. **`tmp_dir` availability in run_cmd**
   - What we know: `run_cmd` should write stderr temp files to `$tmp_dir` so cleanup_handler covers them. But `tmp_dir` is set in hardnormly.sh (main script scope).
   - What's unclear: Should run_cmd read `tmp_dir` as a global (special case), or should callers pass a temp dir argument to run_cmd?
   - Recommendation: Make `tmp_dir` a module-level global in logging.sh (set via `set_tmp_dir` after hardnormly.sh creates it). This is a justified global alongside `_LOG_FILE` and `_DEBUG` — it's infrastructure state, not pipeline state.

3. **shfmt compatibility with nameref**
   - What we know: shfmt formats `local -n ref="$1"` correctly (it's valid bash 4+ syntax).
   - What's unclear: Whether shellcheck flags nameref usage.
   - Recommendation: Verify with `shellcheck -S warning lib/cli.sh` after implementing. If shellcheck warns, add a `# shellcheck disable=SC2034` comment on the nameref line (it sometimes flags namerefs as "unused").

---

## Sources

### Primary (MEDIUM confidence — bash domain has no Context7 equivalent)

- gabrielstaples.com/bash-libraries/ — module structure, `BASH_SOURCE` path resolution, explicit vs auto-load comparison
- wafaicloud.com/blog/crafting-modular-bash-script-libraries/ — include guard pattern, global vs argument passing
- betterdev.blog/minimal-safe-bash-script-template/ — `die()` pattern, trap cleanup, stderr for messaging
- bats-core.readthedocs.io/en/stable/writing-tests.html — BATS `load` for sourced modules, setup patterns

### Secondary (MEDIUM confidence — multiple sources agree)

- sap1ens.com/blog/2017/07/01/bash-scripting-best-practices/ — local variables, naming arguments, BINPATH pattern
- bertvv.github.io/cheat-sheets/Bash.html — general bash best practices (site 403'd but search summary verified)
- citizen428.net/blog/bash-error-handling-with-trap/ — trap ERR limitations
- howtogeek.com/bash-error-handling-patterns-i-use-in-every-script/ — error pattern inventory (set -eE + trap)
- moldstud.com — shell modularity (direct fetch used)

### Tertiary (LOW confidence — search summaries only)

- Google Style Guide mention: "keep scripts under 50 lines" (search result summary, not directly fetched)
- "60% of scripts lacking pipefail encounter silent data loss" (2023 ShellCheck survey mentioned in search result — not verified)

---

## Metadata

**Confidence breakdown:**
- Module structure (lib/ layout, include guards, source pattern): MEDIUM — verified across 3+ credible sources with consistent agreement
- run_cmd pattern: MEDIUM — synthesized from multiple error handling sources; the specific design (stderr to tmp_dir, rich message format) is informed by the project's locked requirements
- Discretion decisions (filter parsing ownership, logging globals, bed/annotate split, no init functions): MEDIUM — based on software engineering principles (SRP, explicit over implicit) applied to bash; no single authoritative "bash SOLID" source exists
- BATS unit test patterns for sourced modules: HIGH — from bats-core official docs

**Research date:** 2026-02-18
**Valid until:** 2026-03-20 (bash best practices are stable; no library version churn)
