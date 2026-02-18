# Phase 1: Infrastructure - Research

**Researched:** 2026-02-18
**Domain:** Bash shell quality tooling (ShellCheck, shfmt, strict mode, GitHub Actions CI, git hooks)
**Confidence:** HIGH

## Summary

Phase 1 establishes the code quality foundation for `hardnormly.sh`. The work covers four distinct areas: (1) ShellCheck static analysis with a `.shellcheckrc` config and zero-warning policy, (2) shfmt formatting enforcement, (3) `set -Eeuo pipefail` strict mode with targeted exception handling and ERR/EXIT traps, and (4) GitHub Actions CI plus pre-commit hooks.

The standard approach for each area is well-established and verified through official documentation. ShellCheck is pre-installed on GitHub ubuntu-latest runners; shfmt requires an explicit install step. The `eval` pipeline in hardnormly.sh is replaceable in Phase 1 without the lib/ modular structure — the pipeline stages map cleanly to an array-of-arrays pattern feeding through `bcftools` stdin. This is assessed as **feasible in Phase 1** (see eval assessment section).

**Primary recommendation:** Use native `shellcheck` and `shfmt` CLI commands directly in CI (not composite marketplace actions) for transparency. Use `set -Eeuo pipefail` with `set -E` for trap inheritance. Replace `eval` with an array-based pipeline that threads stdin through sequential bcftools calls using `-`.

## Standard Stack

### Core

| Tool | Version | Purpose | Why Standard |
|------|---------|---------|--------------|
| ShellCheck | 0.10.x (pre-installed on ubuntu-latest) | Static analysis for shell scripts | De facto standard; pre-installed on GitHub runners; SC codes are widely documented |
| shfmt | 3.8.x (not pre-installed; needs install) | Shell script formatter | Only maintained formatter for bash; mvdan/sh project |
| GitHub Actions | N/A | CI platform | Project already on GitHub; no setup needed |
| bash 5.x | 5.x (ubuntu-latest default) | Script runtime | Already in use; `set -Eeuo pipefail` requires bash |

### Supporting

| Tool/Pattern | Version | Purpose | When to Use |
|-------------|---------|---------|-------------|
| mfinelli/setup-shfmt | v4 | Install shfmt in CI | When you want a specific shfmt version pinned |
| ludeeus/action-shellcheck | v2 | Run ShellCheck with extra options | When you want severity filtering or SARIF output |
| apt-get shfmt | system version | Install shfmt in CI | Acceptable; ubuntu-latest apt provides shfmt 3.x |

### Alternatives Considered

| Instead of | Could Use | Tradeoff |
|------------|-----------|----------|
| Native shellcheck in CI | ludeeus/action-shellcheck | Marketplace action adds complexity/version drift; native is simpler and equally capable |
| Native shfmt in CI | reviewdog/action-shfmt | reviewdog adds PR annotation features but is overkill for lint-only phase |
| Plain git hooks | pre-commit.com framework | Framework requires Python install and is heavier; plain bash hooks are zero-dependency and match the project's bash-first character |

**Installation (CI):**
```bash
# ShellCheck - already pre-installed on ubuntu-latest
shellcheck --version

# shfmt - must be installed (not pre-installed)
sudo apt-get update && sudo apt-get install -y shfmt
# OR use action: mfinelli/setup-shfmt@v4
```

**Installation (local dev):**
```bash
# Debian/Ubuntu
sudo apt-get install shellcheck shfmt

# macOS
brew install shellcheck shfmt

# Conda (already in hardnormly_environment.yml if added)
conda install -c conda-forge shellcheck shfmt
```

## Architecture Patterns

### Recommended Project Structure

```
hardnormly/
├── .shellcheckrc              # ShellCheck config (repo root, checked into git)
├── .github/
│   └── workflows/
│       └── ci.yml             # CI workflow (lint job + test placeholder job)
├── .githooks/
│   └── pre-commit             # Git hook script (plain bash)
├── scripts/
│   └── setup-hooks.sh         # Developer activation script
└── hardnormly.sh              # Main script (modified for strict mode)
```

### Pattern 1: `.shellcheckrc` Configuration

**What:** Central ShellCheck config at repo root, read automatically during CI and local runs.
**When to use:** Always — ensures consistent checking across all environments.

```ini
# Source: https://github.com/koalaman/shellcheck/blob/master/shellcheck.1.md
# .shellcheckrc

# Target shell: bash (ensures bash-specific checks)
shell=bash

# Enable optional checks
# require-double-brackets: warn on [ ] usage in bash (prefer [[ ]])
enable=require-double-brackets
# deprecate-which: warn on `which`, prefer `command -v`
enable=deprecate-which
# check-set-e-suppressed: warn when set -e is suppressed during function calls
enable=check-set-e-suppressed

# Allow sourcing files that shellcheck can't find (common in bash scripts)
external-sources=true
```

**Notes on optional checks:**
- `require-double-brackets`: hardnormly.sh already uses `[[ ]]` consistently — enabling this confirms it and catches any regressions
- `deprecate-which`: no `which` calls currently visible, low impact
- `check-set-e-suppressed`: valuable after adding `set -e`, helps catch cases where strict mode is accidentally bypassed
- Do NOT enable `require-variable-braces` — this would be extremely noisy in the current script (hundreds of `$var` → `${var}` changes) and the user decision specifies zero warnings, not zero-style-violations

### Pattern 2: shfmt Flags

**What:** shfmt flags that match the decisions: tabs, binary ops on new line, case indent.
**When to use:** Apply formatting once to hardnormly.sh, then enforce with `-d` in CI.

```bash
# Source: https://github.com/mvdan/sh/blob/v3.8.0/cmd/shfmt/shfmt.1.scd

# Format (write in place):
shfmt -w -i 0 -bn -ci hardnormly.sh

# Check (exit 1 if diff):
shfmt -d -i 0 -bn -ci hardnormly.sh

# Flags used:
#   -i 0   = tabs (default, explicit)
#   -bn    = binary ops (&&, ||, |) may start a line
#   -ci    = case body indented
#
# Flags NOT used (Claude's discretion decisions):
#   -sr    = space after redirects (off — hardnormly.sh uses >file, <file style)
#   -fn    = function brace on new line (off — standard bash style is same line)
#   -kp    = keep column alignment (off — removes manual padding, too disruptive)
```

**Recommendation on `-sr`:** Do not add `-sr`. The redirect style (`> "$file"` vs `>"$file"`) is already consistent in the script and `-sr` would reformat all redirects, creating a large noisy diff. Keep the current style.

**Recommendation on `-fn`:** Do not add `-fn`. Functions with opening braces on the same line is the POSIX/Google style. `-fn` (brace on new line) is a minority preference. Current script uses same-line braces.

### Pattern 3: `set -Eeuo pipefail` with ERR/EXIT Traps

**What:** Strict mode with `-E` (errtrace) ensures ERR trap fires in functions and subshells. The `-E` flag is the critical addition vs the common `set -euo pipefail`.

```bash
# Source: https://citizen428.net/blog/bash-error-handling-with-trap/
# Place at top of hardnormly.sh, after the shebang and version declaration:

set -Eeuo pipefail

# ERR trap: fires on any command that exits non-zero (with -E, fires inside functions too)
# EXIT trap: fires on script exit for cleanup
# Both reference the SAME handler — EXIT handles cleanup; ERR logs the error first.

err_handler() {
    local exit_code=$?
    local line_number=$1
    local command="${BASH_COMMAND}"
    log_msg "ERROR: Command '${command}' failed with exit code ${exit_code} at line ${line_number}"
    # Note: cleanup_tmp_dir is called by EXIT trap below; don't duplicate here
}

cleanup_handler() {
    # Always clean up temp dir on exit (replaces existing EXIT trap)
    if [[ "${cleanup:-true}" == "true" ]] && [[ -d "${tmp_dir:-}" ]]; then
        rm -rf "$tmp_dir"
    fi
}

trap 'err_handler ${LINENO}' ERR
trap cleanup_handler EXIT
```

**Key variable expansion subtlety:** `${LINENO}` must be expanded in the trap string (single-quoted trap won't work), OR pass it as a parameter. Use the parameter-passing pattern above: `trap 'err_handler ${LINENO}' ERR` — here `${LINENO}` expands at trap-fire time to the correct line number.

**Important:** The existing `trap '[[ $cleanup == true ]] && cleanup_tmp_dir' EXIT` at line 28 must be replaced with the new `cleanup_handler` function. The pattern changes from inline to function-based.

### Pattern 4: Exception Handling Under Strict Mode

**What:** Patterns for commands that legitimately return non-zero under `set -e`.

Under `set -Eeuo pipefail`, any non-zero exit terminates the script. hardnormly.sh has several places where this applies:

**Pattern A — if-condition (automatic bypass):**
```bash
# Bash does NOT trigger set -e for commands used as if-conditions
# These are already correct and need no change:
if grep -q "Warning" "$norm_output"; then   # grep returns 1 on no match — safe here
    log_msg "bcftools norm warnings: ..."
fi
```

**Pattern B — `|| true` for expected failures:**
```bash
# Use when a command may return non-zero and that's acceptable
mysql --user=genome ... 2>/dev/null || {
    log_msg "Error: Failed to create genome file"
    exit 1
}
# But if we WANT to handle the failure ourselves, keep the ||; don't use || true
```

**Pattern C — explicit exit code capture (existing pattern in script):**
```bash
# ALREADY IN SCRIPT — this pattern is correct under set -e because the assignment
# of exit code happens in the same "statement" and set -e doesn't trigger on
# the assignment itself:
norm_exit_code=$?
# However: with set -e, the command before $? capture must NOT have failed
# and triggered set -e first. Better to use:
bcftools norm ... || { log_msg "Error"; exit 1; }
```

**Pattern D — local variable declaration (important gotcha):**
```bash
# WRONG under set -e — 'local' always exits 0, masking the inner command's status:
local result=$(some_command)   # exit code of some_command is LOST

# CORRECT:
local result
result=$(some_command)         # now set -e triggers if some_command fails
```

**Pattern E — `if [[ $? -ne 0 ]]` after command (existing pattern — PROBLEM):**
```bash
# WRONG with set -e — the command at line 401 (eval "$pipeline_cmd") will
# cause set -e to exit before line 404's $? check ever runs:
eval "$pipeline_cmd"
if [[ $? -ne 0 ]]; then       # unreachable!
    log_msg "Error: Failed to filter the VCF."
    exit 1
fi

# CORRECT: use || directly
eval "$pipeline_cmd" || { log_msg "Error: Failed to filter the VCF."; exit 1; }
# Or let ERR trap handle it (no explicit check needed):
eval "$pipeline_cmd"
```

**Specific cases to handle in hardnormly.sh:**
1. Line 190: `mysql ... | grep -v ... | sed ...` — grep/sed return codes may be masked; consider restructuring
2. Line 316: `grep -q "Warning" "$norm_output"` — already in if-condition, safe
3. Line 321: `grep -q "Lines" "$norm_stdout"` — already in if-condition, safe
4. Lines 275, 287, 308, 414, 434: `if [[ $? -ne 0 ]]` pattern — must be replaced with `|| { }` pattern
5. Line 360: `filter_expr=$(echo "$filter_expr" | tr -d '\r\n')` — command substitution inside variable assignment; safe with set -e if subshell exits non-zero? Actually safe — set -e fires if the entire pipeline fails.

### Pattern 5: eval Replacement Assessment

**Feasibility verdict: FEASIBLE in Phase 1.** The eval pipeline is replaceable without the lib/ modular structure, but it requires a design change.

**Current design (eval):**
```bash
# Builds a pipeline string and evals it:
pipeline_cmd="bcftools view $normalized_vcf | bcftools +fill-tags"
# ... conditionally appends: | bcftools filter -s NAME -e 'EXPR'
# ... more conditional appends
eval "$pipeline_cmd"
```

**Problem with eval:** ShellCheck SC2294 will fire. Also, the filter expressions contain single-quoted strings that must be preserved — when building `filter_cmd="bcftools filter -m+ -sNAME -$action '$expr'"`, the single quotes inside the double-quoted string are just literal characters, which eval interprets correctly. BUT: ShellCheck will flag SC2294 AND the security implications.

**Replacement pattern — sequential stdin threading:**

Since `bcftools` accepts `-` as stdin, the pipeline can be replaced by building an array of filter argument-sets and iterating through them, piping via stdin:

```bash
# Source: principle from https://mywiki.wooledge.org/BashFAQ/048
# Each bcftools command reads from the previous step's output

# Step 1: initial pipeline (always runs)
bcftools view "$normalized_vcf" | bcftools +fill-tags > "$tmp_dir/step_current.bcf"

# Build array of filter argument arrays
# Each entry: "filter_name|action|expression" (encoded with separator to avoid subshell issues)
filter_stages=()

if [[ -f "$tmp_dir/merged_include_regions.bed.gz" ]]; then
    filter_stages+=("NOT_IN_INCLUDE_REGION|e|INFO/INCLUDE_REGION!=1")
fi
if [[ -f "$tmp_dir/merged_exclude_regions.bed.gz" ]]; then
    filter_stages+=("IN_EXCLUDE_REGION|e|INFO/EXCLUDE_REGION=1")
fi
for filter in "${filters[@]}"; do
    IFS=" " read -r filter_name filter_action filter_expr <<< "$filter"
    filter_stages+=("${filter_name}|${filter_action}|${filter_expr}")
done
if [[ -n "$filters_file" ]]; then
    while IFS=" " read -r filter_name filter_action filter_expr; do
        filter_expr=$(tr -d '\r\n' <<< "$filter_expr")
        filter_stages+=("${filter_name}|${filter_action}|${filter_expr}")
    done < "$filters_file"
fi

# Apply each filter stage sequentially via temp files
for stage in "${filter_stages[@]}"; do
    IFS="|" read -r stage_name stage_action stage_expr <<< "$stage"
    bcftools filter -m+ -s "$stage_name" -"$stage_action" "$stage_expr" \
        "$tmp_dir/step_current.bcf" -Ob -o "$tmp_dir/step_next.bcf"
    mv "$tmp_dir/step_next.bcf" "$tmp_dir/step_current.bcf"
done

# Final output
if $only_pass; then
    bcftools view -f PASS "$tmp_dir/step_current.bcf" -O"$output_type" -o "$output_vcf"
else
    bcftools view "$tmp_dir/step_current.bcf" -O"$output_type" -o "$output_vcf"
fi
```

**Alternative replacement pattern — inline pipe chaining via stdin:**
bcftools supports stdin via `-`, so each call in a loop can read from a piped process. But in bash, you can't easily build a pipe chain iteratively without eval. The temp-file approach above is the cleanest eval-free replacement.

**Trade-off with temp-file approach:** Disk I/O for each filter stage (potentially 10+ filters). For VCF files, this is usually acceptable since bcftools can use BCF format (binary, fast). Use `-Ob` (BCF output) for intermediate files and `-Oz` (bgzipped VCF) only for final output.

**Conclusion:** Replace eval with the sequential temp-file approach using BCF format for intermediates. This is self-contained, does NOT require lib/ modularization, and eliminates ShellCheck SC2294 plus the security concern. Estimated complexity: moderate — 30-40 lines of restructured code.

### Pattern 6: GitHub Actions CI Structure

**What:** Single `ci.yml` with multiple jobs. Phase 1 scope is lint only.

```yaml
# Source: https://github.com/koalaman/shellcheck/wiki/GitHub-Actions
# .github/workflows/ci.yml
name: CI

on:
  push:
    branches: [main]
  pull_request:
    branches: [main]

jobs:
  lint:
    name: Lint
    runs-on: ubuntu-latest
    steps:
      - name: Checkout
        uses: actions/checkout@v4

      - name: Install shfmt
        run: sudo apt-get update && sudo apt-get install -y shfmt

      - name: Run ShellCheck
        run: shellcheck hardnormly.sh scripts/setup-hooks.sh .githooks/pre-commit

      - name: Check shfmt formatting
        run: shfmt -d -i 0 -bn -ci hardnormly.sh scripts/setup-hooks.sh .githooks/pre-commit

  test:
    name: Test
    runs-on: ubuntu-latest
    needs: []  # Independent of lint; runs in parallel
    steps:
      - name: Checkout
        uses: actions/checkout@v4

      - name: Placeholder (BATS tests — Phase 3)
        run: echo "No tests yet — Phase 3 adds BATS"
```

**Notes:**
- ShellCheck is pre-installed on ubuntu-latest; no install step needed
- shfmt must be installed via apt-get (or use `mfinelli/setup-shfmt@v4` for version pinning)
- `actions/checkout@v4` is current as of 2026-02-18
- The test job runs in parallel with lint (no `needs: [lint]`) — fail fast on either
- Pinning to specific action versions (v4) vs `@master` avoids supply-chain risk

### Pattern 7: Pre-commit Hook

**What:** Plain bash script in `.githooks/pre-commit` that checks staged `.sh` files only.

```bash
#!/bin/bash
# .githooks/pre-commit
# Checks staged shell scripts with ShellCheck and shfmt before commit.

set -euo pipefail

# Get list of staged .sh files (added, copied, or modified)
staged_sh_files=$(git diff --cached --name-only --diff-filter=ACM | grep '\.sh$' || true)

if [[ -z "$staged_sh_files" ]]; then
    exit 0  # No shell files staged; nothing to check
fi

echo "Pre-commit: checking staged shell files..."

failed=0

while IFS= read -r file; do
    echo "  ShellCheck: $file"
    if ! shellcheck "$file"; then
        failed=1
    fi

    echo "  shfmt:      $file"
    if ! shfmt -d -i 0 -bn -ci "$file"; then
        failed=1
    fi
done <<< "$staged_sh_files"

if [[ $failed -ne 0 ]]; then
    echo ""
    echo "Pre-commit hook failed. Fix the issues above before committing."
    echo "Run: shellcheck <file>  and  shfmt -w -i 0 -bn -ci <file>"
    exit 1
fi

echo "Pre-commit: all checks passed."
```

**setup-hooks.sh:**
```bash
#!/bin/bash
# scripts/setup-hooks.sh
# Configures git to use .githooks/ directory for this repository.

set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd "${script_dir}/.." && pwd)"

git -C "$repo_root" config core.hooksPath .githooks
echo "Git hooks configured: .githooks/"
echo "Run 'git config core.hooksPath .githooks' manually if this fails."
```

### Anti-Patterns to Avoid

- **`# shellcheck disable=SC2086`:** User decision mandates zero inline disables. Fix the actual issue.
- **`set -e` without `-E`:** Traps won't fire inside functions. Always use `set -Eeuo pipefail`.
- **`local result=$(cmd)` pattern:** Masks exit code. Always declare `local` separately from assignment.
- **`if [[ $? -ne 0 ]]` after a command:** The command's failure triggers `set -e` before the check runs. Use `cmd || { handle error; }` pattern instead.
- **`eval` on strings that contain array variables:** ShellCheck SC2294. Don't build command strings as text — use arrays or the temp-file pipeline approach.
- **Using `@master` for actions:** Supply-chain risk. Pin to version tags (`@v4`).
- **Checking all `.sh` files on commit (not just staged):** Slow and catches irrelevant files. Use `git diff --cached --name-only`.
- **`grep` exit code under `set -e`:** `grep` returns 1 when no matches found. Use `grep ... || true` when the "no match" case is acceptable, or put in an if-condition.

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| Shell static analysis | Custom lint rules | ShellCheck | Covers 500+ checks, actively maintained, GitHub-aware |
| Shell formatting | Custom awk/sed formatter | shfmt | Only maintained formatter; handles all bash constructs correctly |
| Pipeline composition | Custom eval-based string building | Array-based sequential pipeline with temp files | eval is a ShellCheck violation and security risk; temp-file approach is clean |
| CI shell script validation | Custom docker-based checks | Native shellcheck + shfmt in GitHub Actions | Pre-installed / easily installed; no custom images needed |
| Git hook framework | Pre-commit.com setup | Plain bash `.githooks/pre-commit` | Zero dependencies; fits the project's bash-first character |

**Key insight:** The eval pattern looks like it needs a custom solution, but bcftools' stdin support (`-` input) enables a clean sequential-application pattern using temp BCF files.

## Common Pitfalls

### Pitfall 1: `set -e` Without `-E` (errtrace)

**What goes wrong:** ERR trap fires in the main script body but NOT inside functions or command substitutions. Errors inside `normalize_bed()`, `log_msg()`, etc. silently fail or exit without the trap running.
**Why it happens:** `set -e` only applies errtrace to the current shell by default.
**How to avoid:** Always use `set -Eeuo pipefail` (capital E). Verify with a function that fails: `foo() { false; }; foo` — should fire the ERR trap.
**Warning signs:** ERR trap message doesn't appear when a function fails.

### Pitfall 2: `local` Masking Exit Codes

**What goes wrong:** `local result=$(failing_cmd)` always returns 0 because `local` is a builtin that always succeeds. The failing command's exit code is lost.
**Why it happens:** In bash, `local var=$(cmd)` evaluates to the exit code of `local`, not `cmd`.
**How to avoid:** Declare local separately: `local result; result=$(cmd)`.
**Warning signs:** ShellCheck SC2155 warns about this exact pattern.

### Pitfall 3: `if [[ $? -ne 0 ]]` After Commands with `set -e`

**What goes wrong:** The command fails, `set -e` immediately terminates the script, and the `if [[ $? -ne 0 ]]` check is unreachable.
**Why it happens:** hardnormly.sh currently uses this pattern in multiple places (lines 275, 287, 308, 414, 434).
**How to avoid:** Replace with `cmd || { log_msg "Error..."; exit 1; }` or let the ERR trap handle it.
**Warning signs:** grep for `if \[\[ \$? -ne 0 \]\]` in the script.

### Pitfall 4: `grep` Exit Code 1 Under `set -e`

**What goes wrong:** `grep` returns 1 when no lines match. Under `set -e`, this terminates the script.
**Why it happens:** In hardnormly.sh, `grep -q "Warning" "$norm_output"` (line 316) and `grep -q "Lines" "$norm_stdout"` (line 321) are already in if-conditions (safe). But line 190 has `| grep -v "^chrom" | sed ...` in a mysql pipeline — grep returns 1 if no non-chrom lines exist (unlikely but possible).
**How to avoid:** Put grep in if-conditions, or append `|| true` when "no match" is acceptable.
**Warning signs:** Script exits unexpectedly when processing files with no matches.

### Pitfall 5: shfmt Not Pre-installed on GitHub Runners

**What goes wrong:** CI fails with "shfmt: command not found" even though ShellCheck works fine.
**Why it happens:** ShellCheck IS pre-installed on ubuntu-latest; shfmt is NOT.
**How to avoid:** Add explicit `sudo apt-get install -y shfmt` step before the shfmt check step.
**Warning signs:** Checking the GitHub-hosted runners software list — shfmt is not in the default tool cache.

### Pitfall 6: eval Pipeline Expressions with Single Quotes

**What goes wrong:** The current `pipeline_cmd` string contains single quotes inside double quotes: `"bcftools filter -m+ -s$filter_name -$filter_action '$filter_expr'"`. When replacing eval with arrays, directly adding `"'$expr'"` to an array element creates an element with literal single quotes, which bcftools receives correctly. But developers may introduce quoting bugs.
**Why it happens:** The quoting model changes from "string eval" to "array element is the exact argument".
**How to avoid:** In the array approach, filter expressions go directly as array elements — NO surrounding quotes needed: `bcftools_args+=("-e" "$stage_expr")`. The bcftools `-e` flag takes the expression as its argument; the shell handles quoting when executing `"${cmd[@]}"`.
**Warning signs:** ShellCheck SC2206 (word splitting risk when building arrays from strings).

### Pitfall 7: Pre-commit Hook Not Running After Setup

**What goes wrong:** Developer runs `setup-hooks.sh` but hooks don't fire on commit.
**Why it happens:** `git config core.hooksPath` is repo-local — it only works for the current repo. Hooks must be executable (`chmod +x .githooks/pre-commit`).
**How to avoid:** `setup-hooks.sh` should also run `chmod +x .githooks/*` or the files should be committed with executable bit set.
**Warning signs:** `git config core.hooksPath` not set, or hook file not executable.

## Code Examples

Verified patterns from official sources:

### ShellCheck .shellcheckrc (Verified)

```ini
# Source: https://github.com/koalaman/shellcheck/blob/master/shellcheck.1.md
shell=bash
enable=require-double-brackets
enable=deprecate-which
enable=check-set-e-suppressed
external-sources=true
```

### shfmt Format + Check Commands (Verified)

```bash
# Source: https://github.com/mvdan/sh/blob/v3.8.0/cmd/shfmt/shfmt.1.scd

# Apply formatting:
shfmt -w -i 0 -bn -ci hardnormly.sh

# Check (CI mode — exits 1 with diff if different):
shfmt -d -i 0 -bn -ci hardnormly.sh
```

### ERR Trap with Line Number (Verified)

```bash
# Source: https://phaq.phunsites.net/2010/11/22/trap-errors-exit-codes-and-line-numbers-within-a-bash-script/
# Combined with: https://citizen428.net/blog/bash-error-handling-with-trap/

set -Eeuo pipefail

err_handler() {
    local exit_code=$?
    local line_number=$1
    local failed_command="${BASH_COMMAND}"
    log_msg "ERROR: '${failed_command}' failed (exit ${exit_code}) at line ${line_number}"
}

trap 'err_handler ${LINENO}' ERR
trap cleanup_handler EXIT
```

### Staged Files Pre-commit Check (Verified Pattern)

```bash
# Gets only staged .sh files; grep returns 0 (no exit via set -e) because || true handles no matches
staged_sh_files=$(git diff --cached --name-only --diff-filter=ACM | grep '\.sh$' || true)
```

### bcftools Sequential Filter (No eval) (Derived Pattern)

```bash
# Source: bcftools docs (stdin via -), ShellCheck SC2294 guidance
# Build filter stages as encoded strings to avoid subshell issues

filter_stages=()
# ... populate filter_stages array as shown in Pattern 5 above ...

# Apply each filter sequentially via temp BCF files
current_input="$tmp_dir/step_current.bcf"
bcftools view "$normalized_vcf" | bcftools +fill-tags -Ob -o "$current_input"

for stage in "${filter_stages[@]}"; do
    IFS="|" read -r stage_name stage_action stage_expr <<< "$stage"
    bcftools filter -m+ -s "$stage_name" "-${stage_action}" "$stage_expr" \
        "$current_input" -Ob -o "$tmp_dir/step_next.bcf"
    mv "$tmp_dir/step_next.bcf" "$current_input"
done
```

### GitHub Actions CI Workflow (Verified Structure)

```yaml
# Source: https://github.com/koalaman/shellcheck/wiki/GitHub-Actions
name: CI

on:
  push:
    branches: [main]
  pull_request:
    branches: [main]

jobs:
  lint:
    name: Lint
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4

      - name: Install shfmt
        run: sudo apt-get update && sudo apt-get install -y shfmt

      - name: ShellCheck
        run: shellcheck hardnormly.sh scripts/setup-hooks.sh .githooks/pre-commit

      - name: shfmt format check
        run: shfmt -d -i 0 -bn -ci hardnormly.sh scripts/setup-hooks.sh .githooks/pre-commit

  test:
    name: Test
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4
      - name: Placeholder
        run: echo "BATS tests added in Phase 3"
```

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| `#!/bin/bash` with no strict mode | `set -Eeuo pipefail` + `-E` for errtrace | Standard practice since ~2015, `-E` underused | Catches more errors in functions |
| `set -e` without `-E` | `set -E` (errtrace) | `-E` is the correct companion for ERR traps in functions | Traps now fire inside functions |
| `local var=$(cmd)` | `local var; var=$(cmd)` | Always correct; SC2155 added to ShellCheck ~2019 | Exit code preserved |
| `actions/checkout@v2` | `actions/checkout@v4` | 2023 | Required for newer runner features |
| `eval` for dynamic pipelines | Array + sequential stdin threading | Current recommendation | Eliminates SC2294, security risk |

**Deprecated/outdated:**
- `if [[ $? -ne 0 ]]` after commands: Replace with `|| { }` pattern under strict mode
- `trap '[[ $cleanup == true ]] && cleanup_tmp_dir' EXIT` (inline trap): Replace with function-based trap
- `grep` exit code captured after pipe: Check inside if-condition or use `|| true`

## Open Questions

1. **bcftools `-Ob` intermediate performance impact**
   - What we know: bcftools BCF format is binary and fast; intermediate temp files are in the system temp dir
   - What's unclear: For very large VCF files (multi-GB), the extra disk I/O of the sequential approach vs in-memory piping could matter
   - Recommendation: Proceed with temp-file approach; performance concern is out of scope for Phase 1 (infrastructure phase, not optimization phase)

2. **shfmt version on ubuntu-latest apt**
   - What we know: `sudo apt-get install shfmt` installs from ubuntu-latest apt repos; version is typically 3.x
   - What's unclear: Whether the apt repo version matches the latest shfmt 3.8.x
   - Recommendation: Use apt install for simplicity; if version pinning is needed, use `mfinelli/setup-shfmt@v4` with `shfmt-version: 3.8.0`

3. **ShellCheck findings in hardnormly.sh before fixes**
   - What we know: The script currently has no `set -euo pipefail`, uses `eval`, has `if [[ $? -ne 0 ]]` patterns, and `local timestamp=$(date ...)` (SC2155)
   - What's unclear: Total count of ShellCheck warnings — needs a local `shellcheck hardnormly.sh` run to enumerate
   - Recommendation: Planner should include a task step to run `shellcheck hardnormly.sh` first and enumerate findings before scheduling fix tasks

4. **Whether to add a Makefile (Claude's Discretion)**
   - What we know: A Makefile would provide `make lint`, `make format`, `make setup-hooks` targets for developer convenience
   - Recommendation: Add a minimal Makefile with `lint`, `format`, and `setup-hooks` targets. The project already has shell scripts for everything; Makefile just provides discovery and standard entrypoints. Keep it simple (5-10 targets max).

## Sources

### Primary (HIGH confidence)
- ShellCheck official docs: https://github.com/koalaman/shellcheck/blob/master/shellcheck.1.md — `.shellcheckrc` options, optional checks
- ShellCheck optional checks wiki: https://github.com/koalaman/shellcheck/wiki/Optional — all optional check names with descriptions
- shfmt man page: https://github.com/mvdan/sh/blob/v3.8.0/cmd/shfmt/shfmt.1.scd — all formatting flags
- ShellCheck GitHub Actions wiki: https://github.com/koalaman/shellcheck/wiki/GitHub-Actions — pre-installation status, usage patterns
- SC2294 wiki: https://www.shellcheck.net/wiki/SC2294 — eval alternatives

### Secondary (MEDIUM confidence)
- ERR trap patterns: https://phaq.phunsites.net/2010/11/22/trap-errors-exit-codes-and-line-numbers-within-a-bash-script/ — line number capture pattern (older but verified correct for current bash)
- Bash error handling: https://citizen428.net/blog/bash-error-handling-with-trap/ — `set -Eeuo pipefail` + trap setup; `set -E` for functions
- HowToGeek bash error handling: https://www.howtogeek.com/bash-error-handling-patterns-i-use-in-every-script/ — `|| true` patterns, `set -eE`
- Translucent Computing pipefail: https://translucentcomputing.com/blog/unofficial-bash-strict-mode-pipefail/ — `|| true`, if-condition bypass patterns

### Tertiary (LOW confidence)
- mfinelli/setup-shfmt action: https://github.com/marketplace/actions/setup-shfmt — version and YAML usage (verified from marketplace page)
- action-sh-checker: https://github.com/luizm/action-sh-checker — v0.10.0 (Jan 2026), combined ShellCheck+shfmt action (considered but not recommended for this project)

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH — ShellCheck/shfmt are well-documented; installation status verified
- Architecture patterns: HIGH — `.shellcheckrc`, shfmt flags, CI YAML, hook patterns all from official/primary sources
- Strict mode + traps: HIGH — patterns verified from multiple sources, consistent advice
- eval replacement: MEDIUM — the temp-file approach is derived from principles (bcftools stdin support + array pattern), not directly verified with hardnormly-specific code; feasibility assessment is confident but implementation details need validation during coding
- Pitfalls: HIGH — specific pitfalls identified from actual code analysis of hardnormly.sh (lines referenced)

**Research date:** 2026-02-18
**Valid until:** 2026-04-18 (60 days — ShellCheck/shfmt are stable; CI patterns change slowly)
