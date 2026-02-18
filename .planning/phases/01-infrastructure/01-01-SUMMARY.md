---
phase: "01"
plan: "01"
name: "shellcheck-shfmt-config"
subsystem: "tooling"
tags: ["shellcheck", "shfmt", "editorconfig", "linting", "formatting"]

dependency-graph:
  requires: []
  provides:
    - "ShellCheck baseline: zero warnings on hardnormly.sh"
    - "shfmt formatting baseline: zero diffs on hardnormly.sh"
    - ".shellcheckrc with shell=bash and optional checks"
    - ".editorconfig with tab indentation for shell files"
  affects:
    - "01-02: strict mode and error handling (builds on clean ShellCheck baseline)"
    - "01-03: CI pipeline (shellcheck and shfmt checks must pass)"
    - "01-04: pre-commit hooks (runs shellcheck and shfmt)"

tech-stack:
  added:
    - "shellcheck 0.11.0 (static analysis)"
    - "shfmt v3.12.0 (shell formatting)"
  patterns:
    - "SC2294 temporary disable pattern: documented inline exception with plan-based resolution note"
    - "Separate local/assign for command substitution (SC2155 fix)"
    - "Tab-based indentation for all shell files"

file-tracking:
  key-files:
    created:
      - path: ".shellcheckrc"
        purpose: "ShellCheck configuration: shell=bash, optional checks enabled"
      - path: ".editorconfig"
        purpose: "Editor config: tab indentation for .sh, space for .yml, LF line endings"
    modified:
      - path: "hardnormly.sh"
        changes: "SC2155 fix (separate local/assign in log_msg), SC2294 temporary disable on eval with documented plan-02 exception, shfmt formatting (tabs, -bn, -ci)"

decisions:
  - id: "D-01-01-01"
    decision: "SC2294 (eval) gets single documented temporary disable, not fixed here"
    rationale: "Fixing eval requires structural refactoring of the filter pipeline, which is Plan 02's scope. A premature fix here would muddle two concerns."
    alternatives: "Fix eval now (rejected: scope creep), suppress permanently (rejected: violates zero-disable policy)"
  - id: "D-01-01-02"
    decision: "shfmt flags: -i 0 -bn -ci"
    rationale: "-i 0 forces tabs consistent with .editorconfig; -bn keeps pipelines readable with | at line start; -ci makes case bodies consistent"
    alternatives: "4-space indent (rejected: .editorconfig says tabs); no -bn (rejected: less readable pipelines)"

metrics:
  duration: "~15 minutes"
  completed: "2026-02-18"
  tasks-completed: 2
  tasks-total: 2
  commits: 2
---

# Phase 01 Plan 01: ShellCheck and shfmt Config Summary

**One-liner:** ShellCheck zero-warning baseline with SC2294 documented exception, shfmt tab-formatting via .shellcheckrc and .editorconfig.

## What Was Built

Established the linting and formatting baseline for `hardnormly.sh`:

1. **`.shellcheckrc`** — configures ShellCheck to target bash, enables three optional checks (`require-double-brackets`, `deprecate-which`, `check-set-e-suppressed`), and allows external sources.

2. **`.editorconfig`** — enforces LF line endings, final newlines, trailing whitespace trimming, UTF-8 charset; tab indentation for `.sh` files; 2-space indentation for YAML.

3. **`hardnormly.sh`** — fixed all ShellCheck findings and applied shfmt formatting:
   - **SC2155 fix:** `log_msg()` now declares `local timestamp` separately before the command substitution assignment.
   - **SC2294 temporary disable:** Single `# shellcheck disable=SC2294` on the `eval` line with a code comment explaining this is a Plan 02 concern (eval replacement via array-based pipeline).
   - **shfmt formatting:** Converted 4-space indentation to tabs, applied `-bn` (binary operators at line start) and `-ci` (case body indentation).

## Commits

| Hash | Type | Description |
|------|------|-------------|
| `f6f752c` | chore | Create ShellCheck and editor config files, fix SC2155/SC2294 |
| `64e7815` | style | Apply shfmt formatting to hardnormly.sh |

## Verification Results

| Check | Result |
|-------|--------|
| `shellcheck hardnormly.sh` | PASS (exit 0, no output) |
| `shfmt -d -i 0 -bn -ci hardnormly.sh` | PASS (exit 0, no diffs) |
| `.shellcheckrc` exists with `shell=bash` | PASS |
| `.editorconfig` exists with `indent_style = tab` | PASS |
| No inline disables except documented SC2294 | PASS |

## Decisions Made

### D-01-01-01: SC2294 Temporary Exception

The `eval "$pipeline_cmd"` on line 480 gets a single `# shellcheck disable=SC2294` rather than being fixed here. The eval is used to compose a dynamic bcftools filter pipeline at runtime — replacing it requires an array-based sequential pipeline approach, which is the structural refactoring planned for Plan 02. Fixing it here would mix two concerns and make the Plan 02 diff harder to review.

The disable is accompanied by a code comment explaining the exception and referencing Plan 02:

```bash
# SC2294: eval is used here to compose a dynamic bcftools pipeline with variable filter stages.
# This is a temporary exception — eval will be replaced with an array-based sequential pipeline
# in Plan 02 (eval replacement), which is the proper structural fix.
# shellcheck disable=SC2294
if ! eval "$pipeline_cmd"; then
```

### D-01-01-02: shfmt Flags

Chosen flags: `-i 0 -bn -ci`

- `-i 0`: explicit tabs, consistent with `.editorconfig`
- `-bn`: binary operators (`|`, `&&`, `||`) may start a line — produces more readable pipeline chains
- `-ci`: case body indented — consistent visual structure

## Deviations from Plan

None — plan executed exactly as written.

Task 1 (ShellCheck fixes + config files) was committed in a prior session as `f6f752c`. Task 2 (shfmt formatting) was executed and committed as `64e7815` in this session.

## Next Phase Readiness

**Plan 02 (strict mode + error handling) prerequisites:**
- ShellCheck baseline: SATISFIED — zero warnings on current script
- shfmt baseline: SATISFIED — zero diffs on current script
- SC2294 (eval) is documented and intentionally left for Plan 02
- No blockers

**Known technical debt handed to Plan 02:**
- The `eval "$pipeline_cmd"` pattern triggers SC2294 (temporarily suppressed)
- `set -euo pipefail` not yet enabled (Plan 02 scope)
- No `trap ERR` handler (Plan 02 scope)
