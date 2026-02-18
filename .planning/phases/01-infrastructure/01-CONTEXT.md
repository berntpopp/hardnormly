# Phase 1: Infrastructure - Context

**Gathered:** 2026-02-18
**Status:** Ready for planning

<domain>
## Phase Boundary

Enforce consistent code quality automatically — locally and in CI. This covers linting (ShellCheck), formatting (shfmt), strict shell mode (set -euo pipefail), GitHub Actions CI, and pre-commit hooks. No test data, no BATS tests, no refactoring — those are later phases.

</domain>

<decisions>
## Implementation Decisions

### ShellCheck & shfmt config
- `.shellcheckrc` in repo root — central config, checked into git
- Zero warnings policy: all ShellCheck warnings must be resolved, no inline `# shellcheck disable` directives
- shfmt uses tabs (default style)
- shfmt enforces full formatting: indentation + binary operator position + case body indent

### Strict mode strategy
- `set -euo pipefail` at the top of hardnormly.sh — global, not section-based
- Exception handling pattern: researcher should investigate best practices for handling legitimate non-zero exits under strict mode and pick the best pattern per case
- Global trap on both ERR and EXIT: log failing command/line number and clean up temp dirs
- User prefers replacing eval in Phase 1 rather than deferring to Phase 4; researcher should assess feasibility — if eval replacement requires the lib/ modular structure from Phase 4, keep eval working under strict mode for now and replace in Phase 4 as originally planned

### CI workflow design
- Triggers: push to main + pull requests targeting main
- OS: Ubuntu only (ubuntu-latest)
- Phase 1 scope: lint only (ShellCheck + shfmt). No bcftools/bedtools install, no smoke tests — those come in Phase 3
- Structure: single `.github/workflows/ci.yml` with multiple jobs (lint job, test job placeholder)

### Pre-commit hook scope
- Checks both ShellCheck and shfmt on commit
- Simple git hook script (plain bash in `.githooks/`), not the pre-commit.com framework
- Checks staged `.sh` files only, not all files in repo
- Developer activation: `scripts/setup-hooks.sh` that runs `git config core.hooksPath .githooks`

### Claude's Discretion
- Exact `.shellcheckrc` options beyond shell=bash
- shfmt flags beyond indentation/binary-ops/case (e.g., function opening brace, redirect operators)
- ERR trap implementation details (how to capture line number, formatting of error output)
- CI job naming and step organization
- Whether to add a Makefile or just document commands

</decisions>

<specifics>
## Specific Ideas

- Exception handling under strict mode should be researched — user wants best practices applied, not a one-size-fits-all pattern
- eval replacement preference is strong: user wants it in Phase 1 if feasible, but roadmap has it in Phase 4. Researcher should assess whether it's practical without the lib/ modular structure

</specifics>

<deferred>
## Deferred Ideas

None — discussion stayed within phase scope

</deferred>

---

*Phase: 01-infrastructure*
*Context gathered: 2026-02-18*
