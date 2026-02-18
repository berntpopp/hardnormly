# Phase 4: Refactoring - Context

**Gathered:** 2026-02-18
**Status:** Ready for planning

<domain>
## Phase Boundary

Split `hardnormly.sh` from a monolithic script into focused `lib/` modules with consistent error handling. Replace ad-hoc error patterns with a unified `run_cmd` wrapper. All Phase 3 tests must continue passing — the refactoring is invisible to end-users.

</domain>

<decisions>
## Implementation Decisions

### Module boundaries
- Roadmap lists 8 modules (logging, cli, genome, bed, annotate, normalize, filter, stats) — this is the starting point
- Claude may adapt the module count (merge or split) if research justifies it, with documented reasoning
- bed/annotate split, cross-module dependencies, and global vs. argument passing: Claude researches bash modularization best practices informed by DRY, KISS, SOLID, and modularization principles, then decides

### Error handling pattern
- **Fail fast** — first failure stops the pipeline (keep current behavior)
- **Rich error messages** — include the exact command that failed, its exit code, and captured stderr (upgrade from current terse messages)
- **Retry for network operations only** — UCSC MySQL genome query gets 2-3 retries; all local commands (bcftools, bedtools) fail immediately
- **run_cmd replaces ERR trap** — all external commands go through run_cmd; remove the current err_handler ERR trap entirely

### Filter pipeline redesign
- **Keep the sequential BCF temp file chain** — the array-based filter_stages pattern from Phase 1 stays; no switch to single-pipe architecture
- **fill-tags is internal to filter.sh** — bcftools +fill-tags runs as the first step inside filter.sh, not as a separate orchestrator step
- Filter parsing ownership (cli.sh vs filter.sh): Claude researches separation of concerns best practices for our bash stack, then decides
- PASS filtering + output format ownership: Claude decides based on module cohesion

### Source/load convention
- **Include guards required** — every module uses `[[ -n ${_MODULE_LOADED:-} ]] && return; _MODULE_LOADED=1` pattern
- Init functions vs. ready-on-source: Claude decides based on whether modules actually need initialization
- Module locating strategy (explicit source lines vs. loader function): Claude decides based on maintainability
- Dependency declaration in modules (self-checking vs. orchestrator responsibility): Claude decides

### Claude's Discretion
- bed.sh vs. annotate.sh boundary (header file creation ownership)
- Cross-module call pattern (direct calls, orchestrator mediation, or logging-is-special hybrid)
- Global variables vs. function arguments vs. hybrid approach
- Filter parsing ownership (cli.sh parses vs. filter.sh self-contained)
- PASS filter + output format step ownership
- Init function pattern (ready-on-source vs. explicit init)
- Module loading strategy (explicit source lines vs. auto-loader)
- Module dependency declaration approach

All "Claude's Discretion" items should be decided through research into bash modularization best practices, with findings documented in research output.

</decisions>

<specifics>
## Specific Ideas

- User wants decisions informed by senior developer best practices — research bash modularization patterns, not just "what works"
- Apply DRY, KISS, SOLID principles adapted to bash context
- Modern and maintainable setup is the goal — not minimal-effort migration
- eval-based pipeline is already gone (replaced in Phase 1 with array-based approach) — roadmap's "replace eval" criterion is already met

</specifics>

<deferred>
## Deferred Ideas

None — discussion stayed within phase scope

</deferred>

---

*Phase: 04-refactoring*
*Context gathered: 2026-02-18*
