# State: hardnormly

## Project Reference

See: .planning/PROJECT.md (updated 2026-02-18 after v0.7.0 milestone)

**Core value:** Reliably normalize and filter VCF files for clean variant calls
**Current focus:** Planning next milestone — run `/gsd:new-milestone`

## Current Position

Phase: Complete (v0.7.0 shipped — 5/5 phases, 19/19 plans)
Plan: Not started
Status: Ready to plan next milestone
Last activity: 2026-02-18 — v0.7.0 milestone complete (archived, tagged v0.7.0)

Progress: [███████████████] 100% (v0.7.0 complete — next milestone not yet defined)

## Accumulated Context

### Key Facts

- v0.7.0 shipped: 130 tests, 44/44 requirements, 5 phases, 19 plans
- Tech stack: Bash, bcftools, bedtools, htslib, BATS, GitHub Actions, GNU Make
- 8 lib/ modules: logging, cli, genome, bed, annotate, normalize, filter, stats
- Main orchestrator: hardnormly.sh ~219 lines
- Issue #13 (genome flag) verified closed
- BATS installed via `tests/setup_bats.sh` (clones pinned versions)
- `.gitattributes eol=lf` for *.sh *.bash *.bats (CRLF prevention)

### Pending Todos

- None

### Open Blockers

- `tests/data/generate_expected.sh` has hardcoded PATH (dev-only, deferred)
- `make help` fails on Windows dev (sh.exe @echo issue) — non-blocking, works on CI
- TFWK-04 real-data tests require `ref/hs37d5.fa` — CI skips with HAVE_FULL_REF guard

## Session Continuity

Last session: 2026-02-18
Stopped at: v0.7.0 milestone archived and tagged. Ready for `/gsd:new-milestone`.
Resume file: None
