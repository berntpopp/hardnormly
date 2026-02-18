# Documentation Plan: hardnormly

**Goal:** Create a `docs/` folder with detailed documentation + a concise README that links into it, rendered as a GitHub Pages site via MkDocs Material.

## Audit of Current State

The existing README.md (183 lines) is **stale** — reflects pre-refactoring v0.5:

| Gap | Details |
|-----|---------|
| Missing subcommands | `run-pipeline`, `generate-inclusion-bed`, `generate-exclusion-bed` |
| Missing 12+ flags | `--caller`, `--strip-annotations`, `--auto-index`, `--only-pass`, `--no-cleanup`, `--debug`, `--log-file`, `--tmp-dir`, `--slop`, `--generate-stats`, `--plot-stats`, `--plot-output-dir` |
| Wrong filter examples | Shows bare expressions instead of 3-column format |
| No pipeline diagram | Users can't see the processing flow at a glance |
| No quick start | Jumps straight into usage without orientation |
| No architecture docs | lib/ modular structure undocumented |
| Snakemake bloat | 40% of README but secondary to CLI |

Existing `defaults/filters.md` is **good** — detailed filter reference, keep as-is.

## Documentation Stack

| Layer | File | Purpose |
|-------|------|---------|
| Landing page | `README.md` | Concise: what, install, quick start, link to docs/ |
| Pipeline guide | `docs/pipeline.md` | Step-by-step pipeline logic with diagram |
| Options reference | `docs/options.md` | Complete flag reference, categorized |
| Filter guide | `docs/filters.md` | Filter system, file format, writing custom filters |
| Snakemake guide | `docs/snakemake.md` | Batch processing workflow |
| Examples cookbook | `docs/examples.md` | Real-world recipes, copy-paste ready |
| Architecture | `docs/architecture.md` | Module structure, error handling, contributing |
| Filter reference | `defaults/filters.md` | Already exists — deep-dive on GATK/Freebayes filters |

**Why MkDocs Material:**
- Markdown-native — write docs in plain `.md`, no JSX/React
- One `mkdocs.yml` config file, zero custom code
- GitHub Actions deploys to Pages automatically on push
- Beautiful rendered site with search, navigation, dark mode out of the box
- Python-based (matches bioinformatics ecosystem)
- Trivial to set up: `pip install mkdocs-material`, add 20-line YAML config
- All docs also readable as raw markdown on GitHub (dual access)

## File-by-file Plan

---

### `README.md` (~80-100 lines)

Concise landing page. Users should understand what hardnormly does and run it in under 60 seconds.

```
# hardnormly

Badges: [CI][license][version]

One-liner: VCF normalization and hard filtering toolkit for WES variant processing.

## Quick Start
  - conda install one-liner
  - Minimal 3-line example
  - "What just happened" sentence

## Pipeline at a Glance
  - ASCII flow diagram (compact, 3 lines)
  - Link to docs/pipeline.md for details

## Installation
  - conda (full env)
  - conda (minimal)
  - apt-get (CI)

## Documentation
  Table linking to each docs/ file:
  | Guide | Description |
  | docs/pipeline.md | Pipeline steps and processing logic |
  | docs/options.md | Complete options reference |
  | docs/filters.md | Filter system and custom filters |
  | docs/examples.md | Real-world usage recipes |
  | docs/snakemake.md | Batch processing with Snakemake |
  | docs/architecture.md | Module structure and development |

## License
  MIT
```

---

### `docs/pipeline.md` (~100-120 lines)

The "how it works" guide. Explains what each step does and why.

```
# Pipeline

## Overview Diagram
  Full ASCII diagram with all 7 steps + optional branches

## Steps

### Step 1: Genome File
  - Auto-fetch from UCSC MySQL or use --genome
  - What the genome file is for (bedtools slop)

### Step 2: BED Region Processing
  - Normalize → merge → compress → index
  - Include vs exclude logic
  - Slop padding explained

### Step 3: VCF Annotation
  - INCLUDE_REGION / EXCLUDE_REGION INFO fields
  - How region filters work downstream

### Step 4: Strip Annotations (optional)
  - --strip-annotations removes INFO fields
  - Use case: remove VEP/SnpEff before filtering

### Step 5: VCF Normalization
  - bcftools norm: multiallelic splitting, left-alignment
  - Reference FASTA requirement

### Step 6: Hard Filtering
  - fill-tags (FORMAT/VAF etc.)
  - Sequential filter stages
  - Soft-filter tagging (FILTER column)
  - --only-pass option

### Step 7: Stats & Plots (optional)
  - bcftools stats output
  - plot-vcfstats (non-fatal)
```

---

### `docs/options.md` (~120-150 lines)

Complete flag reference, organized by purpose.

```
# Options Reference

## Subcommands
  run-pipeline (default), generate-inclusion-bed, generate-exclusion-bed

## run-pipeline Options

### Input / Output
  Table: -v, -f, -o with types, defaults, notes

### Region Filtering
  Table: -b, -e, -g, --genome-build, --slop

### Variant Filtering
  Table: --caller, --filters-file, --filters, --strip-annotations, --only-pass

### Stats & Plots
  Table: --generate-stats, --plot-stats, --plot-output-dir

### Runtime
  Table: --auto-index, --tmp-dir, --no-cleanup, --log-file, --debug

## generate-inclusion-bed Options
  Table: -b, -g, -o, --slop, -v, -h

## generate-exclusion-bed Options
  Table: -e, -o, -v, -h
```

---

### `docs/filters.md` (~80-100 lines)

How the filter system works. Complements `defaults/filters.md` (which has the actual filter definitions).

```
# Filters

## How Filtering Works
  - Soft-filter model: tag, don't remove
  - FILTER column in output VCF
  - --only-pass to remove tagged variants

## Filter File Format
  3-column: name action expression
  Table explaining each column
  Example file content

## Using --caller
  --caller gatk → defaults/gatk_filters.txt
  --caller freebayes → defaults/freebayes_filters.txt
  --filters-file overrides --caller

## Inline Filters
  --filters "name action expression" (repeatable)
  Combined with file filters

## Region-Based Filters
  Auto-generated from BED files:
  - NOT_IN_INCLUDE_REGION
  - IN_EXCLUDE_REGION
  Applied before user-defined filters

## fill-tags
  Pipeline runs bcftools +fill-tags before filtering
  Makes FORMAT/VAF, TYPE, etc. available

## Writing Custom Filters
  Link to bcftools expressions docs
  Tips for common patterns
  Link to defaults/filters.md for GATK/Freebayes reference

## Filter Precedence
  1. Region filters (from BED files)
  2. File-based filters (--filters-file or --caller)
  3. Inline filters (--filters)
  Applied sequentially, each adds to FILTER column
```

---

### `docs/examples.md` (~80-100 lines)

Copy-paste recipes for common use cases.

```
# Examples

## Basic GATK Filtering
  The 80% use case — most users start here

## Freebayes with Exclusion Regions
  Second most common workflow

## Custom Filters (inline)
  Ad-hoc filtering without a file

## Strip VEP Annotations Before Filtering
  Clean up INFO field bloat

## Generate Stats and Plots
  Post-filtering QC

## Generate Exclusion BED from Public Sources
  Helper script usage

## Include + Exclude Regions Together
  Combining both BED types

## Merge BED Files (Subcommands)
  generate-inclusion-bed and generate-exclusion-bed

## Debug a Failed Run
  --debug --no-cleanup workflow
```

---

### `docs/snakemake.md` (~60-80 lines)

Batch processing guide, moved from README bloat.

```
# Snakemake Workflow

## Overview
  What the workflow does (parallel VCF processing)

## Setup
  1. Edit config
  2. Create VCF list
  3. Choose profile

## Running
  Dry run, local, SLURM examples

## Workflow Structure
  Directory tree

## Configuration
  Key config.yaml sections
  Link to schema for validation

## Profiles
  Default vs cluster (Charite example)
```

---

### `docs/architecture.md` (~60-80 lines)

For contributors and advanced users.

```
# Architecture

## Module Structure
  Table: 8 lib/ modules with purpose and key functions

## Pipeline Orchestration
  hardnormly.sh as orchestrator (~220 lines)
  Source order matters (logging first)

## Error Handling
  - set -Eeuo pipefail
  - run_cmd pattern (stderr capture, structured error messages)
  - Pipeline guards (|| { error_msg; return 1 })
  - cleanup_handler EXIT trap

## Development
  make lint, make test
  BATS test structure (7 files, 112 tests)
  Adding a new lib/ module checklist

## CI
  GitHub Actions: lint + test
  Pre-commit hooks
```

---

## Writing Principles

1. **README is a billboard** — 80 lines max, link to details
2. **Examples over explanations** — show, don't tell
3. **Categorized tables** — group flags by purpose, not alphabet
4. **Progressive disclosure** — quick start → guide → reference
5. **Copy-paste ready** — all examples should work with the right files
6. **No jargon without context** — define soft-filter, bcftools expression

## Deliverables

| File | Lines | Status |
|------|-------|--------|
| `README.md` | ~80-100 | Rewrite |
| `docs/index.md` | ~60-80 | New (MkDocs landing) |
| `docs/pipeline.md` | ~100-120 | New |
| `docs/options.md` | ~120-150 | New |
| `docs/filters.md` | ~80-100 | New |
| `docs/examples.md` | ~80-100 | New |
| `docs/snakemake.md` | ~60-80 | New (from README) |
| `docs/architecture.md` | ~60-80 | New |
| `mkdocs.yml` | ~30 | New |
| `.github/workflows/docs.yml` | ~25 | New |
| **Total** | **~740-870** | |

## MkDocs Material Setup

### `mkdocs.yml` (~30 lines)

```yaml
site_name: hardnormly
site_description: VCF normalization and hard filtering toolkit
site_url: https://<user>.github.io/hardnormly/
repo_url: https://github.com/<user>/hardnormly
repo_name: hardnormly

theme:
  name: material
  palette:
    - media: "(prefers-color-scheme: light)"
      scheme: default
      primary: teal
      toggle:
        icon: material/brightness-7
        name: Dark mode
    - media: "(prefers-color-scheme: dark)"
      scheme: slate
      primary: teal
      toggle:
        icon: material/brightness-4
        name: Light mode
  features:
    - navigation.sections
    - navigation.expand
    - content.code.copy
    - search.highlight

nav:
  - Home: index.md
  - Pipeline: pipeline.md
  - Options: options.md
  - Filters: filters.md
  - Examples: examples.md
  - Snakemake: snakemake.md
  - Architecture: architecture.md

markdown_extensions:
  - tables
  - admonitions
  - pymdownx.highlight
  - pymdownx.superfences
  - pymdownx.tabbed:
      alternate_style: true
  - toc:
      permalink: true
```

### `docs/index.md`

Symlink or copy of key README content (MkDocs needs `docs/index.md` as entry point).
Include the quick start, pipeline diagram, and docs navigation table.

### `.github/workflows/docs.yml` (~25 lines)

```yaml
name: docs
on:
  push:
    branches: [main]
    paths: ['docs/**', 'mkdocs.yml', 'README.md']

permissions:
  contents: write

jobs:
  deploy:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4
      - uses: actions/setup-python@v5
        with:
          python-version: '3.12'
      - run: pip install mkdocs-material
      - run: mkdocs gh-deploy --force
```

Deploys to `gh-pages` branch on push to main. Only triggers when docs change (path filter).

### GitHub Pages Setup

One-time: Settings → Pages → Source: "Deploy from a branch" → Branch: `gh-pages` / `/ (root)`.

### Additional Deliverables

| File | Purpose |
|------|---------|
| `mkdocs.yml` | MkDocs config |
| `.github/workflows/docs.yml` | Auto-deploy to Pages |
| `docs/index.md` | Landing page (mirrors README essentials) |

## What NOT to Touch

- `defaults/filters.md` — already complete and accurate
- `--help` output in `lib/cli.sh` — already implemented
- `.planning/` docs — internal, not user-facing
