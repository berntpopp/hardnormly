# Options Reference

## Subcommands

| Subcommand | Description |
|------------|-------------|
| `run-pipeline` | Normalize and filter a VCF file (default if omitted) |
| `generate-inclusion-bed` | Merge BED files into a combined inclusion region |
| `generate-exclusion-bed` | Merge BED files into a combined exclusion region |

Running `hardnormly.sh [options]` without a subcommand is equivalent to `hardnormly.sh run-pipeline [options]`.

---

## run-pipeline

### Input / Output

| Flag | Type | Default | Description |
|------|------|---------|-------------|
| `-v, --vcf` | FILE | *required* | Input VCF file (`.vcf` or `.vcf.gz`) |
| `-f, --fasta` | FILE | *required* | Reference FASTA file for normalization |
| `-o, --output` | FILE | stdout | Output VCF file. Extension determines format: `.vcf.gz` writes compressed, `.vcf` writes plain text |

### Region Filtering

| Flag | Type | Default | Description |
|------|------|---------|-------------|
| `-b, --include-bed` | FILE | — | BED file for inclusion regions. Repeatable: multiple files are intersected |
| `-e, --exclude-bed` | FILE | — | BED file for exclusion regions. Repeatable: multiple files are unioned |
| `-g, --genome` | FILE | — | Genome file for `bedtools slop`. Skips UCSC MySQL fetch if provided |
| `--genome-build` | STRING | `hg19` | Genome build for UCSC fetch (e.g., `hg19`, `hg38`) |
| `--slop` | INT | `20` | Region padding in base pairs applied to inclusion BED regions |

### Variant Filtering

| Flag | Type | Default | Description |
|------|------|---------|-------------|
| `--caller` | STRING | — | Auto-select filter file: `gatk` or `freebayes` |
| `--filters-file` | FILE | — | Custom filter file (overrides `--caller` if both provided) |
| `--filters` | STRING | — | Inline filter: `"name action expression"`. Repeatable |
| `--strip-annotations` | STRING | — | Comma-separated INFO fields to remove before filtering (e.g., `INFO/CSQ,INFO/ANN`) |
| `--only-pass` | FLAG | off | Remove filtered variants, keep only PASS in output |

### Stats and Plots

| Flag | Type | Default | Description |
|------|------|---------|-------------|
| `--generate-stats` | FLAG | off | Generate `bcftools stats` file alongside output |
| `--plot-stats` | FLAG | off | Plot stats (requires `--generate-stats` and `--plot-output-dir`) |
| `--plot-output-dir` | DIR | — | Directory for plot output files |

### Runtime

| Flag | Type | Default | Description |
|------|------|---------|-------------|
| `--auto-index` | FLAG | off | Auto-index output if compressed (`.vcf.gz` → `.vcf.gz.tbi`) |
| `--tmp-dir` | DIR | auto | Custom temporary directory (default: system mktemp) |
| `--no-cleanup` | FLAG | off | Preserve temporary directory after run |
| `--log-file` | FILE | — | Write log messages to file instead of stdout |
| `--debug` | FLAG | off | Enable verbose debug output (`set -x` tracing) |

### Info

| Flag | Description |
|------|-------------|
| `--version` | Print version and exit |
| `-h, --help` | Print usage help and exit |

---

## generate-inclusion-bed

Merge one or more BED files into a combined inclusion region file with slop padding.

| Flag | Type | Default | Description |
|------|------|---------|-------------|
| `-b, --include-bed` | FILE | *required* | Input BED file. Repeatable |
| `-g, --genome` | FILE | *required* | Genome file for slop operation |
| `-o, --output` | FILE | *required* | Output merged BED file |
| `--slop` | INT | `20` | Region padding in base pairs |
| `-v, --verbose` | FLAG | off | Show progress messages |
| `-h, --help` | FLAG | — | Show subcommand help |

---

## generate-exclusion-bed

Merge one or more BED files into a combined exclusion region file.

| Flag | Type | Default | Description |
|------|------|---------|-------------|
| `-e, --exclude-bed` | FILE | *required* | Input BED file. Repeatable |
| `-o, --output` | FILE | *required* | Output merged BED file |
| `-v, --verbose` | FLAG | off | Show progress messages |
| `-h, --help` | FLAG | — | Show subcommand help |
