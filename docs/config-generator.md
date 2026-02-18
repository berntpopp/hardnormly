# Config Generator

`scripts/generate_config.sh` is an interactive helper that generates ready-to-use `config.yaml` and `samples.tsv` files for the three pipelines in the Scholl Lab processing stack:

| Pipeline | Direction | Script |
|----------|-----------|--------|
| **sm-alignment** | FASTQ → analysis-ready BAM | [sm-alignment](https://github.com/scholl-lab/sm-alignment) |
| **sm-vcf-calling** | BAM → filtered VCF | [sm-vcf-calling](https://github.com/scholl-lab/sm-vcf-calling) |
| **hardnormly** | VCF → normalized/hard-filtered VCF | this repo |

## Quick Start

```bash
# Interactive — configure all three pipelines, write to ./configs/
scripts/generate_config.sh -o configs/

# Configure only the hardnormly workflow
scripts/generate_config.sh -p hardnormly -o config/

# Non-interactive (accept all defaults, write templates)
scripts/generate_config.sh -p hardnormly -n -o config/
```

## Options

| Flag | Description |
|------|-------------|
| `-p, --pipeline PIPELINE` | Pipeline to configure: `alignment`, `calling`, `hardnormly`, `all` (default: `all`). Repeatable. |
| `-o, --output-dir DIR` | Directory to write generated files (default: current directory) |
| `-n, --non-interactive` | Accept all defaults without prompting |
| `-y, --yes` | Skip the write-confirmation prompt |
| `-h, --help` | Show help and exit |

## Output Files

| File | Target |
|------|--------|
| `sm-alignment-config.yaml` | Copy to `sm-alignment/config/config.yaml` |
| `sm-alignment-samples.tsv` | Fill in, copy to `sm-alignment/config/samples.tsv` |
| `sm-calling-config.yaml` | Copy to `sm-vcf-calling/config/config.yaml` |
| `sm-calling-samples.tsv` | Fill in, copy to `sm-vcf-calling/config/samples.tsv` |
| `hardnormly-config.yaml` | Copy to `config/config.yaml` in this repo |

---

## sm-alignment Settings

Prompted settings for the FASTQ → BAM alignment pipeline:

### Reference

| Setting | Description | Default |
|---------|-------------|---------|
| Reference genome FASTA | Uncompressed reference (required) | — |
| Genome build | `GRCh37` or `GRCh38` | `GRCh38` |
| Known variant sites | VCFs for BQSR (dbSNP, Mills, 1000G) — one per line | — |

### Paths

| Setting | Description | Default |
|---------|-------------|---------|
| Samples TSV | Path to sample metadata file | `config/samples.tsv` |
| FASTQ directory | Directory containing input FASTQ files | — |
| Output directory | Output directory for BAMs | `results/alignment` |
| Log subdirectory | Subdirectory name for log files | `logs` |

### FASTQ / Read Group

| Setting | Description | Default |
|---------|-------------|---------|
| R1 suffix | Filename suffix for read 1 | `_R1_001.fastq.gz` |
| R2 suffix | Filename suffix for read 2 | `_R2_001.fastq.gz` |
| Platform | Sequencing platform for `PL` read group tag | `ILLUMINA` |

### Trimming

| Setting | Description | Default |
|---------|-------------|---------|
| Enable BBDuk trimming | Adapter/quality trimming with BBDuk | `false` |

### QC

| Setting | Description | Default |
|---------|-------------|---------|
| QC master switch | Enable/disable all QC rules | `true` |
| FastQC | Run FastQC on raw FASTQs | `true` |
| samtools stats | Run `samtools stats` on final BAM | `true` |
| samtools flagstat | Run `samtools flagstat` on final BAM | `true` |
| Picard metrics | Run `CollectMultipleMetrics` | `true` |
| Qualimap bamqc | Deep coverage/quality analysis (high memory) | `false` |

### sm-alignment samples.tsv

Each row represents one FASTQ pair. Rows sharing the same `project_sample` value are merged into a single BAM.

| Column | Required | Description |
|--------|----------|-------------|
| `fastq_files_basename` | Yes | Filename prefix before `_R1_001.fastq.gz` |
| `lane` | Yes | Sequencing lane (e.g., `L001`) — used in read group ID |
| `project_sample` | Yes | Logical sample name — rows are merged by this value |
| `mdc_project` | Yes | Project identifier for `PU` read group tag |
| `subfolder` | No | Subdirectory within FASTQ folder |

---

## sm-vcf-calling Settings

### Caller

| Choice | Description |
|--------|-------------|
| `mutect2` | GATK Mutect2 (somatic, tumor-normal, tumor-only, germline) |
| `freebayes` | FreeBayes (germline) |
| `all` | Run both callers |

### Reference & GATK Resources

| Setting | Description | Default |
|---------|-------------|---------|
| Reference genome FASTA | Path to reference | — |
| Genome build | `GRCh37` or `GRCh38` | `GRCh38` |
| Panel of Normals VCF | Required for Mutect2 | — |
| gnomAD AF VCF | Required for Mutect2 | — |
| gnomAD common biallelic VCF | Required for Mutect2 | — |

### Paths

| Setting | Description | Default |
|---------|-------------|---------|
| Samples TSV | Path to sample metadata file | `config/samples.tsv` |
| BAM directory | Directory containing input BAMs | — |
| Output directory | Output directory for VCFs | `results/calling` |
| BAM extension | BAM suffix from sm-alignment | `.merged.dedup.bqsr.bam` |

### Scatter Strategy

| Mode | Description |
|------|-------------|
| `chromosome` | One Mutect2 job per chromosome (fast, recommended) |
| `interval` | Fixed-size interval scatter (configure `count`) |
| `none` | No scatter — single job per sample |

### PureCN (optional)

| Setting | Description | Default |
|---------|-------------|---------|
| Enable PureCN | Copy number analysis using PureCN | `false` |
| Genome | `hg19` or `hg38` | `hg38` |
| Capture bait BED | Required when PureCN is enabled | — |
| normalDB.rds | Pre-built normal database (optional) | — |

### sm-vcf-calling samples.tsv

| Column | Required | Description |
|--------|----------|-------------|
| `sample` | Yes | Unique identifier used in output filenames |
| `tumor_bam` | Yes | BAM basename without extension |
| `normal_bam` | No | Matched normal BAM basename, or `.` if none |
| `analysis_type` | Yes | `tumor_only`, `tumor_normal`, or `germline` |

---

## hardnormly Workflow Settings

### Reference

| Setting | Description | Default |
|---------|-------------|---------|
| Reference FASTA | Path to reference FASTA | — |
| Genome file | Chromosome sizes for `bedtools slop` | `defaults/hg19.genome` |
| Genome build | `GRCh37` or `GRCh38` | `GRCh37` |

### Paths

| Setting | Description | Default |
|---------|-------------|---------|
| VCF list | File with one input VCF path per line | `input/vcfs.txt` |
| Output directory | Output directory for processed VCFs | `results/hardnormly` |

### Regions

| Setting | Description | Default |
|---------|-------------|---------|
| Include BED files | Target capture / exome regions (one per line) | — |
| Exclude BED files | Blacklist / low-complexity regions (one per line) | — |
| Slop | Base pairs to pad include regions | `100` |

!!! tip "Generating exclude BED files"
    Use `scripts/generate_exclusion_bed.sh` or the `generate-exclusion-bed` subcommand to merge public blacklist databases into a single exclusion file.

### Filtering

| Setting | Description | Default |
|---------|-------------|---------|
| Filter source | `caller` preset, `file`, or `none` | `caller` |
| Caller preset | `gatk` or `freebayes` — selects built-in filter set | `gatk` |
| Custom filters file | Path to a 3-column filter TSV (used when source=`file`) | — |
| Strip annotations | Comma-separated `INFO/` fields to remove before filtering (e.g. `INFO/CSQ,INFO/ANN`) | — |
| Only PASS | Remove soft-filtered variants from output | `false` |

See [Filters](filters.md) for the full filter expression format and built-in presets.

### Processing

| Setting | Description | Default |
|---------|-------------|---------|
| Generate stats | Run `bcftools stats` per sample | `true` |
| Auto-index | Index output VCF with `tabix` | `true` |
| Plot stats | Generate plots from stats (requires `plot-vcfstats`) | `false` |
| Plot output directory | Base directory for per-sample plots | `{output_folder}/plots` |

---

## Full Pipeline Example

```bash
# Step 1 — generate all three configs interactively
scripts/generate_config.sh -o project-configs/

# Step 2 — edit sample sheets
nano project-configs/sm-alignment-samples.tsv
nano project-configs/sm-calling-samples.tsv

# Step 3 — run sm-alignment
cp project-configs/sm-alignment-config.yaml sm-alignment/config/config.yaml
cd sm-alignment
snakemake --snakefile workflow/Snakefile --workflow-profile profiles/default

# Step 4 — run sm-vcf-calling (after alignment finishes)
cp project-configs/sm-calling-config.yaml sm-vcf-calling/config/config.yaml
cd ../sm-vcf-calling
snakemake --snakefile workflow/Snakefile --workflow-profile profiles/default

# Step 5 — run hardnormly normalization + filtering
ls results/calling/**/*.vcf.gz > input/vcfs.txt
cp project-configs/hardnormly-config.yaml config/config.yaml
cd ../hardnormly
snakemake --snakefile workflow/Snakefile --configfile config/config.yaml
```
