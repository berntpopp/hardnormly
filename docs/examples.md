# Examples

Copy-paste recipes for common workflows. All examples assume you have `bcftools`, `bedtools`, and `htslib` installed.

## Basic GATK Filtering

The 80% use case — normalize and apply GATK hard filters:

```bash
./hardnormly.sh run-pipeline \
  -v sample.vcf.gz \
  -f reference.fasta \
  --caller gatk \
  -o sample.filtered.vcf.gz
```

Output: soft-filtered VCF with tags like `DPu10het`, `gatkSNPhard`, etc. in the FILTER column.

## Freebayes with Exclusion Regions

Filter Freebayes calls while excluding problematic genomic regions:

```bash
./hardnormly.sh run-pipeline \
  -v sample.freebayes.vcf.gz \
  -f reference.fasta \
  --caller freebayes \
  -e ref/exclude_files/hg19_exclusion.bed \
  -o sample.filtered.vcf.gz
```

Variants inside exclusion regions get the `IN_EXCLUDE_REGION` tag.

## Include and Exclude Regions Together

Use a capture kit BED (include) and a blacklist BED (exclude):

```bash
./hardnormly.sh run-pipeline \
  -v sample.vcf.gz \
  -f reference.fasta \
  -b capture_kit.bed \
  -e blacklist.bed \
  -g hg19.genome \
  --caller gatk \
  --slop 50 \
  -o sample.filtered.vcf.gz
```

Variants outside the capture kit (padded by 50 bp) get `NOT_IN_INCLUDE_REGION`. Variants inside the blacklist get `IN_EXCLUDE_REGION`.

## Custom Inline Filters

Add ad-hoc filters without creating a file:

```bash
./hardnormly.sh run-pipeline \
  -v sample.vcf.gz \
  -f reference.fasta \
  --filters "lowDP e FORMAT/DP<10" \
  --filters "lowQUAL e QUAL<30" \
  --filters "highVAF i FORMAT/VAF>0.2" \
  -o sample.filtered.vcf.gz
```

## Strip VEP Annotations Before Filtering

Clean up INFO field bloat from VEP or SnpEff before filtering:

```bash
./hardnormly.sh run-pipeline \
  -v annotated.vcf.gz \
  -f reference.fasta \
  --strip-annotations INFO/CSQ,INFO/ANN \
  --caller gatk \
  -o clean.filtered.vcf.gz
```

## Keep Only PASS Variants

Remove all filtered variants, outputting only those that passed every filter:

```bash
./hardnormly.sh run-pipeline \
  -v sample.vcf.gz \
  -f reference.fasta \
  --caller gatk \
  --only-pass \
  -o sample.pass_only.vcf.gz
```

## Generate Stats and Plots

Produce QC summaries after filtering:

```bash
./hardnormly.sh run-pipeline \
  -v sample.vcf.gz \
  -f reference.fasta \
  --caller gatk \
  --generate-stats \
  --plot-stats \
  --plot-output-dir qc_plots/ \
  -o sample.filtered.vcf.gz
```

Creates `sample.filtered.stats.txt` and visual plots in `qc_plots/`.

## Generate Exclusion BED from Public Sources

Create a comprehensive exclusion BED file from ENCODE blacklist, segmental duplications, low complexity, and centromere/telomere regions:

```bash
bash scripts/generate_exclusion_bed.sh -b hg19 -o ref/exclude_files/hg19_exclusion.bed -v
bash scripts/generate_exclusion_bed.sh -b hg38 -o ref/exclude_files/hg38_exclusion.bed -v
```

Use the BED matching your reference build. UCSC sources use `chr` contig names; if your VCF/reference uses `1`, `2`, etc., strip the prefix before passing the BED to hardnormly.

## Merge BED Files (Subcommands)

Pre-merge BED files for reuse across multiple samples:

```bash
# Merge inclusion BEDs with 50bp padding
./hardnormly.sh generate-inclusion-bed \
  -b kit_v1.bed \
  -b kit_v2.bed \
  -g hg19.genome \
  --slop 50 \
  -o merged_targets.bed \
  -v

# Merge exclusion BEDs
./hardnormly.sh generate-exclusion-bed \
  -e blacklist.bed \
  -e segdups.bed \
  -o merged_exclusions.bed \
  -v
```

## Debug a Failed Run

When something goes wrong, use debug mode and preserve temp files:

```bash
./hardnormly.sh run-pipeline \
  -v sample.vcf.gz \
  -f reference.fasta \
  --caller gatk \
  --debug \
  --no-cleanup \
  --tmp-dir /tmp/hardnormly-debug \
  --log-file debug.log \
  -o sample.filtered.vcf.gz
```

- `--debug` enables `set -x` tracing (every command is printed)
- `--no-cleanup` preserves intermediate files in the temp directory
- `--log-file` writes all log messages to a file for later review
- `--tmp-dir` uses a predictable path so you can find the files
