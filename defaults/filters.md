# Hard Filter Definitions

hardnormly uses caller-specific filter files to apply hard filters via `bcftools filter`. Each line in a filter file has three space-separated fields:

```
<filter_name> <action> <bcftools_expression>
```

- **filter_name**: Label applied to the FILTER column in the output VCF.
- **action**: `e` (exclude — soft-filter variants matching the expression) or `i` (include — soft-filter variants *not* matching the expression).
- **bcftools_expression**: A valid [bcftools expression](https://samtools.github.io/bcftools/bcftools.html#expressions).

All filters are applied as soft filters (`-m+`), meaning matching variants are tagged in the FILTER column but not removed. Use `--only-pass` to retain only PASS variants after filtering.

Note: The pipeline runs `bcftools +fill-tags` before filtering, so computed tags like `FORMAT/VAF` are available in filter expressions regardless of the caller.

---

## Shared Genotype-Level Filters

These filters are common to both GATK and Freebayes filter sets. They operate on genotype-level fields that are caller-agnostic.

| Filter | Expression | Description |
|--------|-----------|-------------|
| DPu10het | `FORMAT/DP<10 && GT!="hom"` | Exclude heterozygous calls with depth below 10 |
| DPu5hom | `FORMAT/DP<5 && GT=="hom"` | Exclude homozygous calls with depth below 5 |
| VAFu02het | `FORMAT/VAF<0.2 && GT!="hom"` | Exclude het calls with variant allele fraction below 0.2 |
| VAFo08het | `FORMAT/VAF>0.8 && GT!="hom"` | Exclude het calls with variant allele fraction above 0.8 |
| VAFu095hom | `FORMAT/VAF<0.95 && GT=="hom"` | Exclude hom calls with variant allele fraction below 0.95 |

---

## GATK HaplotypeCaller Filters (`gatk_filters.txt`)

Caller-specific filters based on [GATK hard filtering recommendations](https://gatk.broadinstitute.org/hc/en-us/articles/360035890471-Hard-filtering-germline-short-variants).

| Filter | Expression | Description |
|--------|-----------|-------------|
| gatkSNPhard | `TYPE=="SNP" && (AS_FS > 60 \|\| AS_ReadPosRankSum < -8.0 \|\| QUAL < 30.0 \|\| AS_SOR > 3.0 \|\| AS_MQ < 40.0 \|\| AS_MQRankSum < -12.5)` | GATK recommended SNP hard filters: strand bias (Fisher/SOR), read position bias, mapping quality, and minimum QUAL |
| gatkINDELhard | `TYPE=="INDEL" && (AS_FS > 200 \|\| AS_ReadPosRankSum < -20.0 \|\| QUAL < 30.0)` | GATK recommended INDEL hard filters: strand bias, read position bias, and minimum QUAL |

### GATK annotation fields used

- **AS_FS**: Allele-specific Fisher strand bias (phred-scaled p-value)
- **AS_SOR**: Allele-specific strand odds ratio
- **AS_MQ**: Allele-specific root mean square mapping quality
- **AS_MQRankSum**: Allele-specific rank sum test for mapping qualities (ref vs alt)
- **AS_ReadPosRankSum**: Allele-specific rank sum test for read position bias (ref vs alt)
- **QUAL**: Phred-scaled variant quality score

---

## Freebayes Filters (`freebayes_filters.txt`)

Caller-specific filters based on [Freebayes recommended filtering](https://github.com/freebayes/freebayes) criteria from Erik Garrison and community best practices.

| Filter | Expression | Description |
|--------|-----------|-------------|
| lowQUAL | `QUAL<20` | Exclude variants with quality below phred 20 (>1% error probability) |
| QUALperAO | `QUAL / INFO/AO < 10` | Exclude variants where quality per alternate observation is below 10 |
| strandBias | `INFO/SAF==0 \|\| INFO/SAR==0` | Exclude variants where the alt allele is only observed on one strand |
| readPosBias | `INFO/RPR<=1 && INFO/RPL<=1` | Exclude variants where supporting reads are only placed on one side |

### Freebayes annotation fields used

- **QUAL**: Phred-scaled probability that a polymorphism exists at the site
- **INFO/AO**: Number of alternate allele observations
- **INFO/SAF**: Number of alternate observations on the forward strand
- **INFO/SAR**: Number of alternate observations on the reverse strand
- **INFO/RPR**: Number of reads supporting the variant placed to the right
- **INFO/RPL**: Number of reads supporting the variant placed to the left

---

## Usage

Specify the appropriate filter file for your variant caller in the config or on the command line:

```bash
# GATK HaplotypeCaller
./hardnormly.sh -v input.vcf.gz -f ref.fasta --filters-file defaults/gatk_filters.txt ...

# Freebayes
./hardnormly.sh -v input.vcf.gz -f ref.fasta --filters-file defaults/freebayes_filters.txt ...
```

Or in `snakemake/config.yaml`:

```yaml
filters_file: "defaults/gatk_filters.txt"   # for GATK
# filters_file: "defaults/freebayes_filters.txt"  # for Freebayes
```
