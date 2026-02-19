# Filters

hardnormly uses a **soft-filter** model: variants are tagged in the FILTER column, not removed. This preserves all data and lets downstream tools make their own decisions.

## How Filtering Works

1. Each filter evaluates a bcftools expression against every variant
2. Variants that match an **exclude** filter get the filter name added to FILTER
3. Variants that match an **include** filter keep `PASS`; non-matching ones get tagged
4. Multiple filter tags accumulate (e.g., `DPu10het;lowQUAL`)
5. Use `--only-pass` to remove tagged variants from the output

## Filter File Format

Filter files use three space-separated columns per line:

```
<filter_name> <action> <bcftools_expression>
```

| Column | Values | Description |
|--------|--------|-------------|
| `filter_name` | any string | Tag written to the FILTER column (e.g., `DPu10het`) |
| `action` | `e` or `i` | `e` = exclude (tag matching), `i` = include (keep matching) |
| `expression` | bcftools expr | A [bcftools filter expression](https://samtools.github.io/bcftools/bcftools.html#expressions) |

**Example file** (`defaults/gatk_filters.txt`):

```
DPu10het e FORMAT/DP<10 && GT!="hom"
DPu5hom e FORMAT/DP<5 && GT=="hom"
VAFu02het e FORMAT/VAF<0.2 && GT!="hom"
VAFo08het e FORMAT/VAF>0.8 && GT!="hom"
VAFu095hom e FORMAT/VAF<0.95 && GT=="hom"
gatkSNPhard e TYPE=="SNP" && (AS_FS > 60 || AS_ReadPosRankSum < -8.0 || QUAL < 30.0 || AS_SOR > 3.0 || AS_MQ < 40.0 || AS_MQRankSum < -12.5)
gatkINDELhard e TYPE=="INDEL" && (AS_FS > 200 || AS_ReadPosRankSum < -20.0 || QUAL < 30.0)
```

## Using --caller

The `--caller` flag auto-selects a built-in filter file:

| Value | Filter File | Use When |
|-------|-------------|----------|
| `gatk` | `defaults/gatk_filters.txt` | GATK HaplotypeCaller with allele-specific annotations (`AS_FS`, `AS_SOR`, etc.) |
| `gatk-no-as` | `defaults/gatk_filters_no_as.txt` | GATK output without AS annotations (e.g., Varvis-exported VCFs) |
| `freebayes` | `defaults/freebayes_filters.txt` | Freebayes variant calls |

If both `--caller` and `--filters-file` are provided, `--filters-file` takes precedence (with a warning).

## Inline Filters

Use `--filters` to add filters directly on the command line:

```bash
--filters "lowDP e FORMAT/DP<10"
--filters "highQUAL i QUAL>100"
```

The flag is repeatable. Inline filters are applied **before** file-based filters.

## Region-Based Filters

When BED files are provided (`--include-bed` / `--exclude-bed`), two region filters are auto-generated:

| Filter Name | Condition | Meaning |
|-------------|-----------|---------|
| `NOT_IN_INCLUDE_REGION` | `INFO/INCLUDE_REGION!=1` | Variant is outside all inclusion regions |
| `IN_EXCLUDE_REGION` | `INFO/EXCLUDE_REGION==1` | Variant is inside an exclusion region |

These are applied **before** any user-defined filters.

## fill-tags

Before filtering, the pipeline runs `bcftools +fill-tags` to compute derived fields. This means you can use these in filter expressions even if the original VCF doesn't include them:

- `FORMAT/VAF` — Variant allele frequency
- `TYPE` — Variant type (`SNP`, `INDEL`, `MNP`, etc.)
- Other standard bcftools tags

## Filter Precedence

Filters are applied in this order:

1. **Region filters** — From BED file annotations
2. **Inline filters** — From `--filters` flags
3. **File-based filters** — From `--filters-file` or `--caller`

Each stage adds to the FILTER column. A variant can accumulate multiple tags.

## Writing Custom Filters

Tips for writing effective filter expressions:

- Use `TYPE=="SNP"` or `TYPE=="INDEL"` to target specific variant types
- Combine conditions with `&&` (and) and `||` (or)
- Use parentheses for complex logic
- `GT=="hom"` / `GT!="hom"` to distinguish zygosity
- `FORMAT/` prefix for per-sample fields, `INFO/` for site-level fields

See the [bcftools expressions documentation](https://samtools.github.io/bcftools/bcftools.html#expressions) for the full syntax reference.

For detailed explanations of each GATK and Freebayes filter, see [defaults/filters.md](../defaults/filters.md).
