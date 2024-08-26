# hardnormly

**hardnormly** is a toolkit for VCF normalization and hard filtering. It leverages `bcftools` for variant filtering and `bedtools` for region-based operations on genomic intervals. The script provides a streamlined approach for intersecting BED files, applying VCF normalization, and performing hard filtering based on user-defined criteria.

## Features

- **VCF Normalization**: Normalize variants using a reference FASTA file.
- **Region Filtering**: Filter variants based on BED file regions with optional padding.
- **Hard Filtering**: Apply user-defined filters using `bcftools`.
- **Flexible Input**: Provide filters directly via CLI or read them from a file.

## Dependencies

- `bcftools`
- `bedtools`
- `mysql` (for genome file creation if not provided)

## Usage

```bash
./hardnormly.sh -v <input.vcf.gz> -f <reference.fasta> -b <regions1.bed> -b <regions2.bed> -o <output.vcf.gz> [--filters '<filter_expression>'] [--filters-file <filters.txt>]
```

### Example

```bash
./hardnormly.sh \
  -v input.vcf.gz \
  -f reference.fasta \
  -b regions1.bed \
  -b regions2.bed \
  --filters 'FORMAT/DP<20' \
  --filters 'QUAL<300' \
  -o output.vcf.gz
```

### Options

- \`-v, --vcf <input.vcf.gz>\`: Input VCF file (compressed \`.vcf.gz\` format is recommended).
- \`-f, --fasta <reference.fasta>\`: Reference FASTA file for VCF normalization.
- \`-b, --bed <regions.bed>\`: BED file(s) for region-based filtering (you can specify multiple files).
- \`-g, --genome <genome.file>\`: Genome file for \`bedtools slop\` operation (default: \`hg19.genome\`). This is used to apply padding around BED regions.
- \`-o, --output <output.vcf.gz>\`: Output VCF file (compressed \`.vcf.gz\` format is recommended).
- \`--filters <filter_expression>\`: Inline filter expressions for \`bcftools\`. Multiple filters can be specified by repeating this option.
- \`--filters-file <filters.txt>\`: A file containing filter expressions, one per line. This is an alternative to specifying filters on the command line.

### Filter Expression Example

Filter expressions for \`bcftools\` are used to specify conditions for excluding variants. You can specify them directly in the command line or load them from a file.

**Inline Filter Example:**
```bash
--filters 'FORMAT/DP<20' --filters 'QUAL<300'
```

**Filters File Example:**
The file \`filters.txt\` should contain one filter expression per line:
```
FORMAT/DP<20
QUAL<300
FORMAT/VAF<0.2 && GT!="hom"
```

Then use the \`--filters-file\` option:
```bash
--filters-file filters.txt
```

## Genome File Creation

If the genome file for \`bedtools slop\` is not provided, it can be generated using the UCSC MySQL database:
```bash
mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -e "select chrom, size from hg19.chromInfo" | grep -v "^chrom" | sed 's/chr//g' > hg19.genome
```

This file is required for the padding operation performed by \`bedtools slop\`.

### Exclusion BED Files Documentation

#### Purpose:
Exclusion BED files are used to filter out problematic regions of the genome, such as areas with low mappability, high repeat content, or known technical artifacts. These regions are often flagged to avoid false positives in variant calling.

#### Usage:
To use exclusion BED files, specify them with the `--exclude-bed` option. These files define genomic intervals to exclude during the VCF filtering process. The script will use these intervals to annotate variants and filter out those that overlap with exclusion regions.

Example:
```bash
./hardnormly.sh -v input.vcf.gz -f reference.fasta --exclude-bed exclusion.bed -o output.vcf.gz
```

#### Format:
Exclusion BED files follow the standard BED format, which includes:
1. Chromosome
2. Start Position
3. End Position
4. (Optional) Annotation (e.g., "LowMappability")

Example BED file:
```
chr1    10000   10500   LowMappability
chr1    20000   20500   RepeatRegion
```

#### Resources:
- **[Excluderanges Repository](https://github.com/dozmorovlab/excluderanges):** A collection of exclusion BED files for various species and use cases.
- **[UCSC Encode Downloads](https://hgdownload.cse.ucsc.edu/goldenpath/hg19/encodeDCC/wgEncodeMapability):** Provides access to mappability tracks and other exclusion regions for the hg19 genome.

#### Best Practices:
- Ensure your exclusion BED file is properly formatted and indexed (e.g., using `bgzip` and `tabix`).
- Use exclusion files relevant to your specific genome build and analysis goals to improve variant calling accuracy.

## License

This project is licensed under the MIT License.
