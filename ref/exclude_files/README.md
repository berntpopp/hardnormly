# Exclusion BED files

Generated exclusion BED files belong here, but large build-specific BED outputs are not committed.

Generate the file matching your reference build:

```bash
bash scripts/generate_exclusion_bed.sh -b hg19 -o ref/exclude_files/hg19_exclusion.bed -v
bash scripts/generate_exclusion_bed.sh -b hg38 -o ref/exclude_files/hg38_exclusion.bed -v
```

The generator downloads public exclusion regions from ENCODE/Boyle-Lab and UCSC, then sorts and merges them into one BED file.

Use only an exclusion BED that matches both:

- Genome build: `hg19`/GRCh37 vs `hg38`/GRCh38
- Contig naming: `chr1` vs `1`

The downloaded UCSC-style sources use `chr` contig names. If your VCF/reference uses contigs without `chr`, create a matching copy, for example:

```bash
sed 's/^chr//' ref/exclude_files/hg38_exclusion.bed > ref/exclude_files/hg38_exclusion_nochr.bed
```
