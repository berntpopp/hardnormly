#!/bin/bash
# scripts/generate_exclusion_bed.sh
# Downloads and merges exclusion regions from four public genomic databases into
# a single exclusion BED file for use with hardnormly --exclude-bed.
#
# Sources:
#   1. ENCODE Blacklist (Boyle-Lab)
#   2. Segmental Duplications (UCSC genomicSuperDups)
#   3. Low Complexity Regions (UCSC RepeatMasker, Low_complexity class)
#   4. Centromeres/Telomeres (UCSC gap table; cytoBandIdeo fallback for hg38)
#
# Usage: bash scripts/generate_exclusion_bed.sh -b hg19 -o ref/exclude_files/hg19_exclusion.bed

set -Eeuo pipefail

# ---------------------------------------------------------------------------
# Usage
# ---------------------------------------------------------------------------
usage() {
	cat <<EOF
Usage: $(basename "$0") -b <build> -o <output.bed> [-v] [-h]

Download and merge genomic exclusion regions into a single BED file.

Options:
  -b, --build   <build>       Genome build: hg19 or hg38 (required)
  -o, --output  <output.bed>  Output BED file path (required)
  -v, --verbose               Print progress messages to stderr
  -h, --help                  Show this help message and exit

Supported builds: hg19, hg38

Sources merged into the output BED:
  1. ENCODE Blacklist v2 (Boyle-Lab/Blacklist)
  2. Segmental Duplications (UCSC genomicSuperDups)
  3. Low Complexity Regions (UCSC RepeatMasker, Low_complexity class)
  4. Centromeres/Telomeres (UCSC gap table; cytoBandIdeo fallback for hg38)

Example:
  bash scripts/generate_exclusion_bed.sh -b hg19 -o ref/exclude_files/hg19_exclusion.bed
  bash scripts/generate_exclusion_bed.sh -b hg38 -o ref/exclude_files/hg38_exclusion.bed
EOF
}

# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------
build=""
output=""
verbose=false

while [[ $# -gt 0 ]]; do
	case "$1" in
		-b | --build)
			build="$2"
			shift 2
			;;
		-o | --output)
			output="$2"
			shift 2
			;;
		-v | --verbose)
			verbose=true
			shift
			;;
		-h | --help)
			usage
			exit 0
			;;
		*)
			echo "ERROR: Unknown option: $1" >&2
			usage >&2
			exit 1
			;;
	esac
done

# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------
if [[ -z "$build" ]]; then
	echo "ERROR: --build is required (hg19 or hg38)" >&2
	usage >&2
	exit 1
fi

if [[ "$build" != "hg19" && "$build" != "hg38" ]]; then
	echo "ERROR: Unknown genome build '$build'. Supported builds: hg19, hg38" >&2
	exit 1
fi

if [[ -z "$output" ]]; then
	echo "ERROR: --output is required" >&2
	usage >&2
	exit 1
fi

output_dir="$(dirname "$output")"
if [[ ! -d "$output_dir" ]]; then
	echo "ERROR: Output directory does not exist: $output_dir" >&2
	exit 1
fi

# ---------------------------------------------------------------------------
# Temporary directory with cleanup trap
# ---------------------------------------------------------------------------
tmp_dir=$(mktemp -d -t gen_excl_bed-XXXXXXXXXX)

cleanup_handler() {
	local exit_code=$?
	trap - EXIT
	rm -rf "$tmp_dir" 2>/dev/null || true
	exit "$exit_code"
}

trap cleanup_handler EXIT

# ---------------------------------------------------------------------------
# Helper: log progress to stderr in verbose mode
# ---------------------------------------------------------------------------
log_verbose() {
	if [[ "$verbose" == true ]]; then
		echo "$*" >&2
	fi
}

# ---------------------------------------------------------------------------
# Helper: warn if a downloaded BED file is empty
# ---------------------------------------------------------------------------
warn_if_empty() {
	local bed_file="$1"
	local source_name="$2"
	local line_count
	line_count=$(wc -l <"$bed_file")
	if [[ "$line_count" -eq 0 ]]; then
		echo "WARNING: $source_name produced 0 regions — source may be unavailable for $build" >&2
	fi
}

# ---------------------------------------------------------------------------
# Source 1: ENCODE Blacklist v2
# ---------------------------------------------------------------------------
log_verbose "Downloading ENCODE Blacklist v2..."
wget -q -O - \
	"https://github.com/Boyle-Lab/Blacklist/raw/master/lists/${build}-blacklist.v2.bed.gz" \
	| gunzip -c \
	| awk '{OFS="\t"; print $1, $2, $3}' \
	| bedtools sort -i - \
		>"${tmp_dir}/blacklist.bed"
warn_if_empty "${tmp_dir}/blacklist.bed" "ENCODE Blacklist"

# ---------------------------------------------------------------------------
# Source 2: Segmental Duplications (UCSC genomicSuperDups)
# ---------------------------------------------------------------------------
log_verbose "Downloading segmental duplications (UCSC genomicSuperDups)..."
wget -q -O - \
	"https://hgdownload.soe.ucsc.edu/goldenPath/${build}/database/genomicSuperDups.txt.gz" \
	| gunzip -c \
	| grep -v '^#' \
	| awk '{OFS="\t"; print $2, $3, $4}' \
	| bedtools sort -i - \
		>"${tmp_dir}/superdups.bed"
warn_if_empty "${tmp_dir}/superdups.bed" "Segmental Duplications"

# ---------------------------------------------------------------------------
# Source 3: Low Complexity Regions (UCSC RepeatMasker)
# ---------------------------------------------------------------------------
log_verbose "Downloading low complexity regions (UCSC RepeatMasker)..."
wget -q -O - \
	"https://hgdownload.soe.ucsc.edu/goldenPath/${build}/database/rmsk.txt.gz" \
	| gunzip -c \
	| grep -v '^#' \
	| awk '$12 == "Low_complexity" {OFS="\t"; print $6, $7, $8}' \
	| bedtools sort -i - \
		>"${tmp_dir}/lcr.bed"
warn_if_empty "${tmp_dir}/lcr.bed" "Low Complexity Regions"

# ---------------------------------------------------------------------------
# Source 4: Centromeres/Telomeres (UCSC gap table; cytoBandIdeo fallback)
# ---------------------------------------------------------------------------
log_verbose "Downloading centromeres/telomeres (UCSC gap table)..."
wget -q -O - \
	"https://hgdownload.soe.ucsc.edu/goldenPath/${build}/database/gap.txt.gz" \
	| gunzip -c \
	| grep -v '^#' \
	| awk '$8 ~ /centromere|telomere/ {OFS="\t"; print $2, $3, $4}' \
	| bedtools sort -i - \
		>"${tmp_dir}/gaps.bed"

# For hg38, the gap table often has very few centromere entries; fall back to cytoBandIdeo
if [[ "$build" == "hg38" ]]; then
	local_gap_count=$(wc -l <"${tmp_dir}/gaps.bed")
	if [[ "$local_gap_count" -lt 10 ]]; then
		log_verbose "Gap table had $local_gap_count rows; using cytoBandIdeo for centromeres (hg38 fallback)..."
		wget -q -O - \
			"https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/cytoBandIdeo.txt.gz" \
			| gunzip -c \
			| grep -v '^#' \
			| awk '$5 == "acen" {OFS="\t"; print $1, $2, $3}' \
			| bedtools sort -i - \
				>"${tmp_dir}/gaps.bed"
	fi
fi
warn_if_empty "${tmp_dir}/gaps.bed" "Centromeres/Telomeres"

# ---------------------------------------------------------------------------
# Sanity check: fail if ALL sources are empty
# ---------------------------------------------------------------------------
total_regions=0
for bed in "${tmp_dir}/blacklist.bed" "${tmp_dir}/superdups.bed" \
	"${tmp_dir}/lcr.bed" "${tmp_dir}/gaps.bed"; do
	total_regions=$((total_regions + $(wc -l <"$bed")))
done

if [[ "$total_regions" -eq 0 ]]; then
	echo "ERROR: All sources produced 0 regions — downloads may have failed" >&2
	exit 1
fi

# ---------------------------------------------------------------------------
# Merge: sort all sources and merge overlapping intervals
# ---------------------------------------------------------------------------
log_verbose "Merging and sorting all exclusion regions into $output..."
cat "${tmp_dir}/blacklist.bed" "${tmp_dir}/superdups.bed" \
	"${tmp_dir}/lcr.bed" "${tmp_dir}/gaps.bed" \
	| bedtools sort -i - \
	| bedtools merge -i - \
		>"$output"

log_verbose "Done. Output: $output"
