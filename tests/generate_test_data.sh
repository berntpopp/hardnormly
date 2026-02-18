#!/usr/bin/env bash
# generate_test_data.sh — Download and subset public data into tests/data/real/
#
# Usage:
#   ./tests/generate_test_data.sh [--force]
#
# Options:
#   --force    Re-download all files even if they already exist
#
# Dependencies: bcftools, samtools, tabix, wget
# Runtime: several minutes (downloads ~140MB GIAB + 1000G streaming + freebayes tiny)
#
# This script is idempotent: re-running it without --force skips existing files.

set -euo pipefail

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA_DIR="${SCRIPT_DIR}/data"
REAL_DIR="${DATA_DIR}/real"
REGION="22:16000000-16100000"
FORCE=false

# Parse arguments
for arg in "$@"; do
	case "$arg" in
		--force) FORCE=true ;;
		-h | --help)
			grep '^#' "$0" | grep -v '#!/' | sed 's/^# //' | sed 's/^#//'
			exit 0
			;;
		*)
			echo "Unknown argument: $arg" >&2
			echo "Usage: $0 [--force]" >&2
			exit 1
			;;
	esac
done

# ---------------------------------------------------------------------------
# Dependency check
# ---------------------------------------------------------------------------
for tool in bcftools samtools tabix wget; do
	command -v "$tool" >/dev/null 2>&1 || {
		echo "ERROR: '$tool' not found in PATH" >&2
		echo "Install via: conda activate hardnormly" >&2
		exit 1
	}
done

# ---------------------------------------------------------------------------
# Setup
# ---------------------------------------------------------------------------
mkdir -p "${REAL_DIR}"

TMPDIR=$(mktemp -d -t testdata-XXXXXXXXXX)
trap 'rm -rf "$TMPDIR"' EXIT

echo "=== Generating test data ===" >&2
echo "Output directory: ${REAL_DIR}" >&2
echo "Force re-download: ${FORCE}" >&2
echo "" >&2

# ---------------------------------------------------------------------------
# GIAB NA12878 benchmark — VCF and high-confidence BED
# ---------------------------------------------------------------------------
echo "[1/4] GIAB NA12878 benchmark..." >&2

GIAB_BASE="https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/NA12878_HG001/latest/GRCh37"
GIAB_OUT="${REAL_DIR}/giab_NA12878_chr22_16M.vcf.gz"
GIAB_BED="${REAL_DIR}/giab_NA12878_chr22_16M_highconf.bed"

if [[ "$FORCE" == true ]] || [[ ! -f "${GIAB_OUT}" ]]; then
	echo "  Downloading GIAB benchmark VCF (~138MB)..." >&2
	wget -q --show-progress \
		"${GIAB_BASE}/HG001_GRCh37_1_22_v4.2.1_benchmark.vcf.gz" \
		-O "${TMPDIR}/giab.vcf.gz"
	wget -q \
		"${GIAB_BASE}/HG001_GRCh37_1_22_v4.2.1_benchmark.vcf.gz.tbi" \
		-O "${TMPDIR}/giab.vcf.gz.tbi"

	echo "  Subsetting to ${REGION}..." >&2
	bcftools view -r "${REGION}" "${TMPDIR}/giab.vcf.gz" -Oz \
		-o "${GIAB_OUT}"
	bcftools index -t "${GIAB_OUT}"
	echo "  Created: ${GIAB_OUT}" >&2
else
	echo "  Skipping GIAB VCF (exists; use --force to re-download)" >&2
fi

if [[ "$FORCE" == true ]] || [[ ! -f "${GIAB_BED}" ]]; then
	echo "  Downloading GIAB high-confidence BED..." >&2
	wget -q \
		"${GIAB_BASE}/HG001_GRCh37_1_22_v4.2.1_benchmark.bed" \
		-O "${TMPDIR}/giab_benchmark.bed"
	awk '$1=="22" && $3>=16000000 && $2<=16100000' \
		"${TMPDIR}/giab_benchmark.bed" >"${GIAB_BED}"
	echo "  Created: ${GIAB_BED}" >&2
else
	echo "  Skipping GIAB BED (exists; use --force to re-download)" >&2
fi

# ---------------------------------------------------------------------------
# 1000 Genomes Phase 3 — single-sample subsets via remote tabix streaming
# ---------------------------------------------------------------------------
echo "" >&2
echo "[2/4] 1000 Genomes Phase 3 samples..." >&2

KG_URL="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"

for SAMPLE in NA12878 NA19247 HG00096; do
	OUTFILE="${REAL_DIR}/1kg_${SAMPLE}_chr22_16M.vcf.gz"
	if [[ "$FORCE" == true ]] || [[ ! -f "$OUTFILE" ]]; then
		echo "  Streaming ${SAMPLE} from 1000G (remote tabix)..." >&2
		tabix -h "${KG_URL}" "${REGION}" \
			| bcftools view -s "${SAMPLE}" --min-ac 1 \
				-Oz -o "${OUTFILE}"
		bcftools index -t "${OUTFILE}"
		echo "  Created: ${OUTFILE}" >&2
	else
		echo "  Skipping 1kg_${SAMPLE} (exists; use --force to re-download)" >&2
	fi
done

# ---------------------------------------------------------------------------
# Freebayes tiny test data — authentic Freebayes output with real INFO fields
# ---------------------------------------------------------------------------
echo "" >&2
echo "[3/4] Freebayes tiny test data..." >&2

FB_RAW="https://raw.githubusercontent.com/freebayes/freebayes/master/test/tiny"
FB_VCF="${REAL_DIR}/freebayes_tiny.vcf.gz"
FB_REF="${REAL_DIR}/freebayes_tiny_ref.fa"

if [[ "$FORCE" == true ]] || [[ ! -f "${FB_VCF}" ]]; then
	echo "  Downloading freebayes tiny VCF..." >&2
	wget -q "${FB_RAW}/q.vcf.gz" -O "${FB_VCF}"
	wget -q "${FB_RAW}/q.vcf.gz.tbi" -O "${FB_VCF}.tbi"
	echo "  Created: ${FB_VCF}" >&2
else
	echo "  Skipping freebayes_tiny.vcf.gz (exists; use --force to re-download)" >&2
fi

if [[ "$FORCE" == true ]] || [[ ! -f "${FB_REF}" ]]; then
	echo "  Downloading freebayes tiny reference FASTA..." >&2
	wget -q "${FB_RAW}/q.fa" -O "${FB_REF}"
	wget -q "${FB_RAW}/q.fa.fai" -O "${FB_REF}.fai"
	echo "  Created: ${FB_REF}" >&2
else
	echo "  Skipping freebayes_tiny_ref.fa (exists; use --force to re-download)" >&2
fi

# ---------------------------------------------------------------------------
# chr22 reference subset — for real-data integration tests
# ---------------------------------------------------------------------------
echo "" >&2
echo "[4/4] chr22:16M-16.1M reference subset..." >&2

CHR22_REF="${REAL_DIR}/chr22_16M_ref.fa"

if [[ "$FORCE" == true ]] || [[ ! -f "${CHR22_REF}" ]]; then
	# Prefer a local hs37d5.fa if available (avoids large download)
	if [[ -f "ref/hs37d5.fa" ]]; then
		echo "  Using local ref/hs37d5.fa..." >&2
		HS37D5="ref/hs37d5.fa"
		samtools faidx "${HS37D5}" "22:16000000-16100000" \
			| sed '1s/.*/\>22/' >"${CHR22_REF}"
	elif [[ -f "${TMPDIR}/hs37d5.fa" ]]; then
		echo "  Using previously downloaded hs37d5.fa..." >&2
		samtools faidx "${TMPDIR}/hs37d5.fa" "22:16000000-16100000" \
			| sed '1s/.*/\>22/' >"${CHR22_REF}"
	else
		echo "  Downloading hs37d5 reference (~3.4GB compressed — this will take a while)..." >&2
		echo "  TIP: Copy ref/hs37d5.fa to avoid re-downloading in future runs." >&2
		HS37D5_URL="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz"
		wget -q --show-progress "${HS37D5_URL}" -O "${TMPDIR}/hs37d5.fa.gz"
		echo "  Decompressing..." >&2
		gunzip "${TMPDIR}/hs37d5.fa.gz"
		samtools faidx "${TMPDIR}/hs37d5.fa"
		samtools faidx "${TMPDIR}/hs37d5.fa" "22:16000000-16100000" \
			| sed '1s/.*/\>22/' >"${CHR22_REF}"
	fi
	samtools faidx "${CHR22_REF}"
	echo "  Created: ${CHR22_REF}" >&2
	echo "  Contig header: $(head -1 "${CHR22_REF}")" >&2
else
	echo "  Skipping chr22_16M_ref.fa (exists; use --force to re-download)" >&2
fi

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
echo "" >&2
echo "=== Test data generation complete ===" >&2
echo "" >&2
echo "Files in ${REAL_DIR}:" >&2
ls -lh "${REAL_DIR}/" >&2

echo "" >&2
echo "Quick validation:" >&2
for vcf in "${REAL_DIR}"/*.vcf.gz; do
	count=$(bcftools view -H "$vcf" 2>/dev/null | wc -l)
	echo "  $(basename "$vcf"): ${count} records" >&2
done
