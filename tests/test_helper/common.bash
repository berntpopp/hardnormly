# tests/test_helper/common.bash
# Shared BATS helper sourced by every test file via _common_setup()

_common_setup() {
	load 'test_helper/bats-support/load'
	load 'test_helper/bats-assert/load'
	load 'test_helper/bats-file/load'

	# Compute repo root from test file location so tests run from any directory
	REPO_ROOT="$(cd "$(dirname "$BATS_TEST_FILENAME")/.." >/dev/null 2>&1 && pwd)"
	export REPO_ROOT

	# Canonical paths to key test assets
	HARDNORMLY="${REPO_ROOT}/hardnormly.sh"
	TEST_DATA="${REPO_ROOT}/tests/data"
	SYNTH="${TEST_DATA}/synthetic"
	REAL_DATA="${TEST_DATA}/real"
	EXPECTED="${TEST_DATA}/expected"
	export HARDNORMLY TEST_DATA SYNTH REAL_DATA EXPECTED

	# Per-test temp directory provided by BATS
	TEST_TMP="${BATS_TEST_TMPDIR}"
	export TEST_TMP
}

# Skip test if bioinformatics tools are not available
_require_tools() {
	local missing=()
	for tool in bcftools bedtools bgzip tabix; do
		command -v "$tool" >/dev/null 2>&1 || missing+=("$tool")
	done
	if [[ ${#missing[@]} -gt 0 ]]; then
		skip "Required tools not found: ${missing[*]}"
	fi
}

# Extract FILTER field for a specific chromosome position from a VCF
# Usage: _get_filter <vcf> <chrom> <pos>
# For chr22 VCFs: _get_filter "$vcf" 22 10200
_get_filter() {
	local vcf="$1"
	local chrom="$2"
	local pos="$3"
	bcftools query -r "${chrom}:${pos}-${pos}" -f '%FILTER\n' "$vcf"
}
