#!/usr/bin/env bats
# tests/07-lib-modules.bats — Unit tests for pipeline lib/ modules (REFR-13)
#
# Covers: bed.sh, annotate.sh, normalize.sh, filter.sh, stats.sh, genome.sh
# All tests use bioinformatics tools (bcftools/bedtools/bgzip/tabix).
# Sources all 8 lib/ modules in setup().
#
# Test data:
#   tests/data/synthetic/minimal.vcf.gz     — single variant chr22:10001 (GT only)
#   tests/data/synthetic/gatk_sample_A.vcf.gz — 16 GATK variants with INFO fields
#   tests/data/synthetic/mini_ref.fa        — all-N chr22:1-100001 reference
#   tests/data/include_regions.bed          — chr22 include regions
#   tests/data/hg19_chr22.genome            — chr22 genome file

bats_require_minimum_version 1.5.0

setup() {
	load 'test_helper/common'
	_common_setup
	_require_tools
	source "${REPO_ROOT}/lib/logging.sh"
	source "${REPO_ROOT}/lib/cli.sh"
	source "${REPO_ROOT}/lib/genome.sh"
	source "${REPO_ROOT}/lib/bed.sh"
	source "${REPO_ROOT}/lib/annotate.sh"
	source "${REPO_ROOT}/lib/normalize.sh"
	source "${REPO_ROOT}/lib/filter.sh"
	source "${REPO_ROOT}/lib/stats.sh"
	set_log_file ""
	set_debug "false"
	set_tmp_dir "${BATS_TEST_TMPDIR}"
}

# ===========================================================================
# lib/bed.sh — normalize_bed, create_header_file, compress_index_bed
# ===========================================================================

@test "normalize_bed adds annotation column to BED file" {
	local input="${BATS_TEST_TMPDIR}/input.bed"
	local output="${BATS_TEST_TMPDIR}/output.bed"
	printf '22\t10000\t20000\n' >"$input"
	normalize_bed "$input" "1" "$output"
	local cols
	cols=$(awk '{print NF}' "$output")
	# Must have 4 columns (chrom, start, end, annotation)
	[ "$cols" -eq 4 ]
}

@test "normalize_bed annotation column contains expected value" {
	local input="${BATS_TEST_TMPDIR}/input.bed"
	local output="${BATS_TEST_TMPDIR}/output.bed"
	printf '22\t10000\t20000\n' >"$input"
	normalize_bed "$input" "MYREGION" "$output"
	local annot
	annot=$(awk '{print $4}' "$output")
	[ "$annot" = "MYREGION" ]
}

@test "normalize_bed sorts multiple regions by coordinate" {
	local input="${BATS_TEST_TMPDIR}/input.bed"
	local output="${BATS_TEST_TMPDIR}/output.bed"
	printf '22\t50000\t60000\n22\t10000\t20000\n' >"$input"
	normalize_bed "$input" "1" "$output"
	local first_start
	first_start=$(awk 'NR==1{print $2}' "$output")
	# First region after sort must be the one with lower start coordinate
	[ "$first_start" -eq 10000 ]
}

@test "create_header_file writes VCF INFO header line" {
	local hdr="${BATS_TEST_TMPDIR}/test.hdr"
	create_header_file "MY_FIELD" "Test description" "$hdr"
	[ -f "$hdr" ]
	grep -q "##INFO=" "$hdr"
}

@test "create_header_file header contains correct field ID" {
	local hdr="${BATS_TEST_TMPDIR}/test.hdr"
	create_header_file "INCLUDE_REGION" "Inclusion region flag" "$hdr"
	grep -q "ID=INCLUDE_REGION" "$hdr"
}

@test "create_header_file header has Type=Integer" {
	local hdr="${BATS_TEST_TMPDIR}/test.hdr"
	create_header_file "MYFIELD" "desc" "$hdr"
	grep -q "Type=Integer" "$hdr"
}

@test "compress_index_bed creates .gz and .tbi files" {
	local bed="${BATS_TEST_TMPDIR}/regions.bed"
	printf '22\t10000\t20000\t1\n' >"$bed"
	compress_index_bed "$bed"
	[ -f "${bed}.gz" ]
	[ -f "${bed}.gz.tbi" ]
}

# ===========================================================================
# lib/annotate.sh — annotate_vcf_with_regions
# ===========================================================================

@test "annotate_vcf_with_regions adds INFO field to VCF" {
	# Prepare a normalized + compressed include BED
	local bed="${BATS_TEST_TMPDIR}/include.bed"
	local hdr="${BATS_TEST_TMPDIR}/include.hdr"
	local output="${BATS_TEST_TMPDIR}/annotated.vcf.gz"

	# minimal.vcf.gz has a variant at 22:10001, which falls in include region 10000-30000
	printf '22\t10000\t30000\t1\n' >"$bed"
	bgzip -f "$bed"
	tabix -p bed "${bed}.gz"
	create_header_file "INCLUDE_REGION" "Inclusion region" "$hdr"

	annotate_vcf_with_regions \
		"${SYNTH}/minimal.vcf.gz" \
		"${bed}.gz" \
		"$hdr" \
		"INCLUDE_REGION" \
		"$output"

	# Header must declare the new INFO field
	bcftools view -h "$output" | grep -q "INCLUDE_REGION"
}

@test "annotate_vcf_with_regions sets INFO field to 1 for overlapping variant" {
	local bed="${BATS_TEST_TMPDIR}/include.bed"
	local hdr="${BATS_TEST_TMPDIR}/include.hdr"
	local output="${BATS_TEST_TMPDIR}/annotated.vcf.gz"

	# minimal.vcf.gz variant at 22:10001; region covers it
	printf '22\t10000\t30000\t1\n' >"$bed"
	bgzip -f "$bed"
	tabix -p bed "${bed}.gz"
	create_header_file "INCLUDE_REGION" "Inclusion region" "$hdr"

	annotate_vcf_with_regions \
		"${SYNTH}/minimal.vcf.gz" \
		"${bed}.gz" \
		"$hdr" \
		"INCLUDE_REGION" \
		"$output"

	local info
	info=$(bcftools query -f '%INFO/INCLUDE_REGION\n' "$output")
	[ "$info" = "1" ]
}

# ===========================================================================
# lib/normalize.sh — normalize_vcf
# ===========================================================================

@test "normalize_vcf produces output file" {
	local output="${BATS_TEST_TMPDIR}/normalized.vcf.gz"
	normalize_vcf \
		"${SYNTH}/minimal.vcf.gz" \
		"${SYNTH}/mini_ref.fa" \
		"$output" \
		"${BATS_TEST_TMPDIR}"
	[ -f "$output" ]
}

@test "normalize_vcf output is valid VCF with header" {
	local output="${BATS_TEST_TMPDIR}/normalized.vcf.gz"
	normalize_vcf \
		"${SYNTH}/minimal.vcf.gz" \
		"${SYNTH}/mini_ref.fa" \
		"$output" \
		"${BATS_TEST_TMPDIR}"
	run bcftools view -h "$output"
	assert_success
	assert_output --partial "#CHROM"
}

@test "normalize_vcf preserves variant count for already-normalized VCF" {
	local output="${BATS_TEST_TMPDIR}/normalized.vcf.gz"
	normalize_vcf \
		"${SYNTH}/minimal.vcf.gz" \
		"${SYNTH}/mini_ref.fa" \
		"$output" \
		"${BATS_TEST_TMPDIR}"
	local count
	count=$(bcftools view -H "$output" | wc -l)
	# minimal.vcf.gz has exactly 1 variant
	[ "$count" -eq 1 ]
}

@test "normalize_vcf splits multiallelic sites" {
	local output="${BATS_TEST_TMPDIR}/normalized.vcf.gz"
	normalize_vcf \
		"${SYNTH}/multiallelic.vcf.gz" \
		"${SYNTH}/mini_ref.fa" \
		"$output" \
		"${BATS_TEST_TMPDIR}"
	local count
	count=$(bcftools view -H "$output" | wc -l)
	# multiallelic site should be split into multiple records
	[ "$count" -gt 1 ]
}

# ===========================================================================
# lib/filter.sh — init_filter_pipeline, apply_filter_stages, write_filtered_output
# ===========================================================================

@test "init_filter_pipeline creates filter_current.bcf" {
	init_filter_pipeline "${SYNTH}/minimal.vcf.gz" "${BATS_TEST_TMPDIR}"
	[ -f "${BATS_TEST_TMPDIR}/filter_current.bcf" ]
}

@test "apply_filter_stages applies a single exclude filter correctly" {
	init_filter_pipeline "${SYNTH}/minimal.vcf.gz" "${BATS_TEST_TMPDIR}"
	# minimal.vcf.gz: DP=20, GT=0/1 — DPu10het triggers on DP<10; DP=20 should NOT trigger
	apply_filter_stages "${BATS_TEST_TMPDIR}" "DPtest|e|FORMAT/DP<5"
	# Variant at 22:10001 with DP=20 should remain PASS
	local filter
	filter=$(bcftools query -f '%FILTER\n' "${BATS_TEST_TMPDIR}/filter_current.bcf")
	[ "$filter" = "PASS" ]
}

@test "apply_filter_stages tags variant when filter expression matches" {
	init_filter_pipeline "${SYNTH}/minimal.vcf.gz" "${BATS_TEST_TMPDIR}"
	# minimal.vcf.gz DP=20, GT=0/1 — exclude on DP<30 should tag the variant
	apply_filter_stages "${BATS_TEST_TMPDIR}" "lowDP|e|FORMAT/DP<30"
	local filter
	filter=$(bcftools query -f '%FILTER\n' "${BATS_TEST_TMPDIR}/filter_current.bcf")
	[ "$filter" = "lowDP" ]
}

@test "apply_filter_stages chains multiple filters" {
	init_filter_pipeline "${SYNTH}/minimal.vcf.gz" "${BATS_TEST_TMPDIR}"
	# Two filters — both conditions apply
	apply_filter_stages "${BATS_TEST_TMPDIR}" \
		"filterA|e|FORMAT/DP<30" \
		"filterB|e|FORMAT/DP<25"
	local filter
	filter=$(bcftools query -f '%FILTER\n' "${BATS_TEST_TMPDIR}/filter_current.bcf")
	# Both filter tags should be present (semicolon-delimited)
	[[ "$filter" == *"filterA"* ]]
	[[ "$filter" == *"filterB"* ]]
}

@test "write_filtered_output produces compressed VCF file" {
	init_filter_pipeline "${SYNTH}/minimal.vcf.gz" "${BATS_TEST_TMPDIR}"
	local output="${BATS_TEST_TMPDIR}/output.vcf.gz"
	write_filtered_output "${BATS_TEST_TMPDIR}" "$output" "false" "false"
	[ -f "$output" ]
}

@test "write_filtered_output produces valid VCF with variants" {
	init_filter_pipeline "${SYNTH}/minimal.vcf.gz" "${BATS_TEST_TMPDIR}"
	local output="${BATS_TEST_TMPDIR}/output.vcf.gz"
	write_filtered_output "${BATS_TEST_TMPDIR}" "$output" "false" "false"
	local count
	count=$(bcftools view -H "$output" | wc -l)
	[ "$count" -eq 1 ]
}

@test "write_filtered_output applies PASS filter when only_pass is true" {
	init_filter_pipeline "${SYNTH}/minimal.vcf.gz" "${BATS_TEST_TMPDIR}"
	# Tag the variant with a filter
	apply_filter_stages "${BATS_TEST_TMPDIR}" "lowDP|e|FORMAT/DP<30"
	local output="${BATS_TEST_TMPDIR}/output.vcf.gz"
	write_filtered_output "${BATS_TEST_TMPDIR}" "$output" "true" "false"
	# Variant was tagged, so PASS filter removes it — output should be empty
	local count
	count=$(bcftools view -H "$output" | wc -l)
	[ "$count" -eq 0 ]
}

# ===========================================================================
# lib/stats.sh — generate_stats
# ===========================================================================

@test "generate_stats produces a stats file" {
	local stats="${BATS_TEST_TMPDIR}/output.stats"
	generate_stats "${SYNTH}/minimal.vcf.gz" "$stats"
	[ -f "$stats" ]
}

@test "generate_stats output contains SN section header" {
	local stats="${BATS_TEST_TMPDIR}/output.stats"
	generate_stats "${SYNTH}/minimal.vcf.gz" "$stats"
	grep -q "^SN" "$stats"
}

@test "generate_stats output reports at least one variant" {
	local stats="${BATS_TEST_TMPDIR}/output.stats"
	generate_stats "${SYNTH}/minimal.vcf.gz" "$stats"
	# SN section contains number of records
	local records
	records=$(grep "^SN.*number of records" "$stats" | awk '{print $NF}')
	[ "$records" -ge 1 ]
}

# ===========================================================================
# lib/genome.sh — create_genome_file (existence check only, no network call)
# ===========================================================================

@test "create_genome_file function is defined" {
	# Verify the function exists without making a network call
	declare -f create_genome_file >/dev/null
}
