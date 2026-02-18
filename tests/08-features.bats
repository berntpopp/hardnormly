#!/usr/bin/env bats
# tests/08-features.bats — Functional tests for Phase 5 features
#
# FEAT-06: --caller flag end-to-end (resolves to filter file, produces correct FILTER tags)
# FEAT-08: --strip-annotations flag (removes INFO fields from output)
# FEAT-03: generate-inclusion-bed subcommand (produces valid merged BED output)
# FEAT-04: generate-exclusion-bed subcommand (produces valid merged BED output)
#
# All tests require bioinformatics tools (bcftools/bedtools/bgzip/tabix).

setup() {
	load 'test_helper/common'
	_common_setup
	_require_tools

	OUTPUT_VCF="${BATS_TEST_TMPDIR}/output.vcf.gz"
	export OUTPUT_VCF
}

# ===========================================================================
# FEAT-06 — --caller flag resolves to filter file and produces correct tags
# ===========================================================================

@test "FEAT-06: --caller gatk produces same FILTER tags as --filters-file gatk_filters.txt" {
	local out_caller="${BATS_TEST_TMPDIR}/out_caller.vcf.gz"
	local out_file="${BATS_TEST_TMPDIR}/out_file.vcf.gz"

	# Run with --caller gatk
	"$HARDNORMLY" \
		-v "${SYNTH}/gatk_sample_A.vcf.gz" \
		-f "${SYNTH}/mini_ref.fa" \
		-g "${TEST_DATA}/hg19_chr22.genome" \
		-b "${TEST_DATA}/include_regions.bed" \
		-e "${TEST_DATA}/exclude_regions.bed" \
		--caller gatk \
		-o "$out_caller"

	# Run with explicit --filters-file
	"$HARDNORMLY" \
		-v "${SYNTH}/gatk_sample_A.vcf.gz" \
		-f "${SYNTH}/mini_ref.fa" \
		-g "${TEST_DATA}/hg19_chr22.genome" \
		-b "${TEST_DATA}/include_regions.bed" \
		-e "${TEST_DATA}/exclude_regions.bed" \
		--filters-file "${REPO_ROOT}/defaults/gatk_filters.txt" \
		-o "$out_file"

	# Both must produce identical FILTER columns
	local filters_caller filters_file
	filters_caller=$(bcftools query -f '%CHROM\t%POS\t%FILTER\n' "$out_caller")
	filters_file=$(bcftools query -f '%CHROM\t%POS\t%FILTER\n' "$out_file")

	[[ "$filters_caller" = "$filters_file" ]]
}

@test "FEAT-06: --caller gatk output contains gatkSNPhard FILTER tag" {
	run "$HARDNORMLY" \
		-v "${SYNTH}/gatk_sample_A.vcf.gz" \
		-f "${SYNTH}/mini_ref.fa" \
		-g "${TEST_DATA}/hg19_chr22.genome" \
		--caller gatk \
		-o "${OUTPUT_VCF}"

	assert_success

	# At least one variant should have gatkSNPhard tag (AS_FS=65 at 10201)
	local tags
	tags=$(bcftools query -f '%FILTER\n' "${OUTPUT_VCF}")
	[[ "$tags" == *"gatkSNPhard"* ]]
}

@test "FEAT-06: --caller freebayes output contains lowQUAL FILTER tag" {
	run "$HARDNORMLY" \
		-v "${SYNTH}/freebayes_sample_A.vcf.gz" \
		-f "${SYNTH}/mini_ref.fa" \
		-g "${TEST_DATA}/hg19_chr22.genome" \
		--caller freebayes \
		-o "${OUTPUT_VCF}"

	assert_success

	# At least one variant should have lowQUAL tag (QUAL=15 at 10201)
	local tags
	tags=$(bcftools query -f '%FILTER\n' "${OUTPUT_VCF}")
	[[ "$tags" == *"lowQUAL"* ]]
}

# ===========================================================================
# FEAT-08 — --strip-annotations removes INFO fields from output
# ===========================================================================

@test "FEAT-08: --strip-annotations removes specified INFO field from output" {
	run "$HARDNORMLY" \
		-v "${SYNTH}/gatk_sample_A.vcf.gz" \
		-f "${SYNTH}/mini_ref.fa" \
		-g "${TEST_DATA}/hg19_chr22.genome" \
		--strip-annotations "INFO/AS_FS" \
		-o "${OUTPUT_VCF}"

	assert_success

	# AS_FS should NOT appear in the output header
	local header
	header=$(bcftools view -h "${OUTPUT_VCF}")
	[[ "$header" != *"ID=AS_FS"* ]]
}

@test "FEAT-08: --strip-annotations preserves other INFO fields" {
	run "$HARDNORMLY" \
		-v "${SYNTH}/gatk_sample_A.vcf.gz" \
		-f "${SYNTH}/mini_ref.fa" \
		-g "${TEST_DATA}/hg19_chr22.genome" \
		--strip-annotations "INFO/AS_FS" \
		-o "${OUTPUT_VCF}"

	assert_success

	# AS_SOR should still be in the output header (not stripped)
	bcftools view -h "${OUTPUT_VCF}" | grep -q "ID=AS_SOR"
}

@test "FEAT-08: --strip-annotations with multiple fields removes all of them" {
	run "$HARDNORMLY" \
		-v "${SYNTH}/gatk_sample_A.vcf.gz" \
		-f "${SYNTH}/mini_ref.fa" \
		-g "${TEST_DATA}/hg19_chr22.genome" \
		--strip-annotations "INFO/AS_FS,INFO/AS_SOR" \
		-o "${OUTPUT_VCF}"

	assert_success

	local header
	header=$(bcftools view -h "${OUTPUT_VCF}")
	[[ "$header" != *"ID=AS_FS"* ]]
	[[ "$header" != *"ID=AS_SOR"* ]]
}

# ===========================================================================
# FEAT-03 — generate-inclusion-bed produces valid merged BED output
# ===========================================================================

@test "FEAT-03: generate-inclusion-bed produces non-empty BED file" {
	local output="${BATS_TEST_TMPDIR}/inclusion.bed"

	run "$HARDNORMLY" generate-inclusion-bed \
		-b "${TEST_DATA}/include_regions.bed" \
		-g "${TEST_DATA}/hg19_chr22.genome" \
		-o "$output"

	assert_success
	[[ -f "$output" ]]

	local lines
	lines=$(wc -l <"$output")
	[[ "$lines" -gt 0 ]]
}

@test "FEAT-03: generate-inclusion-bed output has 3+ columns and sorted coordinates" {
	local output="${BATS_TEST_TMPDIR}/inclusion.bed"

	"$HARDNORMLY" generate-inclusion-bed \
		-b "${TEST_DATA}/include_regions.bed" \
		-g "${TEST_DATA}/hg19_chr22.genome" \
		-o "$output"

	# Each line must have at least 3 tab-separated columns
	local min_cols
	min_cols=$(awk -F'\t' '{print NF}' "$output" | sort -n | head -1)
	[[ "$min_cols" -ge 3 ]]

	# Coordinates must be sorted (start of each line <= start of next)
	local sorted
	sorted=$(awk -F'\t' '{print $2}' "$output" | sort -n)
	local actual
	actual=$(awk -F'\t' '{print $2}' "$output")
	[[ "$sorted" = "$actual" ]]
}

@test "FEAT-03: generate-inclusion-bed applies slop to coordinates" {
	local output="${BATS_TEST_TMPDIR}/inclusion.bed"

	"$HARDNORMLY" generate-inclusion-bed \
		-b "${TEST_DATA}/include_regions.bed" \
		-g "${TEST_DATA}/hg19_chr22.genome" \
		--slop 20 \
		-o "$output"

	# include_regions.bed has 22:10000-30000 as first region
	# With slop=20, output should have 22:9980-30020
	local first_start first_end
	first_start=$(awk -F'\t' 'NR==1{print $2}' "$output")
	first_end=$(awk -F'\t' 'NR==1{print $3}' "$output")
	[[ "$first_start" -eq 9980 ]]
	[[ "$first_end" -eq 30020 ]]
}

# ===========================================================================
# FEAT-04 — generate-exclusion-bed produces valid merged BED output
# ===========================================================================

@test "FEAT-04: generate-exclusion-bed produces non-empty BED file" {
	local output="${BATS_TEST_TMPDIR}/exclusion.bed"

	run "$HARDNORMLY" generate-exclusion-bed \
		-e "${TEST_DATA}/exclude_regions.bed" \
		-o "$output"

	assert_success
	[[ -f "$output" ]]

	local lines
	lines=$(wc -l <"$output")
	[[ "$lines" -gt 0 ]]
}

@test "FEAT-04: generate-exclusion-bed output preserves region coordinates" {
	local output="${BATS_TEST_TMPDIR}/exclusion.bed"

	"$HARDNORMLY" generate-exclusion-bed \
		-e "${TEST_DATA}/exclude_regions.bed" \
		-o "$output"

	# exclude_regions.bed has 22:20000-25000 and 22:80000-85000
	local first_start first_end
	first_start=$(awk -F'\t' 'NR==1{print $2}' "$output")
	first_end=$(awk -F'\t' 'NR==1{print $3}' "$output")
	[[ "$first_start" -eq 20000 ]]
	[[ "$first_end" -eq 25000 ]]
}

@test "FEAT-04: generate-exclusion-bed merges multiple BED files" {
	local bed_a="${BATS_TEST_TMPDIR}/excl_a.bed"
	local bed_b="${BATS_TEST_TMPDIR}/excl_b.bed"
	local output="${BATS_TEST_TMPDIR}/merged_excl.bed"

	printf '22\t10000\t20000\n' >"$bed_a"
	printf '22\t30000\t40000\n' >"$bed_b"

	run "$HARDNORMLY" generate-exclusion-bed \
		-e "$bed_a" \
		-e "$bed_b" \
		-o "$output"

	assert_success

	# Merged output should contain regions from both files
	local lines
	lines=$(wc -l <"$output")
	[[ "$lines" -ge 2 ]]
}
