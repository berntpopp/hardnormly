#!/usr/bin/env bats
# tests/04-integration.bats — Integration, regression, genome flag, and edge case tests
#
# TFWK-04: Full pipeline integration on real 1000G data
# TFWK-05: Regression test against committed golden file
# TFWK-06: -g/--genome flag verification (issue #13)
# TFWK-07: Edge cases — empty VCF, no filters, no BED
#
# Each test runs a complete independent pipeline invocation (no setup_file).
# The full reference (ref/hs37d5.fa) is used for real data tests because the
# 1000G variants are at chr22:16M+ and the synthetic mini_ref.fa only covers
# chr22:1-100001.
#
# Real data:  tests/data/real/1kg_NA12878_chr22_16M.vcf.gz  (20 variants, GT only)
# Regression: tests/data/synthetic/gatk_sample_A.vcf.gz vs
#             tests/data/expected/gatk_sample_A_filtered.vcf.gz
# Edge cases: tests/data/synthetic/empty_variants.vcf.gz

setup() {
	load 'test_helper/common'
	_common_setup
	_require_tools

	# Full reference genome covering chr22 completely (required for 16M+ variants)
	FULL_REF="${REPO_ROOT}/ref/hs37d5.fa"
	export FULL_REF

	# Convenience: output VCF in per-test temp dir
	OUTPUT_VCF="${BATS_TEST_TMPDIR}/output.vcf.gz"
	export OUTPUT_VCF

	# Skip real-data tests if full reference is absent (e.g. CI without ref/)
	if [[ ! -f "${FULL_REF}" ]]; then
		HAVE_FULL_REF=false
	else
		HAVE_FULL_REF=true
	fi
	export HAVE_FULL_REF
}

# ===========================================================================
# TFWK-04 — Integration: full pipeline on real 1000G data
# ===========================================================================

@test "TFWK-04: full pipeline on real 1000G NA12878 data produces non-empty output" {
	if [[ "${HAVE_FULL_REF}" == "false" ]]; then
		skip "Full reference ref/hs37d5.fa not available"
	fi

	run "$HARDNORMLY" \
		-v "${REAL_DATA}/1kg_NA12878_chr22_16M.vcf.gz" \
		-f "${FULL_REF}" \
		-g "${TEST_DATA}/hg19_chr22.genome" \
		-o "${OUTPUT_VCF}"

	assert_success

	# Output VCF must exist and contain at least one variant
	local count
	count=$(bcftools view -H "${OUTPUT_VCF}" | wc -l)
	[[ "$count" -gt 0 ]]
}

@test "TFWK-04: real data output is valid bgzip-compressed VCF" {
	if [[ "${HAVE_FULL_REF}" == "false" ]]; then
		skip "Full reference ref/hs37d5.fa not available"
	fi

	run "$HARDNORMLY" \
		-v "${REAL_DATA}/1kg_NA12878_chr22_16M.vcf.gz" \
		-f "${FULL_REF}" \
		-g "${TEST_DATA}/hg19_chr22.genome" \
		-o "${OUTPUT_VCF}"

	assert_success

	# bcftools view -h must succeed on a valid compressed VCF
	run bcftools view -h "${OUTPUT_VCF}"
	assert_success
	assert_output --partial "#CHROM"
}

@test "TFWK-04: pipeline with GIAB high-confidence BED includes variants in region" {
	if [[ "${HAVE_FULL_REF}" == "false" ]]; then
		skip "Full reference ref/hs37d5.fa not available"
	fi

	run "$HARDNORMLY" \
		-v "${REAL_DATA}/giab_NA12878_chr22_16M.vcf.gz" \
		-f "${FULL_REF}" \
		-g "${TEST_DATA}/hg19_chr22.genome" \
		-b "${REAL_DATA}/giab_NA12878_chr22_16M_highconf.bed" \
		-o "${OUTPUT_VCF}"

	assert_success

	# Output must have at least one variant (GIAB VCF has 7 variants, all in highconf BED)
	local count
	count=$(bcftools view -H "${OUTPUT_VCF}" | wc -l)
	[[ "$count" -gt 0 ]]
}

# ===========================================================================
# TFWK-05 — Regression: compare against committed golden file
# ===========================================================================

@test "TFWK-05: gatk_sample_A filter results match committed golden file" {
	run "$HARDNORMLY" \
		-v "${SYNTH}/gatk_sample_A.vcf.gz" \
		-f "${SYNTH}/mini_ref.fa" \
		-b "${TEST_DATA}/include_regions.bed" \
		-e "${TEST_DATA}/exclude_regions.bed" \
		--filters-file "${REPO_ROOT}/defaults/gatk_filters.txt" \
		-g "${TEST_DATA}/hg19_chr22.genome" \
		-o "${OUTPUT_VCF}"

	assert_success

	# Extract CHROM/POS/FILTER from actual and expected for comparison
	local actual_filters expected_filters
	actual_filters=$(bcftools query -f '%CHROM\t%POS\t%FILTER\n' "${OUTPUT_VCF}")
	expected_filters=$(bcftools query -f '%CHROM\t%POS\t%FILTER\n' "${EXPECTED}/gatk_sample_A_filtered.vcf.gz")

	if [[ "$actual_filters" != "$expected_filters" ]]; then
		local diff_count first_diffs
		diff_count=$(diff <(echo "$expected_filters") <(echo "$actual_filters") | grep -c '^[<>]')
		first_diffs=$(diff <(echo "$expected_filters") <(echo "$actual_filters") | head -20)
		echo "Regression failure: ${diff_count} changed lines" >&2
		echo "First differences (expected=< actual=>):" >&2
		echo "$first_diffs" >&2
		return 1
	fi
}

@test "TFWK-05: regression output has same variant count as golden file" {
	run "$HARDNORMLY" \
		-v "${SYNTH}/gatk_sample_A.vcf.gz" \
		-f "${SYNTH}/mini_ref.fa" \
		-b "${TEST_DATA}/include_regions.bed" \
		-e "${TEST_DATA}/exclude_regions.bed" \
		--filters-file "${REPO_ROOT}/defaults/gatk_filters.txt" \
		-g "${TEST_DATA}/hg19_chr22.genome" \
		-o "${OUTPUT_VCF}"

	assert_success

	local actual_count expected_count
	actual_count=$(bcftools view -H "${OUTPUT_VCF}" | wc -l)
	expected_count=$(bcftools view -H "${EXPECTED}/gatk_sample_A_filtered.vcf.gz" | wc -l)

	[[ "$actual_count" -eq "$expected_count" ]]
}

# ===========================================================================
# TFWK-06 — Genome flag: verify -g and --genome work (issue #13)
# ===========================================================================

@test "TFWK-06: -g short flag provides genome file and skips UCSC query" {
	run "$HARDNORMLY" \
		-v "${SYNTH}/gatk_sample_A.vcf.gz" \
		-f "${SYNTH}/mini_ref.fa" \
		-b "${TEST_DATA}/include_regions.bed" \
		-g "${TEST_DATA}/hg19_chr22.genome" \
		--filters-file "${REPO_ROOT}/defaults/gatk_filters.txt" \
		-o "${OUTPUT_VCF}"

	assert_success
	# Log must confirm the provided genome file was used (not generated)
	assert_output --partial "Using provided genome file"
}

@test "TFWK-06: --genome long flag provides genome file and skips UCSC query" {
	run "$HARDNORMLY" \
		-v "${SYNTH}/gatk_sample_A.vcf.gz" \
		-f "${SYNTH}/mini_ref.fa" \
		-b "${TEST_DATA}/include_regions.bed" \
		--genome "${TEST_DATA}/hg19_chr22.genome" \
		--filters-file "${REPO_ROOT}/defaults/gatk_filters.txt" \
		-o "${OUTPUT_VCF}"

	assert_success
	assert_output --partial "Using provided genome file"
}

@test "TFWK-06: -g and --genome produce identical output" {
	local out_short="${BATS_TEST_TMPDIR}/out_short.vcf.gz"
	local out_long="${BATS_TEST_TMPDIR}/out_long.vcf.gz"

	"$HARDNORMLY" \
		-v "${SYNTH}/gatk_sample_A.vcf.gz" \
		-f "${SYNTH}/mini_ref.fa" \
		-b "${TEST_DATA}/include_regions.bed" \
		-g "${TEST_DATA}/hg19_chr22.genome" \
		--filters-file "${REPO_ROOT}/defaults/gatk_filters.txt" \
		-o "$out_short"

	"$HARDNORMLY" \
		-v "${SYNTH}/gatk_sample_A.vcf.gz" \
		-f "${SYNTH}/mini_ref.fa" \
		-b "${TEST_DATA}/include_regions.bed" \
		--genome "${TEST_DATA}/hg19_chr22.genome" \
		--filters-file "${REPO_ROOT}/defaults/gatk_filters.txt" \
		-o "$out_long"

	# Both must produce the same CHROM/POS/FILTER results
	local filters_short filters_long
	filters_short=$(bcftools query -f '%CHROM\t%POS\t%FILTER\n' "$out_short")
	filters_long=$(bcftools query -f '%CHROM\t%POS\t%FILTER\n' "$out_long")

	[[ "$filters_short" = "$filters_long" ]]
}

# ===========================================================================
# TFWK-07 — Edge cases
# ===========================================================================

@test "TFWK-07: empty VCF input does not crash the pipeline" {
	run "$HARDNORMLY" \
		-v "${SYNTH}/empty_variants.vcf.gz" \
		-f "${SYNTH}/mini_ref.fa" \
		-g "${TEST_DATA}/hg19_chr22.genome" \
		-o "${OUTPUT_VCF}"

	assert_success

	# Output must be a valid VCF (header present)
	run bcftools view -h "${OUTPUT_VCF}"
	assert_success
	assert_output --partial "#CHROM"
}

@test "TFWK-07: empty VCF input produces output with zero variants" {
	run "$HARDNORMLY" \
		-v "${SYNTH}/empty_variants.vcf.gz" \
		-f "${SYNTH}/mini_ref.fa" \
		-g "${TEST_DATA}/hg19_chr22.genome" \
		-o "${OUTPUT_VCF}"

	assert_success

	local count
	count=$(bcftools view -H "${OUTPUT_VCF}" | wc -l)
	[[ "$count" -eq 0 ]]
}

@test "TFWK-07: running without filters file does not crash the pipeline" {
	# No --filters-file and no --filters provided
	run "$HARDNORMLY" \
		-v "${SYNTH}/gatk_sample_A.vcf.gz" \
		-f "${SYNTH}/mini_ref.fa" \
		-g "${TEST_DATA}/hg19_chr22.genome" \
		-b "${TEST_DATA}/include_regions.bed" \
		-e "${TEST_DATA}/exclude_regions.bed" \
		-o "${OUTPUT_VCF}"

	assert_success

	# Output must still have variants (region filters apply but no hard filters)
	local count
	count=$(bcftools view -H "${OUTPUT_VCF}" | wc -l)
	[[ "$count" -gt 0 ]]
}

@test "TFWK-07: running without BED files does not crash the pipeline" {
	# No -b and no -e provided
	run "$HARDNORMLY" \
		-v "${SYNTH}/gatk_sample_A.vcf.gz" \
		-f "${SYNTH}/mini_ref.fa" \
		-g "${TEST_DATA}/hg19_chr22.genome" \
		--filters-file "${REPO_ROOT}/defaults/gatk_filters.txt" \
		-o "${OUTPUT_VCF}"

	assert_success

	# Output must still have variants (only hard filters apply, no region filters)
	local count
	count=$(bcftools view -H "${OUTPUT_VCF}" | wc -l)
	[[ "$count" -gt 0 ]]
}

@test "TFWK-07: running without BED files or filters succeeds and preserves all variants" {
	# Bare minimum: only VCF + FASTA + genome + output
	run "$HARDNORMLY" \
		-v "${SYNTH}/gatk_sample_A.vcf.gz" \
		-f "${SYNTH}/mini_ref.fa" \
		-g "${TEST_DATA}/hg19_chr22.genome" \
		-o "${OUTPUT_VCF}"

	assert_success

	# All 16 variants must be present (no filtering applied)
	local count
	count=$(bcftools view -H "${OUTPUT_VCF}" | wc -l)
	[[ "$count" -eq 16 ]]
}
