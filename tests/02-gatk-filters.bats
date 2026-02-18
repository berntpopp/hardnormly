#!/usr/bin/env bats
# tests/02-gatk-filters.bats — GATK filter unit tests (TFWK-02)
#
# Verifies each of the 7 GATK filters tags exactly the variants it should,
# and leaves control variants untagged. Pipeline runs once per file via
# setup_file(); individual @test blocks query the output VCF.
#
# Synthetic VCF: tests/data/synthetic/gatk_sample_A.vcf.gz
# Filters: defaults/gatk_filters.txt (7 filters)
#
# Expected FILTER assignments (from generate_expected.sh baseline):
#   POS 10101 — PASS           (control: good SNP, all GATK checks pass)
#   POS 10201 — gatkSNPhard    (AS_FS=65 > 60)
#   POS 10301 — gatkSNPhard    (AS_SOR=4 > 3)
#   POS 10401 — gatkSNPhard    (AS_MQ=35 < 40)
#   POS 10501 — gatkSNPhard    (AS_ReadPosRankSum=-10 < -8)
#   POS 10601 — gatkSNPhard    (AS_MQRankSum=-15 < -12.5)
#   POS 10701 — gatkSNPhard    (QUAL=20 < 30)
#   POS 20101 — IN_EXCLUDE_REGION                     (good INDEL in exclude region)
#   POS 20201 — IN_EXCLUDE_REGION;gatkINDELhard       (AS_FS=205 > 200)
#   POS 20301 — IN_EXCLUDE_REGION;gatkINDELhard       (AS_ReadPosRankSum=-22 < -20)
#   POS 20401 — IN_EXCLUDE_REGION;gatkINDELhard       (QUAL=25 < 30)
#   POS 30101 — NOT_IN_INCLUDE_REGION;DPu10het        (DP=8 < 10, het)
#   POS 30201 — NOT_IN_INCLUDE_REGION;DPu5hom         (DP=3 < 5, hom)
#   POS 40101 — VAFu02het                             (VAF=2/20=0.1 < 0.2, het; in include 40000-60000)
#   POS 50101 — VAFo08het                             (VAF=18/20=0.9 > 0.8, het)
#   POS 60101 — NOT_IN_INCLUDE_REGION;VAFu095hom      (outside include 40000-60000, hom VAF=18/20=0.9 < 0.95)

setup_file() {
	load 'test_helper/common'
	_common_setup
	_require_tools

	export OUTPUT_VCF="${BATS_FILE_TMPDIR}/gatk_filtered.vcf.gz"

	"$HARDNORMLY" \
		-v "${SYNTH}/gatk_sample_A.vcf.gz" \
		-f "${SYNTH}/mini_ref.fa" \
		-b "${TEST_DATA}/include_regions.bed" \
		-e "${TEST_DATA}/exclude_regions.bed" \
		--filters-file "${REPO_ROOT}/defaults/gatk_filters.txt" \
		-g "${TEST_DATA}/hg19_chr22.genome" \
		-o "$OUTPUT_VCF" \
		--auto-index
}

setup() {
	load 'test_helper/common'
	_common_setup
}

# ---------------------------------------------------------------------------
# Control variant — must PASS all filters
# ---------------------------------------------------------------------------

@test "GATK control SNP at 10101 is PASS" {
	local filter
	filter=$(_get_filter "$OUTPUT_VCF" 22 10101)
	[[ "$filter" = "PASS" ]]
}

# ---------------------------------------------------------------------------
# gatkSNPhard — individual trigger conditions
# ---------------------------------------------------------------------------

@test "gatkSNPhard triggered by AS_FS=65 at 10201" {
	local filter
	filter=$(_get_filter "$OUTPUT_VCF" 22 10201)
	[[ "$filter" = "gatkSNPhard" ]]
}

@test "gatkSNPhard triggered by AS_SOR=4 at 10301" {
	local filter
	filter=$(_get_filter "$OUTPUT_VCF" 22 10301)
	[[ "$filter" = "gatkSNPhard" ]]
}

@test "gatkSNPhard triggered by AS_MQ=35 at 10401" {
	local filter
	filter=$(_get_filter "$OUTPUT_VCF" 22 10401)
	[[ "$filter" = "gatkSNPhard" ]]
}

@test "gatkSNPhard triggered by AS_ReadPosRankSum=-10 at 10501" {
	local filter
	filter=$(_get_filter "$OUTPUT_VCF" 22 10501)
	[[ "$filter" = "gatkSNPhard" ]]
}

@test "gatkSNPhard triggered by AS_MQRankSum=-15 at 10601" {
	local filter
	filter=$(_get_filter "$OUTPUT_VCF" 22 10601)
	[[ "$filter" = "gatkSNPhard" ]]
}

@test "gatkSNPhard triggered by QUAL=20 at 10701" {
	local filter
	filter=$(_get_filter "$OUTPUT_VCF" 22 10701)
	[[ "$filter" = "gatkSNPhard" ]]
}

# ---------------------------------------------------------------------------
# Region-based filters — INDEL in exclude region
# ---------------------------------------------------------------------------

@test "INDEL in exclude region at 20101 gets IN_EXCLUDE_REGION only (no hard filter)" {
	local filter
	filter=$(_get_filter "$OUTPUT_VCF" 22 20101)
	[[ "$filter" = "IN_EXCLUDE_REGION" ]]
}

# ---------------------------------------------------------------------------
# gatkINDELhard — combined with IN_EXCLUDE_REGION (variants are in exclude region)
# ---------------------------------------------------------------------------

@test "gatkINDELhard triggered by AS_FS=205 at 20201 (also IN_EXCLUDE_REGION)" {
	local filter
	filter=$(_get_filter "$OUTPUT_VCF" 22 20201)
	[[ "$filter" = "IN_EXCLUDE_REGION;gatkINDELhard" ]]
}

@test "gatkINDELhard triggered by AS_ReadPosRankSum=-22 at 20301 (also IN_EXCLUDE_REGION)" {
	local filter
	filter=$(_get_filter "$OUTPUT_VCF" 22 20301)
	[[ "$filter" = "IN_EXCLUDE_REGION;gatkINDELhard" ]]
}

@test "gatkINDELhard triggered by QUAL=25 at 20401 (also IN_EXCLUDE_REGION)" {
	local filter
	filter=$(_get_filter "$OUTPUT_VCF" 22 20401)
	[[ "$filter" = "IN_EXCLUDE_REGION;gatkINDELhard" ]]
}

# ---------------------------------------------------------------------------
# DPu10het — low depth heterozygous filter
# ---------------------------------------------------------------------------

@test "DPu10het triggered by DP=8 het at 30101 (also NOT_IN_INCLUDE_REGION)" {
	local filter
	filter=$(_get_filter "$OUTPUT_VCF" 22 30101)
	[[ "$filter" = "NOT_IN_INCLUDE_REGION;DPu10het" ]]
}

# ---------------------------------------------------------------------------
# DPu5hom — low depth homozygous filter
# ---------------------------------------------------------------------------

@test "DPu5hom triggered by DP=3 hom at 30201 (also NOT_IN_INCLUDE_REGION)" {
	local filter
	filter=$(_get_filter "$OUTPUT_VCF" 22 30201)
	[[ "$filter" = "NOT_IN_INCLUDE_REGION;DPu5hom" ]]
}

# ---------------------------------------------------------------------------
# VAF filters
# ---------------------------------------------------------------------------

@test "VAFu02het triggered by low VAF het at 40101" {
	local filter
	filter=$(_get_filter "$OUTPUT_VCF" 22 40101)
	[[ "$filter" = "VAFu02het" ]]
}

@test "VAFo08het triggered by high VAF het at 50101" {
	local filter
	filter=$(_get_filter "$OUTPUT_VCF" 22 50101)
	[[ "$filter" = "VAFo08het" ]]
}

@test "VAFu095hom triggered by low VAF hom at 60101 (also NOT_IN_INCLUDE_REGION)" {
	local filter
	filter=$(_get_filter "$OUTPUT_VCF" 22 60101)
	[[ "$filter" = "NOT_IN_INCLUDE_REGION;VAFu095hom" ]]
}

# ---------------------------------------------------------------------------
# Negative assertions — control must not be tagged with any filter name
# ---------------------------------------------------------------------------

@test "GATK control SNP at 10101 is not tagged with gatkSNPhard" {
	local filter
	filter=$(_get_filter "$OUTPUT_VCF" 22 10101)
	[[ "$filter" != *"gatkSNPhard"* ]]
}

@test "GATK control SNP at 10101 is not tagged with DPu10het" {
	local filter
	filter=$(_get_filter "$OUTPUT_VCF" 22 10101)
	[[ "$filter" != *"DPu10het"* ]]
}
