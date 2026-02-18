#!/usr/bin/env bats
# tests/03-freebayes-filters.bats — Freebayes filter unit tests (TFWK-03)
#
# Verifies each of the 9 Freebayes filters tags exactly the variants it should,
# and leaves control variants untagged. Pipeline runs once per file via
# setup_file(); individual @test blocks query the output VCF.
#
# Synthetic VCF: tests/data/synthetic/freebayes_sample_A.vcf.gz
# Filters: defaults/freebayes_filters.txt (9 filters)
#
# Expected FILTER assignments (from generate_expected.sh baseline):
#   POS 10101 — PASS                                  (control: good SNP, all checks pass)
#   POS 10201 — lowQUAL;QUALperAO                    (QUAL=15 < 20; QUAL/AO=0.75 < 10)
#   POS 10301 — QUALperAO                            (QUAL/AO=50/8=6.25 < 10; QUAL=50 >= 20)
#   POS 10401 — strandBias                           (SAF=0)
#   POS 10501 — strandBias                           (SAR=0)
#   POS 10601 — readPosBias                          (RPR=1 <= 1 AND RPL=1 <= 1)
#   POS 20101 — IN_EXCLUDE_REGION;DPu10het           (DP=8 < 10, het; in exclude region)
#   POS 20201 — IN_EXCLUDE_REGION;DPu5hom            (DP=3 < 5, hom; in exclude region)
#   POS 30101 — NOT_IN_INCLUDE_REGION;VAFu02het;readPosBias  (outside include; VAF=0.1 het; RPR=1&&RPL=1)
#   POS 40101 — VAFo08het                            (VAF=18/20=0.9 > 0.8, het; in include 40000-60000)
#   POS 50101 — VAFu095hom                           (VAF=18/20=0.9 < 0.95, hom)
#   POS 60101 — NOT_IN_INCLUDE_REGION                (outside include 40000-60000; no other filter)

setup_file() {
	load 'test_helper/common'
	_common_setup
	_require_tools

	export OUTPUT_VCF="${BATS_FILE_TMPDIR}/freebayes_filtered.vcf.gz"

	"$HARDNORMLY" \
		-v "${SYNTH}/freebayes_sample_A.vcf.gz" \
		-f "${SYNTH}/mini_ref.fa" \
		-b "${TEST_DATA}/include_regions.bed" \
		-e "${TEST_DATA}/exclude_regions.bed" \
		--filters-file "${REPO_ROOT}/defaults/freebayes_filters.txt" \
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

@test "Freebayes control SNP at 10101 is PASS" {
	local filter
	filter=$(_get_filter "$OUTPUT_VCF" 22 10101)
	[[ "$filter" = "PASS" ]]
}

# ---------------------------------------------------------------------------
# lowQUAL — QUAL < 20
# ---------------------------------------------------------------------------

@test "lowQUAL triggered by QUAL=15 at 10201 (combined with QUALperAO)" {
	local filter
	filter=$(_get_filter "$OUTPUT_VCF" 22 10201)
	[[ "$filter" = "lowQUAL;QUALperAO" ]]
}

# ---------------------------------------------------------------------------
# QUALperAO — QUAL / INFO/AO < 10
# ---------------------------------------------------------------------------

@test "QUALperAO triggered by QUAL/AO=6.25 at 10301 (QUAL=50 not lowQUAL)" {
	local filter
	filter=$(_get_filter "$OUTPUT_VCF" 22 10301)
	[[ "$filter" = "QUALperAO" ]]
}

@test "10201 has both lowQUAL and QUALperAO in combined tag" {
	local filter
	filter=$(_get_filter "$OUTPUT_VCF" 22 10201)
	[[ "$filter" == *"lowQUAL"* ]] && [[ "$filter" == *"QUALperAO"* ]]
}

# ---------------------------------------------------------------------------
# strandBias — INFO/SAF==0 || INFO/SAR==0
# ---------------------------------------------------------------------------

@test "strandBias triggered by SAF=0 at 10401" {
	local filter
	filter=$(_get_filter "$OUTPUT_VCF" 22 10401)
	[[ "$filter" = "strandBias" ]]
}

@test "strandBias triggered by SAR=0 at 10501" {
	local filter
	filter=$(_get_filter "$OUTPUT_VCF" 22 10501)
	[[ "$filter" = "strandBias" ]]
}

# ---------------------------------------------------------------------------
# readPosBias — INFO/RPR<=1 && INFO/RPL<=1
# ---------------------------------------------------------------------------

@test "readPosBias triggered by RPR=1 and RPL=1 at 10601" {
	local filter
	filter=$(_get_filter "$OUTPUT_VCF" 22 10601)
	[[ "$filter" = "readPosBias" ]]
}

# ---------------------------------------------------------------------------
# DPu10het — low depth heterozygous filter (combined with IN_EXCLUDE_REGION)
# ---------------------------------------------------------------------------

@test "DPu10het triggered by DP=8 het at 20101 (also IN_EXCLUDE_REGION)" {
	local filter
	filter=$(_get_filter "$OUTPUT_VCF" 22 20101)
	[[ "$filter" = "IN_EXCLUDE_REGION;DPu10het" ]]
}

# ---------------------------------------------------------------------------
# DPu5hom — low depth homozygous filter (combined with IN_EXCLUDE_REGION)
# ---------------------------------------------------------------------------

@test "DPu5hom triggered by DP=3 hom at 20201 (also IN_EXCLUDE_REGION)" {
	local filter
	filter=$(_get_filter "$OUTPUT_VCF" 22 20201)
	[[ "$filter" = "IN_EXCLUDE_REGION;DPu5hom" ]]
}

# ---------------------------------------------------------------------------
# VAFu02het + readPosBias + NOT_IN_INCLUDE_REGION — multi-filter variant
# ---------------------------------------------------------------------------

@test "30101 has NOT_IN_INCLUDE_REGION, VAFu02het, and readPosBias combined" {
	local filter
	filter=$(_get_filter "$OUTPUT_VCF" 22 30101)
	[[ "$filter" = "NOT_IN_INCLUDE_REGION;VAFu02het;readPosBias" ]]
}

@test "30101 contains VAFu02het tag (VAF=2/20=0.1 < 0.2, het)" {
	local filter
	filter=$(_get_filter "$OUTPUT_VCF" 22 30101)
	[[ "$filter" == *"VAFu02het"* ]]
}

@test "30101 contains readPosBias tag (RPR=1 and RPL=1)" {
	local filter
	filter=$(_get_filter "$OUTPUT_VCF" 22 30101)
	[[ "$filter" == *"readPosBias"* ]]
}

# ---------------------------------------------------------------------------
# VAFo08het — high VAF heterozygous filter
# ---------------------------------------------------------------------------

@test "VAFo08het triggered by high VAF het at 40101 (in include region)" {
	local filter
	filter=$(_get_filter "$OUTPUT_VCF" 22 40101)
	[[ "$filter" = "VAFo08het" ]]
}

# ---------------------------------------------------------------------------
# VAFu095hom — low VAF homozygous filter
# ---------------------------------------------------------------------------

@test "VAFu095hom triggered by VAF=0.9 hom at 50101" {
	local filter
	filter=$(_get_filter "$OUTPUT_VCF" 22 50101)
	[[ "$filter" = "VAFu095hom" ]]
}

# ---------------------------------------------------------------------------
# Region-only filter — NOT_IN_INCLUDE_REGION with no other trigger
# ---------------------------------------------------------------------------

@test "60101 outside include region gets NOT_IN_INCLUDE_REGION only" {
	local filter
	filter=$(_get_filter "$OUTPUT_VCF" 22 60101)
	[[ "$filter" = "NOT_IN_INCLUDE_REGION" ]]
}

# ---------------------------------------------------------------------------
# Negative assertions — control must not be tagged with any filter name
# ---------------------------------------------------------------------------

@test "Freebayes control SNP at 10101 is not tagged with lowQUAL" {
	local filter
	filter=$(_get_filter "$OUTPUT_VCF" 22 10101)
	[[ "$filter" != *"lowQUAL"* ]]
}

@test "Freebayes control SNP at 10101 is not tagged with strandBias" {
	local filter
	filter=$(_get_filter "$OUTPUT_VCF" 22 10101)
	[[ "$filter" != *"strandBias"* ]]
}

@test "Freebayes control SNP at 10101 is not tagged with readPosBias" {
	local filter
	filter=$(_get_filter "$OUTPUT_VCF" 22 10101)
	[[ "$filter" != *"readPosBias"* ]]
}
