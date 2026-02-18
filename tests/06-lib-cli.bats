#!/usr/bin/env bats
# tests/06-lib-cli.bats — Unit tests for lib/cli.sh (REFR-13)
#
# Verifies CLI functions in isolation. Sources lib/logging.sh then lib/cli.sh.
# No bioinformatics tools required.
#
# Functions covered: show_help, show_version, parse_filter_args

setup() {
	load 'test_helper/common'
	_common_setup
	source "${REPO_ROOT}/lib/logging.sh"
	set_log_file ""
	set_debug "false"
	set_tmp_dir "${BATS_TEST_TMPDIR}"
	source "${REPO_ROOT}/lib/cli.sh"
}

# ---------------------------------------------------------------------------
# show_help — outputs usage text and exits 1
# ---------------------------------------------------------------------------

@test "show_help outputs Usage: text" {
	run show_help
	assert_output --partial "Usage:"
}

@test "show_help outputs options section" {
	run show_help
	assert_output --partial "Options (run-pipeline):"
}

@test "show_help exits with non-zero status" {
	run show_help
	assert_failure
}

@test "show_help mentions required -v flag" {
	run show_help
	assert_output --partial "-v"
}

@test "show_help mentions required -f flag" {
	run show_help
	assert_output --partial "-f"
}

# ---------------------------------------------------------------------------
# show_version — outputs version string and exits 0
# ---------------------------------------------------------------------------

@test "show_version outputs Version: prefix" {
	run show_version "1.2.3"
	assert_output --partial "Version:"
}

@test "show_version includes the version number passed in" {
	run show_version "0.6.0"
	assert_output --partial "0.6.0"
}

@test "show_version exits successfully" {
	run show_version "1.0.0"
	assert_success
}

# ---------------------------------------------------------------------------
# parse_filter_args — parses inline and file-based filters into array
# ---------------------------------------------------------------------------

@test "parse_filter_args parses a single inline filter into name|action|expr" {
	local stages=()
	parse_filter_args stages "" "myFilter e FORMAT/DP<10"
	[ "${#stages[@]}" -eq 1 ]
	[ "${stages[0]}" = "myFilter|e|FORMAT/DP<10" ]
}

@test "parse_filter_args parses multiple inline filters" {
	local stages=()
	parse_filter_args stages "" \
		"filterA e FORMAT/DP<10" \
		"filterB i FORMAT/GQ>30"
	[ "${#stages[@]}" -eq 2 ]
	[ "${stages[0]}" = "filterA|e|FORMAT/DP<10" ]
	[ "${stages[1]}" = "filterB|i|FORMAT/GQ>30" ]
}

@test "parse_filter_args parses filters from a file" {
	local filter_file="${BATS_TEST_TMPDIR}/filters.txt"
	printf 'fileFilter e FORMAT/DP<5\n' >"$filter_file"
	local stages=()
	parse_filter_args stages "$filter_file"
	[ "${#stages[@]}" -eq 1 ]
	[ "${stages[0]}" = "fileFilter|e|FORMAT/DP<5" ]
}

@test "parse_filter_args combines inline and file-based filters" {
	local filter_file="${BATS_TEST_TMPDIR}/filters.txt"
	printf 'fileFilter e FORMAT/GQ>20\n' >"$filter_file"
	local stages=()
	parse_filter_args stages "$filter_file" "inlineFilter i FORMAT/DP>10"
	# Inline comes first, then file
	[ "${#stages[@]}" -eq 2 ]
	[ "${stages[0]}" = "inlineFilter|i|FORMAT/DP>10" ]
	[ "${stages[1]}" = "fileFilter|e|FORMAT/GQ>20" ]
}

@test "parse_filter_args handles empty inputs gracefully" {
	local stages=()
	parse_filter_args stages ""
	[ "${#stages[@]}" -eq 0 ]
}

@test "parse_filter_args strips carriage returns from file-based filters" {
	local filter_file="${BATS_TEST_TMPDIR}/crlf_filters.txt"
	# Write a filter with Windows-style CRLF line endings
	printf 'crlfFilter e FORMAT/DP<10\r\n' >"$filter_file"
	local stages=()
	parse_filter_args stages "$filter_file"
	[ "${#stages[@]}" -eq 1 ]
	# Expression must NOT contain a carriage return
	local expr="${stages[0]}"
	[[ "$expr" != *$'\r'* ]]
}

@test "parse_filter_args parses multi-word bcftools expression correctly" {
	local stages=()
	parse_filter_args stages "" 'gatkSNP e TYPE=="SNP" && AS_FS > 60'
	[ "${#stages[@]}" -eq 1 ]
	[ "${stages[0]}" = 'gatkSNP|e|TYPE=="SNP" && AS_FS > 60' ]
}

@test "parse_filter_args parses the real gatk_filters.txt file" {
	local stages=()
	parse_filter_args stages "${REPO_ROOT}/defaults/gatk_filters.txt"
	# gatk_filters.txt has 7 non-empty filter lines
	[ "${#stages[@]}" -eq 7 ]
}
