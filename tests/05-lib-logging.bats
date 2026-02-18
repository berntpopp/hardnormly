#!/usr/bin/env bats
# tests/05-lib-logging.bats — Unit tests for lib/logging.sh (REFR-13)
#
# Verifies logging functions in isolation without any bioinformatics tools.
# Sources only lib/logging.sh in setup().
#
# Functions covered: log_msg, debug_msg, error_msg, run_cmd
# Note: run_cmd_with_retry is implicitly covered via run_cmd; not retested here.

bats_require_minimum_version 1.5.0

setup() {
	load 'test_helper/common'
	_common_setup
	# Reset logging state before each test
	source "${REPO_ROOT}/lib/logging.sh"
	set_log_file ""
	set_debug "false"
	set_tmp_dir "${BATS_TEST_TMPDIR}"
}

# ---------------------------------------------------------------------------
# log_msg — writes to stdout or file
# ---------------------------------------------------------------------------

@test "log_msg writes timestamped message to stdout when no log file is set" {
	run log_msg "hello world"
	assert_success
	assert_output --partial "hello world"
}

@test "log_msg output contains a timestamp bracket" {
	run log_msg "test"
	assert_success
	assert_output --partial "["
}

@test "log_msg writes to file when log file is set" {
	local log="${BATS_TEST_TMPDIR}/test.log"
	set_log_file "$log"
	log_msg "written to file"
	# stdout should be empty
	run log_msg "second message"
	# File must contain our messages
	grep -q "written to file" "$log"
	grep -q "second message" "$log"
}

@test "log_msg does not write to stdout when log file is set" {
	local log="${BATS_TEST_TMPDIR}/test.log"
	set_log_file "$log"
	run log_msg "silent stdout"
	# run captures stdout; should be empty
	assert_output ""
}

# ---------------------------------------------------------------------------
# debug_msg — only outputs when debug is true
# ---------------------------------------------------------------------------

@test "debug_msg outputs nothing when debug is false" {
	set_debug "false"
	run debug_msg "should not appear"
	assert_success
	assert_output ""
}

@test "debug_msg outputs message when debug is true" {
	set_debug "true"
	run debug_msg "debug output"
	assert_success
	assert_output --partial "debug output"
}

@test "debug_msg includes [DEBUG] prefix when debug is true" {
	set_debug "true"
	run debug_msg "my debug message"
	assert_success
	assert_output --partial "[DEBUG]"
}

# ---------------------------------------------------------------------------
# error_msg — writes to stderr
# ---------------------------------------------------------------------------

@test "error_msg writes to stderr (stdout is empty)" {
	run --separate-stderr error_msg "something went wrong"
	# stdout must be empty (message goes to stderr only)
	assert_output ""
}

@test "error_msg message appears on stderr" {
	run --separate-stderr error_msg "critical failure"
	# shellcheck disable=SC2154
	[[ "$stderr" == *"critical failure"* ]]
}

@test "error_msg message includes ERROR: prefix" {
	run --separate-stderr error_msg "bad input"
	# shellcheck disable=SC2154
	[[ "$stderr" == *"ERROR:"* ]]
}

# ---------------------------------------------------------------------------
# run_cmd — command execution with error capture
# ---------------------------------------------------------------------------

@test "run_cmd succeeds for a valid command" {
	run run_cmd true
	assert_success
}

@test "run_cmd returns non-zero for failing command" {
	run run_cmd false
	assert_failure
}

@test "run_cmd cleans up stderr temp file on success" {
	local before after
	before=$(find "${BATS_TEST_TMPDIR}" -maxdepth 1 -name 'runcmd-*' | wc -l)
	run_cmd true
	after=$(find "${BATS_TEST_TMPDIR}" -maxdepth 1 -name 'runcmd-*' | wc -l)
	# No temp file should remain after successful run_cmd
	[[ "$after" -eq "$before" ]]
}

@test "run_cmd cleans up stderr temp file on failure" {
	local before after
	before=$(find "${BATS_TEST_TMPDIR}" -maxdepth 1 -name 'runcmd-*' | wc -l)
	run_cmd false || true
	after=$(find "${BATS_TEST_TMPDIR}" -maxdepth 1 -name 'runcmd-*' | wc -l)
	# No temp file should remain after failed run_cmd
	[[ "$after" -eq "$before" ]]
}

@test "run_cmd reports command name in error output on failure" {
	run --separate-stderr run_cmd false
	assert_failure
	# shellcheck disable=SC2154
	[[ "$stderr" == *"false"* ]]
}

@test "run_cmd captures and reports stderr content on failure" {
	# Use a shell command that writes to stderr and fails
	run --separate-stderr run_cmd bash -c 'echo "custom error output" >&2; exit 1'
	assert_failure
	# shellcheck disable=SC2154
	[[ "$stderr" == *"custom error output"* ]]
}

@test "run_cmd passes arguments to command correctly" {
	local out_file="${BATS_TEST_TMPDIR}/touch_target"
	run run_cmd touch "$out_file"
	assert_success
	[[ -f "$out_file" ]]
}
