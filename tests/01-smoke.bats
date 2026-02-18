#!/usr/bin/env bats
# tests/01-smoke.bats — CLI argument validation (TFWK-01)
# These tests run WITHOUT bioinformatics tools (no _require_tools call)

setup() {
	load 'test_helper/common'
	_common_setup
}

@test "--help shows usage and exits non-zero" {
	run "$HARDNORMLY" --help
	assert_failure
	assert_output --partial "Usage:"
}

@test "--version prints version string and exits 0" {
	run "$HARDNORMLY" --version
	assert_success
	assert_output --partial "Version:"
}

@test "no arguments exits non-zero" {
	run "$HARDNORMLY"
	assert_failure
}

@test "missing -v flag exits non-zero" {
	run "$HARDNORMLY" -f /dev/null
	assert_failure
}

@test "missing -f flag exits non-zero" {
	run "$HARDNORMLY" -v /dev/null
	assert_failure
}

@test "unknown parameter exits non-zero" {
	run "$HARDNORMLY" --bogus-flag
	assert_failure
	assert_output --partial "Unknown parameter"
}

@test "run-pipeline --help shows usage and exits non-zero" {
	run "$HARDNORMLY" run-pipeline --help
	assert_failure
	assert_output --partial "Usage:"
}

@test "run-pipeline --version prints version and exits 0" {
	run "$HARDNORMLY" run-pipeline --version
	assert_success
	assert_output --partial "Version:"
}

@test "unknown subcommand exits non-zero" {
	run "$HARDNORMLY" bogus-subcommand
	assert_failure
	assert_output --partial "Unknown subcommand"
}

@test "--help output lists subcommands" {
	run "$HARDNORMLY" --help
	assert_failure
	assert_output --partial "run-pipeline"
	assert_output --partial "generate-inclusion-bed"
	assert_output --partial "generate-exclusion-bed"
}

@test "--help output shows --caller flag" {
	run "$HARDNORMLY" --help
	assert_failure
	assert_output --partial "--caller"
}

@test "--help output shows --strip-annotations flag" {
	run "$HARDNORMLY" --help
	assert_failure
	assert_output --partial "--strip-annotations"
}
