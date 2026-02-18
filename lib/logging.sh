#!/bin/bash
# lib/logging.sh — Logging functions and command execution wrapper
# Provides: log_msg, debug_msg, error_msg, run_cmd, run_cmd_with_retry
# Setters:  set_log_file, set_debug, set_tmp_dir

[[ -n "${_LIB_LOGGING_LOADED:-}" ]] && return 0
readonly _LIB_LOGGING_LOADED=1

# Global logging state — justified exception (logging config needed by all callers)
_LOG_FILE=""
_DEBUG=false
_TMP_DIR=""

# set_log_file — configure the log file path
set_log_file() {
	_LOG_FILE="$1"
}

# set_debug — configure debug mode ("true"/"false")
set_debug() {
	_DEBUG="$1"
}

# set_tmp_dir — configure the temp directory for run_cmd stderr files
set_tmp_dir() {
	_TMP_DIR="$1"
}

# log_msg — write a timestamped message to log file or stdout
log_msg() {
	local timestamp
	timestamp=$(date +"%Y-%m-%d %H:%M:%S")
	if [[ -n "$_LOG_FILE" ]]; then
		echo "[$timestamp] $1" >>"$_LOG_FILE"
	else
		echo "[$timestamp] $1"
	fi
}

# debug_msg — write a DEBUG-prefixed message (only when debug mode is on)
debug_msg() {
	[[ "$_DEBUG" == "true" ]] || return 0
	log_msg "[DEBUG] $1"
}

# error_msg — write an ERROR-prefixed message to stderr (and log file if set)
error_msg() {
	local timestamp
	timestamp=$(date +"%Y-%m-%d %H:%M:%S")
	# Always write to stderr to honor the function contract
	echo "[$timestamp] ERROR: $1" >&2
	# Additionally append to log file when configured
	if [[ -n "$_LOG_FILE" ]]; then
		echo "[$timestamp] ERROR: $1" >>"$_LOG_FILE"
	fi
}

# run_cmd — execute a command and capture stderr; report details on failure
# Usage: run_cmd <cmd> [args...]
run_cmd() {
	local stderr_file
	# Use _TMP_DIR if set, otherwise fall back to mktemp default (/tmp)
	stderr_file=$(mktemp "${_TMP_DIR:+${_TMP_DIR}/}runcmd-XXXXXX")
	local cmd_name="$1"
	local exit_code=0

	"$@" 2>"$stderr_file" || exit_code=$?

	if [[ "$exit_code" -ne 0 ]]; then
		local stderr_content
		stderr_content=$(cat "$stderr_file")
		rm -f "$stderr_file"
		error_msg "'${cmd_name}' failed (exit ${exit_code})"
		error_msg "Command: $*"
		if [[ -n "$stderr_content" ]]; then
			error_msg "Stderr: ${stderr_content}"
		fi
		return "$exit_code"
	fi

	rm -f "$stderr_file"
	return 0
}

# run_cmd_with_retry — retry wrapper for unreliable network operations (e.g. UCSC MySQL)
# Usage: run_cmd_with_retry <max_attempts> <cmd> [args...]
run_cmd_with_retry() {
	local max_attempts="$1"
	shift
	local attempt=1

	while [[ "$attempt" -le "$max_attempts" ]]; do
		if run_cmd "$@"; then
			return 0
		fi
		if [[ "$attempt" -lt "$max_attempts" ]]; then
			log_msg "Attempt ${attempt}/${max_attempts} failed for '$1', retrying in 2s..."
			sleep 2
		fi
		attempt=$((attempt + 1))
	done

	error_msg "All ${max_attempts} attempts failed for '$1'"
	return 1
}
