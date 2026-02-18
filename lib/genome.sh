#!/bin/bash
# lib/genome.sh — Genome file creation via UCSC MySQL
# Provides: create_genome_file
# Requires: lib/logging.sh (log_msg, debug_msg, error_msg)

[[ -n "${_LIB_GENOME_LOADED:-}" ]] && return 0
readonly _LIB_GENOME_LOADED=1

# create_genome_file — fetch chromosome sizes from UCSC MySQL and write a genome file
# Usage: create_genome_file <genome_build> <tmp_dir>
# Prints the path to the created genome file on stdout.
# Returns 1 on failure after max_attempts.
# Uses a manual retry loop because this function captures stdout for its return
# value — run_cmd's stderr redirect would conflict with that.
create_genome_file() {
	local genome_build="$1"
	local tmp_dir="$2"
	local genome_file="${tmp_dir}/${genome_build}.genome"
	local max_attempts=3
	local attempt=1
	local mysql_stderr
	log_msg "Creating genome file: ${genome_file} for genome build: ${genome_build}"
	while [[ "$attempt" -le "$max_attempts" ]]; do
		mysql_stderr=$(mktemp "${tmp_dir}/mysql-stderr-XXXXXX")
		if mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A \
			-e "select chrom, size from ${genome_build}.chromInfo" 2>"$mysql_stderr" \
			| grep -v "^chrom" \
			| sed 's/chr//g' >"$genome_file"; then
			rm -f "$mysql_stderr"
			debug_msg "Genome file created: ${genome_file}"
			printf '%s' "$genome_file"
			return 0
		fi
		local stderr_content
		stderr_content=$(cat "$mysql_stderr")
		rm -f "$mysql_stderr"
		if [[ -n "$stderr_content" ]]; then
			error_msg "mysql stderr: ${stderr_content}"
		fi
		if [[ "$attempt" -lt "$max_attempts" ]]; then
			log_msg "Genome query attempt ${attempt}/${max_attempts} failed. Retrying..."
			sleep 2
		fi
		attempt=$((attempt + 1))
	done
	error_msg "Failed to create genome file after ${max_attempts} attempts"
	return 1
}
