#!/bin/bash
# lib/genome.sh — Genome file creation via UCSC MySQL
# Provides: create_genome_file
# Requires: lib/logging.sh (log_msg, debug_msg, error_msg, run_cmd_with_retry)

[[ -n "${_LIB_GENOME_LOADED:-}" ]] && return 0
readonly _LIB_GENOME_LOADED=1

# create_genome_file — fetch chromosome sizes from UCSC MySQL and write a genome file
# Usage: create_genome_file <genome_build> <tmp_dir>
# Prints the path to the created genome file on stdout.
# Returns 1 on failure after max_attempts.
create_genome_file() {
	local genome_build="$1"
	local tmp_dir="$2"
	local genome_file="${tmp_dir}/${genome_build}.genome"
	local max_attempts=3
	local attempt=1
	log_msg "Creating genome file: ${genome_file} for genome build: ${genome_build}"
	while [[ "$attempt" -le "$max_attempts" ]]; do
		if mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A \
			-e "select chrom, size from ${genome_build}.chromInfo" 2>/dev/null \
			| grep -v "^chrom" \
			| sed 's/chr//g' >"$genome_file"; then
			debug_msg "Genome file created: ${genome_file}"
			printf '%s' "$genome_file"
			return 0
		fi
		log_msg "Genome query attempt ${attempt}/${max_attempts} failed. Retrying..."
		((attempt++))
		sleep 2
	done
	error_msg "Failed to create genome file after ${max_attempts} attempts"
	return 1
}
