#!/bin/bash
# lib/stats.sh — VCF stats generation and plotting via bcftools stats / plot-vcfstats
# Provides: generate_stats, plot_stats_output
# Requires: lib/logging.sh (log_msg, debug_msg, run_cmd)

[[ -n "${_LIB_STATS_LOADED:-}" ]] && return 0
readonly _LIB_STATS_LOADED=1

# generate_stats — run bcftools stats on a VCF and write output to a stats file
# Usage: generate_stats <output_vcf> <stats_output>
# Returns 1 on failure.
generate_stats() {
	local output_vcf="$1"
	local stats_output="$2"

	run_cmd bcftools stats "$output_vcf" >"$stats_output"
}

# plot_stats_output — run plot-vcfstats and log its output
# Usage: plot_stats_output <stats_output> <plot_output_dir> <tmp_dir>
# Returns 1 on failure.
plot_stats_output() {
	local stats_output="$1"
	local plot_output_dir="$2"
	local tmp_dir="$3"

	local plot_output
	plot_output=$(mktemp "${tmp_dir}/plot-vcfstats-XXXXXX")

	plot-vcfstats "$stats_output" -p "$plot_output_dir" >"$plot_output" 2>&1 \
		|| {
			log_msg "Error: Failed to plot stats."
			while IFS= read -r line; do
				log_msg "Plot-vcfstats output: $line"
			done <"$plot_output"
			rm -f "$plot_output"
			return 1
		}

	while IFS= read -r line; do
		log_msg "Plot-vcfstats output: $line"
	done <"$plot_output"

	log_msg "Plots saved to $plot_output_dir"
	rm -f "$plot_output"
}
