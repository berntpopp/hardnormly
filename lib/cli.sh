#!/bin/bash
# lib/cli.sh — Argument parsing, validation, and filter parsing
# Provides: show_help, show_version, parse_args, validate_args, parse_filter_args
# Requires: lib/logging.sh (log_msg used in validate_args)

[[ -n "${_LIB_CLI_LOADED:-}" ]] && return 0
readonly _LIB_CLI_LOADED=1

# show_help — print usage message and exit with code 1
show_help() {
	echo "Usage: $0 -v <vcf_file> -f <fasta_file> [-o <output_vcf>] [options]"
	echo ""
	echo "Options:"
	echo "  -v, --vcf            Input VCF file (required). The variant call format file to be processed."
	echo "  -f, --fasta          Reference FASTA file for normalization (required). The reference genome sequence in FASTA format."
	echo "  -b, --include-bed    BED file(s) for inclusion. Specifies regions to include. Can specify multiple BED files."
	echo "  -e, --exclude-bed    BED file(s) for exclusion. Specifies regions to exclude. Can specify multiple BED files."
	echo "  -g, --genome         Genome file for slop operation. A file defining chromosome sizes for applying padding. If provided, it skips genome file generation."
	echo "  --genome-build       Genome build to use for UCSC MySQL query (default: hg19). If no genome file is provided, this will fetch chromosome sizes."
	echo "  --slop               Slop size for region padding (default: 20bp). Adds padding to the BED regions during processing."
	echo "  -o, --output         Output VCF file. If not specified, the result will be sent to stdout."
	echo "  --filters            Inline bcftools filter expression. You can specify multiple filters in the format: filter_name action expression."
	echo "  --filters-file       File containing bcftools filter expressions. Each line should be in the format: filter_name action expression."
	echo "  --only-pass          Filter to retain only variants with a PASS status in the VCF."
	echo "  --generate-stats     Generate a statistics file from the output VCF using bcftools stats."
	echo "  --plot-stats         Plot the stats file using plot-vcfstats. Requires --generate-stats."
	echo "  --plot-output-dir    Directory to save the plots. Required if --plot-stats is set."
	echo "  --tmp-dir            Temporary directory to use. By default, a unique directory is created using mktemp."
	echo "  --no-cleanup         Do not clean up the temporary directory after execution. Useful for debugging."
	echo "  --log-file           File to write logs to. If not provided, logs will be written to stdout."
	echo "  --auto-index         Automatically index the output VCF (if compressed). Adds '-W' to bcftools view."
	echo "  --debug              Enable debug mode. Prints all executed commands and detailed messages for troubleshooting."
	echo "  --version            Display the script version."
	echo "  -h, --help           Display this help message."
	exit 1
}

# show_version — print version string and exit 0
# Usage: show_version "$version"
show_version() {
	echo "Version: $1"
	exit 0
}

# parse_args — parse command-line arguments and set caller-scope variables
# Caller must declare all variables before calling parse_args.
# Variables set: include_bed_files, exclude_bed_files, filters (arrays)
#   filters_file, fasta_file, vcf_file, output_vcf, genome_file, genome_build,
#   log_file, plot_output_dir, tmp_dir, debug, cleanup, only_pass,
#   generate_stats, plot_stats, auto_index, slop, version
# shellcheck disable=SC2034,SC2154
parse_args() {
	while [[ "$#" -gt 0 ]]; do
		case "$1" in
			-v | --vcf)
				[[ -z "$2" || "$2" == -* ]] && {
					echo "Error: Argument for $1 is missing"
					show_help
				}
				vcf_file="$2"
				shift
				;;
			-f | --fasta)
				[[ -z "$2" || "$2" == -* ]] && {
					echo "Error: Argument for $1 is missing"
					show_help
				}
				fasta_file="$2"
				shift
				;;
			-b | --include-bed)
				[[ -z "$2" || "$2" == -* ]] && {
					echo "Error: Argument for $1 is missing"
					show_help
				}
				include_bed_files+=("$2")
				shift
				;;
			-e | --exclude-bed)
				[[ -z "$2" || "$2" == -* ]] && {
					echo "Error: Argument for $1 is missing"
					show_help
				}
				exclude_bed_files+=("$2")
				shift
				;;
			-g | --genome)
				[[ -z "$2" || "$2" == -* ]] && {
					echo "Error: Argument for $1 is missing"
					show_help
				}
				genome_file="$2"
				shift
				;;
			--genome-build)
				[[ -z "$2" || "$2" == -* ]] && {
					echo "Error: Argument for $1 is missing"
					show_help
				}
				genome_build="$2"
				shift
				;;
			--slop)
				[[ -z "$2" || "$2" == -* ]] && {
					echo "Error: Argument for $1 is missing"
					show_help
				}
				slop="$2"
				shift
				;;
			-o | --output)
				[[ -z "$2" || "$2" == -* ]] && {
					echo "Error: Argument for $1 is missing"
					show_help
				}
				output_vcf="$2"
				shift
				;;
			--filters)
				[[ -z "$2" || "$2" == -* ]] && {
					echo "Error: Argument for $1 is missing"
					show_help
				}
				filters+=("$2")
				shift
				;;
			--filters-file)
				[[ -z "$2" || "$2" == -* ]] && {
					echo "Error: Argument for $1 is missing"
					show_help
				}
				filters_file="$2"
				shift
				;;
			--plot-output-dir)
				[[ -z "$2" || "$2" == -* ]] && {
					echo "Error: Argument for $1 is missing"
					show_help
				}
				plot_output_dir="$2"
				shift
				;;
			--only-pass)
				only_pass=true
				;;
			--generate-stats)
				generate_stats=true
				;;
			--plot-stats)
				plot_stats=true
				;;
			--tmp-dir)
				[[ -z "$2" || "$2" == -* ]] && {
					echo "Error: Argument for $1 is missing"
					show_help
				}
				tmp_dir="$2"
				shift
				;;
			--no-cleanup)
				cleanup=false
				;;
			--log-file)
				[[ -z "$2" || "$2" == -* ]] && {
					echo "Error: Argument for $1 is missing"
					show_help
				}
				log_file="$2"
				shift
				;;
			--auto-index)
				auto_index=true
				;;
			--debug)
				debug=true
				;;
			--version)
				show_version "$version"
				;;
			-h | --help)
				show_help
				;;
			*)
				echo "Unknown parameter: $1"
				show_help
				;;
		esac
		shift
	done
}

# validate_args — validate parsed arguments, exit with help on failure
# Requires log_msg (from lib/logging.sh, sourced before lib/cli.sh)
validate_args() {
	# --plot-stats requires both --generate-stats and --plot-output-dir
	if [[ "$plot_stats" == "true" ]] && { [[ -z "$plot_output_dir" ]] || [[ "$generate_stats" != "true" ]]; }; then
		log_msg "Error: --plot-stats requires both --generate-stats and --plot-output-dir."
		exit 1
	fi

	# Required parameters
	if [[ -z "$vcf_file" || -z "$fasta_file" ]]; then
		log_msg "Error: Missing required parameters."
		show_help
	fi
}

# parse_filter_args — parse inline and file-based filters into a named array
# Usage: parse_filter_args <array_name> <filters_file> [inline_filter...]
# Each inline_filter is "name action expression" (space-separated triple).
# Each entry appended to <array_name> as "name|action|expression".
parse_filter_args() {
	# shellcheck disable=SC2178
	local -n _stages_ref="$1"
	local _filters_file="$2"
	shift 2
	local _inline_filters=("$@")

	# Inline filters
	local _filter _filter_name _filter_action _filter_expr
	for _filter in "${_inline_filters[@]}"; do
		_filter_name=""
		_filter_action=""
		_filter_expr=""
		IFS=" " read -r _filter_name _filter_action _filter_expr <<<"$_filter"
		_stages_ref+=("${_filter_name}|${_filter_action}|${_filter_expr}")
	done

	# File-based filters
	if [[ -n "$_filters_file" ]]; then
		while IFS=" " read -r _filter_name _filter_action _filter_expr; do
			_filter_expr=$(tr -d '\r\n' <<<"$_filter_expr")
			_stages_ref+=("${_filter_name}|${_filter_action}|${_filter_expr}")
		done <"$_filters_file"
	fi
}
