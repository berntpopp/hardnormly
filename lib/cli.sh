#!/bin/bash
# lib/cli.sh — Argument parsing, validation, and filter parsing
# Provides: show_help, show_version, parse_args, validate_args, parse_filter_args,
#           show_help_generate_inclusion_bed, show_help_generate_exclusion_bed
# Requires: lib/logging.sh (log_msg used in validate_args)

[[ -n "${_LIB_CLI_LOADED:-}" ]] && return 0
readonly _LIB_CLI_LOADED=1

# show_help — print compact subcommand-aware usage and exit with code 1
show_help() {
	cat <<'HELP'
Usage: hardnormly.sh <subcommand> [options]
       hardnormly.sh [options]           (legacy, implies run-pipeline)

Subcommands:
  run-pipeline             Normalize and filter a VCF file (default)
  generate-inclusion-bed   Merge BED files into a combined inclusion region
  generate-exclusion-bed   Merge BED files into a combined exclusion region

Options (run-pipeline):
  -v, --vcf FILE           Input VCF file (required)
  -f, --fasta FILE         Reference FASTA file (required)
  -o, --output FILE        Output VCF file (default: stdout)
  -b, --include-bed FILE   Include BED file (repeatable)
  -e, --exclude-bed FILE   Exclude BED file (repeatable)
  -g, --genome FILE        Genome file for slop (skips UCSC fetch)
  --genome-build BUILD     Genome build for UCSC fetch (default: hg19)
  --slop N                 Region padding in bp (default: 20)
  --caller CALLER          Auto-select filter file: gatk, gatk-no-as, freebayes
  --filters-file FILE      Filter expression file
  --filters EXPR           Inline filter: "name action expression"
  --strip-annotations LIST Remove INFO fields before filtering (e.g. INFO/CSQ,INFO/ANN)
  --only-pass              Keep only PASS variants in output
  --generate-stats         Generate bcftools stats file
  --plot-stats             Plot stats (requires --generate-stats, --plot-output-dir)
  --plot-output-dir DIR    Directory for plot output
  --auto-index             Auto-index compressed output VCF
  --tmp-dir DIR            Custom temp directory
  --no-cleanup             Preserve temp directory after run
  --log-file FILE          Log to file instead of stdout
  --debug                  Enable verbose debug output
  --version                Show version
  -h, --help               Show this help

Filter file format (--filters-file):
  Each line: <filter_name> <e|i> <bcftools_expression>
  e = exclude (soft-filter), i = include (keep matching)

Example:
  hardnormly.sh run-pipeline -v input.vcf.gz -f ref.fasta -o output.vcf.gz

Run 'hardnormly.sh <subcommand> --help' for subcommand-specific options.
HELP
	exit 1
}

# show_version — print version string and exit 0
# Usage: show_version "$version"
show_version() {
	echo "Version: $1"
	exit 0
}

# show_help_generate_inclusion_bed — help for generate-inclusion-bed subcommand
show_help_generate_inclusion_bed() {
	cat <<'HELP'
Usage: hardnormly.sh generate-inclusion-bed [options]

Merge one or more BED files into a combined inclusion region file.
Applies bedtools intersect (for multiple files) and slop padding.

Options:
  -b, --include-bed FILE   Input BED file (required, repeatable)
  -g, --genome FILE        Genome file for slop operation (required)
  -o, --output FILE        Output merged BED file (required)
  --slop N                 Region padding in bp (default: 20)
  -v, --verbose            Show progress messages
  -h, --help               Show this help
HELP
	exit 0
}

# show_help_generate_exclusion_bed — help for generate-exclusion-bed subcommand
show_help_generate_exclusion_bed() {
	cat <<'HELP'
Usage: hardnormly.sh generate-exclusion-bed [options]

Merge one or more BED files into a combined exclusion region file.
Uses bedtools multiinter for union of all regions.

Options:
  -e, --exclude-bed FILE   Input BED file (required, repeatable)
  -o, --output FILE        Output merged BED file (required)
  -v, --verbose            Show progress messages
  -h, --help               Show this help
HELP
	exit 0
}

# parse_args — parse command-line arguments and set caller-scope variables
# Caller must declare all variables before calling parse_args.
# Variables set: include_bed_files, exclude_bed_files, filters (arrays)
#   filters_file, fasta_file, vcf_file, output_vcf, genome_file, genome_build,
#   log_file, plot_output_dir, tmp_dir, debug, cleanup, only_pass,
#   generate_stats, plot_stats, auto_index, slop, version,
#   caller, strip_annotations
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
			--caller)
				[[ -z "$2" || "$2" == -* ]] && {
					echo "Error: Argument for $1 is missing"
					show_help
				}
				caller="$2"
				shift
				;;
			--strip-annotations)
				[[ -z "$2" || "$2" == -* ]] && {
					echo "Error: Argument for $1 is missing"
					show_help
				}
				strip_annotations="$2"
				shift
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
