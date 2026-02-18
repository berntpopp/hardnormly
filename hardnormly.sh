#!/bin/bash

# Script version
version="0.7.0"

set -Eeuo pipefail

# Resolve script directory for sourcing lib/ modules
_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Source library modules (logging first — cli.sh depends on log_msg)
source "${_SCRIPT_DIR}/lib/logging.sh"
source "${_SCRIPT_DIR}/lib/cli.sh"
source "${_SCRIPT_DIR}/lib/genome.sh"
source "${_SCRIPT_DIR}/lib/bed.sh"
source "${_SCRIPT_DIR}/lib/annotate.sh"
source "${_SCRIPT_DIR}/lib/normalize.sh"
source "${_SCRIPT_DIR}/lib/filter.sh"
source "${_SCRIPT_DIR}/lib/stats.sh"

# Default values for parameters
include_bed_files=()
exclude_bed_files=()
filters=()
filters_file=""
fasta_file=""
vcf_file=""
output_vcf=""
genome_file=""
genome_build="hg19"                           # Default genome build
debug=false                                   # Debug mode disabled by default
log_file=""                                   # Default to no log file
tmp_dir=$(mktemp -d -t hardnormly-XXXXXXXXXX) # Use mktemp for a unique tmp directory
cleanup=true                                  # Default to cleaning up the temporary directory
slop=20                                       # Default slop value (in base pairs)
only_pass=false                               # Option to filter only PASS variants
generate_stats=false                          # Option to generate stats file
plot_stats=false                              # Option to plot the stats
plot_output_dir=""
auto_index=false     # New option: auto-index output if compressed
caller=""            # --caller flag (gatk or freebayes)
strip_annotations="" # --strip-annotations comma-separated INFO fields

# Cleanup handler — always runs on EXIT, preserving the original exit code
cleanup_handler() {
	local exit_code=$?
	trap - EXIT # Prevent recursive re-entry if exit is called within this handler
	if [[ "${cleanup:-true}" == "true" ]] && [[ -d "${tmp_dir:-}" ]]; then
		if [[ "${debug:-false}" == "true" ]]; then
			find "$tmp_dir" -type f -print -delete 2>/dev/null || true
		else
			rm -rf "$tmp_dir" 2>/dev/null || true
		fi
	fi
	exit "$exit_code"
}

trap cleanup_handler EXIT

# cmd_generate_inclusion_bed — merge BED files into a combined inclusion region file
cmd_generate_inclusion_bed() {
	local bed_files=()
	local genome_file=""
	local output_file=""
	local slop=20
	local verbose=false

	while [[ "$#" -gt 0 ]]; do
		case "$1" in
			-b | --include-bed)
				bed_files+=("$2")
				shift
				;;
			-g | --genome)
				genome_file="$2"
				shift
				;;
			-o | --output)
				output_file="$2"
				shift
				;;
			--slop)
				slop="$2"
				shift
				;;
			-v | --verbose)
				verbose=true
				;;
			-h | --help)
				show_help_generate_inclusion_bed
				;;
			*)
				echo "Error: Unknown option '$1'" >&2
				show_help_generate_inclusion_bed
				;;
		esac
		shift
	done

	# Validate required args
	if [[ "${#bed_files[@]}" -eq 0 ]]; then
		echo "Error: At least one -b/--include-bed file is required" >&2
		exit 1
	fi
	if [[ -z "$genome_file" ]]; then
		echo "Error: -g/--genome file is required for slop operation" >&2
		exit 1
	fi
	if [[ -z "$output_file" ]]; then
		echo "Error: -o/--output is required" >&2
		exit 1
	fi

	# Configure logging for verbose mode
	if [[ "$verbose" == true ]]; then
		set_log_file ""
		set_debug false
	fi

	# Create temp dir for intermediate files
	local tmp_dir
	tmp_dir=$(mktemp -d -t hardnormly-gen-inc-XXXXXXXXXX)
	# shellcheck disable=SC2064
	trap "rm -rf '$tmp_dir'" EXIT

	# Normalize each BED file
	local normalized_files=()
	local bed_file normalized_file
	for bed_file in "${bed_files[@]}"; do
		normalized_file="$tmp_dir/$(basename "$bed_file").normalized.bed"
		normalize_bed "$bed_file" "1" "$normalized_file"
		normalized_files+=("$normalized_file")
		[[ "$verbose" == true ]] && log_msg "Normalized: $bed_file"
	done

	# Merge with slop
	merge_include_beds "$tmp_dir/merged.bed" "$slop" "$genome_file" "${normalized_files[@]}"
	[[ "$verbose" == true ]] && log_msg "Merged ${#bed_files[@]} BED file(s) with slop=${slop}bp"

	# Copy to output
	cp "$tmp_dir/merged.bed" "$output_file"
	[[ "$verbose" == true ]] && log_msg "Output written to: $output_file"
	return 0
}

# cmd_generate_exclusion_bed — merge BED files into a combined exclusion region file
cmd_generate_exclusion_bed() {
	local bed_files=()
	local output_file=""
	local verbose=false

	while [[ "$#" -gt 0 ]]; do
		case "$1" in
			-e | --exclude-bed)
				bed_files+=("$2")
				shift
				;;
			-o | --output)
				output_file="$2"
				shift
				;;
			-v | --verbose)
				verbose=true
				;;
			-h | --help)
				show_help_generate_exclusion_bed
				;;
			*)
				echo "Error: Unknown option '$1'" >&2
				show_help_generate_exclusion_bed
				;;
		esac
		shift
	done

	# Validate required args
	if [[ "${#bed_files[@]}" -eq 0 ]]; then
		echo "Error: At least one -e/--exclude-bed file is required" >&2
		exit 1
	fi
	if [[ -z "$output_file" ]]; then
		echo "Error: -o/--output is required" >&2
		exit 1
	fi

	# Configure logging for verbose mode
	if [[ "$verbose" == true ]]; then
		set_log_file ""
		set_debug false
	fi

	# Create temp dir for intermediate files
	local tmp_dir
	tmp_dir=$(mktemp -d -t hardnormly-gen-exc-XXXXXXXXXX)
	# shellcheck disable=SC2064
	trap "rm -rf '$tmp_dir'" EXIT

	# Normalize each BED file
	local normalized_files=()
	local bed_file normalized_file
	for bed_file in "${bed_files[@]}"; do
		normalized_file="$tmp_dir/$(basename "$bed_file").normalized.bed"
		normalize_bed "$bed_file" "1" "$normalized_file"
		normalized_files+=("$normalized_file")
		[[ "$verbose" == true ]] && log_msg "Normalized: $bed_file"
	done

	# Merge exclusion regions
	merge_exclude_beds "$tmp_dir/merged.bed" "${normalized_files[@]}"
	[[ "$verbose" == true ]] && log_msg "Merged ${#bed_files[@]} BED file(s)"

	# Copy to output
	cp "$tmp_dir/merged.bed" "$output_file"
	[[ "$verbose" == true ]] && log_msg "Output written to: $output_file"
	return 0
}

# Subcommand dispatcher — routes to subcommand handler or falls through to run-pipeline
_subcommand="${1:-}"
case "$_subcommand" in
	run-pipeline)
		shift
		parse_args "$@"
		;;
	generate-inclusion-bed)
		shift
		cmd_generate_inclusion_bed "$@"
		exit 0
		;;
	generate-exclusion-bed)
		shift
		cmd_generate_exclusion_bed "$@"
		exit 0
		;;
	"" | -h | --help)
		show_help
		;;
	--version)
		show_version "$version"
		;;
	-*)
		# Legacy mode: first arg is a flag, route to run-pipeline implicitly
		parse_args "$@"
		;;
	*)
		echo "Error: Unknown subcommand '${_subcommand}'" >&2
		echo "Run 'hardnormly.sh --help' for usage information." >&2
		exit 1
		;;
esac

# Configure logging module with parsed values
set_log_file "$log_file"
set_debug "$debug"

# Resolve --caller to a filter file path (if set)
if [[ -n "$caller" ]]; then
	case "$caller" in
		gatk)
			_caller_file="${_SCRIPT_DIR}/defaults/gatk_filters.txt"
			;;
		freebayes)
			_caller_file="${_SCRIPT_DIR}/defaults/freebayes_filters.txt"
			;;
		*)
			echo "Error: Unknown --caller '$caller'. Valid values: gatk, freebayes" >&2
			exit 1
			;;
	esac
	if [[ -n "$filters_file" ]]; then
		log_msg "Warning: Both --caller and --filters-file provided; --filters-file takes precedence."
	else
		filters_file="$_caller_file"
	fi
fi

# Validate parsed arguments (delegates to lib/cli.sh)
validate_args

# Enable debugging if --debug flag is set
if [[ "$debug" == "true" ]]; then
	set -x # Enable command tracing (prints every command)
fi

# Create the temporary directory if it doesn't exist
mkdir -p "$tmp_dir"
set_tmp_dir "$tmp_dir"
debug_msg "Temporary directory set to: $tmp_dir"

# Step 1: Create the genome file (if not provided)
if [[ -z "$genome_file" ]]; then
	log_msg "No genome file provided. Creating genome file for genome build: $genome_build"
	genome_file=$(
		set -e
		create_genome_file "$genome_build" "$tmp_dir"
	)
else
	log_msg "Using provided genome file: $genome_file"
fi

# Step 2: Normalize and process BED files for inclusion and exclusion
normalized_include_bed_files=()
normalized_exclude_bed_files=()

# Normalize include BED files
for bed_file in "${include_bed_files[@]}"; do
	normalized_file="$tmp_dir/$(basename "$bed_file").normalized.bed"
	normalize_bed "$bed_file" "1" "$normalized_file"
	normalized_include_bed_files+=("$normalized_file")
done

# Normalize exclude BED files
for bed_file in "${exclude_bed_files[@]}"; do
	normalized_file="$tmp_dir/$(basename "$bed_file").normalized.bed"
	normalize_bed "$bed_file" "1" "$normalized_file"
	normalized_exclude_bed_files+=("$normalized_file")
done

# Merge and compress/index include BED files
if [[ "${#normalized_include_bed_files[@]}" -gt 0 ]]; then
	merge_include_beds \
		"$tmp_dir/merged_include_regions.bed" \
		"$slop" \
		"$genome_file" \
		"${normalized_include_bed_files[@]}"
	compress_index_bed "$tmp_dir/merged_include_regions.bed"
	log_msg "Inclusion BED files normalized, merged, and indexed: $tmp_dir/merged_include_regions.bed.gz"
else
	log_msg "No inclusion BED files provided; skipping inclusion annotation."
fi

# Merge and compress/index exclude BED files
if [[ "${#normalized_exclude_bed_files[@]}" -gt 0 ]]; then
	merge_exclude_beds \
		"$tmp_dir/merged_exclude_regions.bed" \
		"${normalized_exclude_bed_files[@]}"
	compress_index_bed "$tmp_dir/merged_exclude_regions.bed"
	log_msg "Exclusion BED files normalized, merged, and indexed: $tmp_dir/merged_exclude_regions.bed.gz"
else
	log_msg "No exclusion BED files provided; skipping exclusion annotation."
fi

# Step 3: Create header files for INFO fields with the correct format
if [[ -f "$tmp_dir/merged_include_regions.bed.gz" ]]; then
	create_header_file "INCLUDE_REGION" "Included region" "$tmp_dir/include_regions.hdr"
fi
if [[ -f "$tmp_dir/merged_exclude_regions.bed.gz" ]]; then
	create_header_file "EXCLUDE_REGION" "Excluded region" "$tmp_dir/exclude_regions.hdr"
fi

# Step 4: Annotate the VCF file with the BED regions
log_msg "Annotating VCF with BED regions..."

# Annotate with inclusion regions if the file exists
if [[ -f "$tmp_dir/merged_include_regions.bed.gz" ]]; then
	annotate_vcf_with_regions \
		"$vcf_file" \
		"$tmp_dir/merged_include_regions.bed.gz" \
		"$tmp_dir/include_regions.hdr" \
		"INCLUDE_REGION" \
		"$tmp_dir/temp_include_annotated.vcf.gz"
	vcf_file="$tmp_dir/temp_include_annotated.vcf.gz"
	debug_msg "Annotated VCF with inclusion regions: $vcf_file"
else
	debug_msg "Skipping annotation with inclusion regions because the file does not exist."
fi

# Annotate with exclusion regions if the file exists
if [[ -f "$tmp_dir/merged_exclude_regions.bed.gz" ]]; then
	annotate_vcf_with_regions \
		"$vcf_file" \
		"$tmp_dir/merged_exclude_regions.bed.gz" \
		"$tmp_dir/exclude_regions.hdr" \
		"EXCLUDE_REGION" \
		"$tmp_dir/temp_exclude_annotated.vcf.gz"
	vcf_file="$tmp_dir/temp_exclude_annotated.vcf.gz"
	debug_msg "Annotated VCF with exclusion regions: $vcf_file"
else
	debug_msg "Skipping annotation with exclusion regions because the file does not exist."
fi

# Step 4.5: Strip specified INFO annotations (if --strip-annotations was provided)
if [[ -n "$strip_annotations" ]]; then
	log_msg "Stripping annotations: $strip_annotations"
	strip_vcf_annotations "$vcf_file" "$strip_annotations" "$tmp_dir/temp_stripped.vcf.gz"
	vcf_file="$tmp_dir/temp_stripped.vcf.gz"
	debug_msg "Stripped annotations from VCF: $vcf_file"
fi

# Step 5: Normalize the VCF file and write to an intermediate file
normalized_vcf="$tmp_dir/normalized.vcf.gz"
normalize_vcf "$vcf_file" "$fasta_file" "$normalized_vcf" "$tmp_dir"

# Step 6: Apply filters to the normalized VCF
log_msg "Filtering the normalized VCF file..."

# Initialize the filter pipeline with fill-tags applied to a BCF temp file
init_filter_pipeline "$normalized_vcf" "$tmp_dir"

# Build filter stages array — each entry encodes: "name|action|expression"
filter_stages=()

# Region-based filters (depend on BED processing results above)
if [[ -f "$tmp_dir/merged_include_regions.bed.gz" ]]; then
	filter_stages+=("NOT_IN_INCLUDE_REGION|e|INFO/INCLUDE_REGION!=1")
fi
if [[ -f "$tmp_dir/merged_exclude_regions.bed.gz" ]]; then
	filter_stages+=("IN_EXCLUDE_REGION|e|INFO/EXCLUDE_REGION=1")
fi

# Parse inline and file-based filters into filter_stages (delegates to lib/cli.sh)
parse_filter_args filter_stages "$filters_file" "${filters[@]}"

# Apply each filter stage sequentially via BCF temp files
apply_filter_stages "$tmp_dir" "${filter_stages[@]}"

# Write final output (format detection, optional PASS filter, optional auto-index)
write_filtered_output "$tmp_dir" "$output_vcf" "$only_pass" "$auto_index"

# Step 7: Generate stats file if the --generate-stats option is set and output_vcf is provided
if [[ "$generate_stats" == "true" ]] && [[ -n "$output_vcf" ]]; then
	stats_output="${output_vcf%.vcf.gz}.stats.txt"
	debug_msg "Generating stats file: $stats_output"
	generate_stats "$output_vcf" "$stats_output"
	log_msg "Stats file saved to $stats_output"

	# If plotting is requested (non-fatal: plot failure does not abort the pipeline)
	if [[ "$plot_stats" == "true" ]]; then
		log_msg "Plotting stats to $plot_output_dir"
		_plot_rc=0
		# shellcheck disable=SC2310
		if ! plot_stats_output "$stats_output" "$plot_output_dir" "$tmp_dir"; then
			_plot_rc=1
		fi
		if [[ "$_plot_rc" -ne 0 ]]; then
			log_msg "Warning: plot-vcfstats failed; pipeline continues."
		fi
	fi
else
	debug_msg "Stats generation skipped (either --generate-stats was not set or no output file provided)."
fi

if [[ -n "$output_vcf" ]]; then
	log_msg "Filtered VCF saved to $output_vcf"
fi
