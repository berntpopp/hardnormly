#!/bin/bash

# Script version
version="0.6.0"

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
auto_index=false # New option: auto-index output if compressed

# Error handler — fires on any command failure due to set -E (errtrace)
err_handler() {
	local exit_code=$?
	local line_number=$1
	local failed_command="${BASH_COMMAND}"
	echo "ERROR: '${failed_command}' failed (exit ${exit_code}) at line ${line_number}" >&2
}

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

trap 'err_handler ${LINENO}' ERR
trap cleanup_handler EXIT

# Parse command-line arguments (delegates to lib/cli.sh)
parse_args "$@"

# Configure logging module with parsed values
set_log_file "$log_file"
set_debug "$debug"

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

	# If plotting is requested
	if [[ "$plot_stats" == "true" ]]; then
		log_msg "Plotting stats to $plot_output_dir"
		plot_stats_output "$stats_output" "$plot_output_dir" "$tmp_dir"
	fi
else
	debug_msg "Stats generation skipped (either --generate-stats was not set or no output file provided)."
fi

if [[ -n "$output_vcf" ]]; then
	log_msg "Filtered VCF saved to $output_vcf"
fi
