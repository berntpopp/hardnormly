#!/bin/bash

# Script version
version="0.5.0"

# Default values for parameters
include_bed_files=()
exclude_bed_files=()
filters=()
filters_file=""
fasta_file=""
vcf_file=""
output_vcf=""
genome_file=""
genome_build="hg19"  # Default genome build
debug=false  # Debug mode disabled by default
log_file=""  # Default to no log file
tmp_dir=$(mktemp -d -t hardnormly-XXXXXXXXXX)  # Use mktemp for a unique tmp directory
cleanup=true  # Default to cleaning up the temporary directory
slop=20  # Default slop value (in base pairs)
only_pass=false  # Option to filter only PASS variants
generate_stats=false  # Option to generate stats file
plot_stats=false  # Option to plot the stats
plot_output_dir=""

# Ensure the temporary directory is cleaned up on exit or error
trap '[[ $cleanup == true ]] && cleanup_tmp_dir' EXIT

# Function to clean up the temporary directory
cleanup_tmp_dir() {
    if $cleanup; then
        debug_msg "Starting cleanup of temporary directory: $tmp_dir"
        if $debug; then
            find "$tmp_dir" -type f -print -delete
        else
            rm -rf "$tmp_dir"
        fi
        debug_msg "Cleanup of temporary directory completed."
    fi
}

# Function to display help message
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
    echo "  --debug              Enable debug mode. Prints all executed commands and detailed messages for troubleshooting."
    echo "  --version            Display the script version."
    echo "  -h, --help           Display this help message."
    exit 1
}

# Function to log messages
log_msg() {
    local timestamp=$(date +"%Y-%m-%d %H:%M:%S")
    if [[ -n "$log_file" ]]; then
        echo "[$timestamp] $1" >> "$log_file"
    else
        echo "[$timestamp] $1"
    fi
}

# Function to print debug messages
debug_msg() {
    if $debug; then
        log_msg "[DEBUG] $1"
    fi
}

# Parsing command-line arguments
while [[ "$#" -gt 0 ]]; do
    case "$1" in
        -v|--vcf)
            [[ -z "$2" || "$2" == -* ]] && { echo "Error: Argument for $1 is missing"; show_help; }
            vcf_file="$2"; shift ;;
        -f|--fasta)
            [[ -z "$2" || "$2" == -* ]] && { echo "Error: Argument for $1 is missing"; show_help; }
            fasta_file="$2"; shift ;;
        -b|--include-bed)
            [[ -z "$2" || "$2" == -* ]] && { echo "Error: Argument for $1 is missing"; show_help; }
            include_bed_files+=("$2"); shift ;;
        -e|--exclude-bed)
            [[ -z "$2" || "$2" == -* ]] && { echo "Error: Argument for $1 is missing"; show_help; }
            exclude_bed_files+=("$2"); shift ;;
        -g|--genome)
            [[ -z "$2" || "$2" == -* ]] && { echo "Error: Argument for $1 is missing"; show_help; }
            genome_file="$2"; shift ;;
        --genome-build)
            [[ -z "$2" || "$2" == -* ]] && { echo "Error: Argument for $1 is missing"; show_help; }
            genome_build="$2"; shift ;;
        --slop)
            [[ -z "$2" || "$2" == -* ]] && { echo "Error: Argument for $1 is missing"; show_help; }
            slop="$2"; shift ;;
        -o|--output)
            [[ -z "$2" || "$2" == -* ]] && { echo "Error: Argument for $1 is missing"; show_help; }
            output_vcf="$2"; shift ;;
        --filters)
            [[ -z "$2" || "$2" == -* ]] && { echo "Error: Argument for $1 is missing"; show_help; }
            filters+=("$2"); shift ;;
        --filters-file)
            [[ -z "$2" || "$2" == -* ]] && { echo "Error: Argument for $1 is missing"; show_help; }
            filters_file="$2"; shift ;;
        --plot-output-dir)
            [[ -z "$2" || "$2" == -* ]] && { echo "Error: Argument for $1 is missing"; show_help; }
            plot_output_dir="$2"; shift ;;
        --only-pass) only_pass=true ;;
        --generate-stats) generate_stats=true ;;
        --plot-stats) plot_stats=true ;;
        --tmp-dir)
            [[ -z "$2" || "$2" == -* ]] && { echo "Error: Argument for $1 is missing"; show_help; }
            tmp_dir="$2"; shift ;;
        --no-cleanup) cleanup=false ;;
        --log-file)
            [[ -z "$2" || "$2" == -* ]] && { echo "Error: Argument for $1 is missing"; show_help; }
            log_file="$2"; shift ;;
        --debug) debug=true ;;
        --version) echo "Version: $version"; exit 0 ;;
        -h|--help) show_help ;;
        *) echo "Unknown parameter: $1"; show_help ;;
    esac
    shift
done

# Validate --plot-stats requires both --generate-stats and --plot-output-dir
if $plot_stats && { [[ -z "$plot_output_dir" ]] || ! $generate_stats; }; then
    log_msg "Error: --plot-stats requires both --generate-stats and --plot-output-dir."
    exit 1
fi

# Enable debugging if --debug flag is set
if $debug; then
    set -x  # Enable command tracing (prints every command)
fi

# Check required parameters
if [[ -z "$vcf_file" || -z "$fasta_file" ]]; then
    log_msg "Error: Missing required parameters."
    show_help
fi

# Create the temporary directory if it doesn't exist
mkdir -p "$tmp_dir"
debug_msg "Temporary directory set to: $tmp_dir"

# Function to normalize BED files
normalize_bed() {
    local bed_file="$1"
    local annotation="$2"
    local output_file="$3"
    debug_msg "Normalizing BED file: $bed_file with annotation: $annotation"
    awk -v annot="$annotation" '{OFS="\t"; print $1, $2, $3, annot}' "$bed_file" | bedtools sort -i - > "$output_file"
    debug_msg "Normalized BED file written to: $output_file"
}

# Step 1: Create the genome file (if not provided)
if [[ -z "$genome_file" ]]; then
    genome_file="$tmp_dir/${genome_build}.genome"
    log_msg "No genome file provided. Creating genome file: $genome_file for genome build: $genome_build"
    mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -e "select chrom, size from ${genome_build}.chromInfo" | grep -v "^chrom" | sed 's/chr//g' > "$genome_file" || { log_msg "Error: Failed to create genome file"; exit 1; }
    debug_msg "Genome file created: $genome_file"
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
    normalize_bed "$bed_file" "exclude" "$normalized_file"
    normalized_exclude_bed_files+=("$normalized_file")
done

# Combine normalized exclude BED files using bedtools multiinter
if [[ "${#normalized_exclude_bed_files[@]}" -gt 1 ]]; then
    log_msg "Combining normalized exclusion BED files..."
    bedtools multiinter -i "${normalized_exclude_bed_files[@]}" | \
    bedtools sort -i - | \
    awk '{OFS="\t"; print $1, $2, $3, "1"}' > "$tmp_dir/merged_exclude_regions.bed"
    debug_msg "Combined exclusion regions written to: $tmp_dir/merged_exclude_regions.bed"
elif [[ "${#normalized_exclude_bed_files[@]}" -eq 1 ]]; then
    cp "${normalized_exclude_bed_files[0]}" "$tmp_dir/merged_exclude_regions.bed"
    debug_msg "Single exclusion BED file copied to: $tmp_dir/merged_exclude_regions.bed"
fi

# Combine normalized include BED files using bedtools intersect and apply slop
if [[ "${#normalized_include_bed_files[@]}" -gt 1 ]]; then
    log_msg "Intersecting and padding normalized inclusion BED files..."
    bedtools intersect -a "${normalized_include_bed_files[0]}" -b "${normalized_include_bed_files[@]:1}" | \
    bedtools sort -i - | \
    bedtools slop -b "$slop" -g "$genome_file" > "$tmp_dir/merged_include_regions.bed"
    debug_msg "Combined inclusion regions written to: $tmp_dir/merged_include_regions.bed"
elif [[ "${#normalized_include_bed_files[@]}" -eq 1 ]]; then
    log_msg "Padding single normalized inclusion BED file..."
    bedtools slop -b "$slop" -g "$genome_file" -i "${normalized_include_bed_files[0]}" > "$tmp_dir/merged_include_regions.bed"
    debug_msg "Single inclusion BED file padded and written to: $tmp_dir/merged_include_regions.bed"
fi

# Compress and index the merged BED files for inclusion
if [[ -f "$tmp_dir/merged_include_regions.bed" ]]; then
    bgzip -f "$tmp_dir/merged_include_regions.bed"  # Force overwrite
    tabix -p bed "$tmp_dir/merged_include_regions.bed.gz"
    log_msg "Inclusion BED files normalized, merged, and indexed: $tmp_dir/merged_include_regions.bed.gz"
    debug_msg "Inclusion BED files compressed and indexed."
else
    log_msg "No inclusion BED files provided; skipping inclusion annotation."
fi

# Compress and index the merged BED files for exclusion
if [[ -f "$tmp_dir/merged_exclude_regions.bed" ]]; then
    bgzip -f "$tmp_dir/merged_exclude_regions.bed"  # Force overwrite
    tabix -p bed "$tmp_dir/merged_exclude_regions.bed.gz"
    log_msg "Exclusion BED files normalized, merged, and indexed: $tmp_dir/merged_exclude_regions.bed.gz"
    debug_msg "Exclusion BED files compressed and indexed."
else
    log_msg "No exclusion BED files provided; skipping exclusion annotation."
fi

# Step 3: Create header files for INFO fields with the correct format
if [[ -f "$tmp_dir/merged_include_regions.bed.gz" ]]; then
    echo '##INFO=<ID=INCLUDE_REGION,Number=1,Type=Integer,Description="Included region">' > "$tmp_dir/include_regions.hdr"
    debug_msg "Created include_regions.hdr file: $tmp_dir/include_regions.hdr"
fi
if [[ -f "$tmp_dir/merged_exclude_regions.bed.gz" ]]; then
    echo '##INFO=<ID=EXCLUDE_REGION,Number=1,Type=Integer,Description="Excluded region">' > "$tmp_dir/exclude_regions.hdr"
    debug_msg "Created exclude_regions.hdr file: $tmp_dir/exclude_regions.hdr"
fi

# Step 4: Annotate the VCF file with the BED regions
log_msg "Annotating VCF with BED regions..."

# Annotate with inclusion regions if the file exists
if [[ -f "$tmp_dir/merged_include_regions.bed.gz" ]]; then
    bcftools annotate -a "$tmp_dir/merged_include_regions.bed.gz" -h "$tmp_dir/include_regions.hdr" -c CHROM,FROM,TO,INCLUDE_REGION "$vcf_file" -Oz -o "$tmp_dir/temp_include_annotated.vcf.gz"
    if [[ $? -ne 0 ]]; then
        log_msg "Error: Failed to annotate VCF with inclusion regions."
        exit 1
    fi
    vcf_file="$tmp_dir/temp_include_annotated.vcf.gz"
    debug_msg "Annotated VCF with inclusion regions: $vcf_file"
else
    debug_msg "Skipping annotation with inclusion regions because the file does not exist."
fi

# Annotate with exclusion regions if the file exists
if [[ -f "$tmp_dir/merged_exclude_regions.bed.gz" ]]; then
    bcftools annotate -a "$tmp_dir/merged_exclude_regions.bed.gz" -h "$tmp_dir/exclude_regions.hdr" -c CHROM,FROM,TO,EXCLUDE_REGION "$vcf_file" -Oz -o "$tmp_dir/temp_exclude_annotated.vcf.gz"
    if [[ $? -ne 0 ]]; then
        log_msg "Error: Failed to annotate VCF with exclusion regions."
        exit 1
    fi
    vcf_file="$tmp_dir/temp_exclude_annotated.vcf.gz"
    debug_msg "Annotated VCF with exclusion regions: $vcf_file"
else
    debug_msg "Skipping annotation with exclusion regions because the file does not exist."
fi

# Step 5: Normalize the VCF file and write to an intermediate file
normalized_vcf="$tmp_dir/normalized.vcf.gz"
norm_output=$(mktemp)
norm_stdout=$(mktemp)

# Run the command and capture both stdout and stderr
bcftools norm -m-any --force -a --atom-overlaps . -W -f "$fasta_file" "$vcf_file" -Oz -o "$normalized_vcf" 2> "$norm_output" 1> "$norm_stdout"
norm_exit_code=$?

# Check if there are any warnings or errors in the stderr output
if [[ $norm_exit_code -ne 0 ]]; then
    log_msg "Error: Failed to normalize the VCF."
    log_msg "bcftools norm error details: $(cat "$norm_output")"
    rm -f "$norm_output" "$norm_stdout"
    exit 1
fi

# Log any warnings (even if the command succeeded)
if grep -q "Warning" "$norm_output"; then
    log_msg "bcftools norm warnings: $(cat "$norm_output")"
fi

# Log the summary line from stdout (e.g., Lines total/split/joined/realigned/skipped)
if grep -q "Lines" "$norm_stdout"; then
    log_msg "bcftools norm summary: $(grep 'Lines' "$norm_stdout")"
fi

# Clean up the temporary files
rm -f "$norm_output" "$norm_stdout"

debug_msg "Normalized VCF written to: $normalized_vcf"

# Step 6: Apply filters to the normalized VCF
log_msg "Filtering the normalized VCF file..."

# Start building the pipeline command
pipeline_cmd="bcftools view $normalized_vcf | bcftools +fill-tags"

# Apply filters based on the inclusion and exclusion annotations
if [[ -f "$tmp_dir/merged_include_regions.bed.gz" ]]; then
    filter_cmd="bcftools filter -s NOT_IN_INCLUDE_REGION -m+ -e 'INFO/INCLUDE_REGION!=1'"
    pipeline_cmd="$pipeline_cmd | $filter_cmd"
    debug_msg "Applied filter for inclusion regions: $filter_cmd"
fi

if [[ -f "$tmp_dir/merged_exclude_regions.bed.gz" ]]; then
    filter_cmd="bcftools filter -s IN_EXCLUDE_REGION -m+ -e 'INFO/EXCLUDE_REGION=1'"
    pipeline_cmd="$pipeline_cmd | $filter_cmd"
    debug_msg "Applied filter for exclusion regions: $filter_cmd"
fi

# Apply inline filters
for filter in "${filters[@]}"; do
    IFS=" " read -r filter_name filter_action filter_expr <<< "$filter"
    filter_cmd="bcftools filter -m+ -s$filter_name -$filter_action '$filter_expr'"
    pipeline_cmd="$pipeline_cmd | $filter_cmd"
    debug_msg "Applied inline filter: $filter_cmd"
done

# Apply filters from file, if provided
if [[ -n "$filters_file" ]]; then
    while IFS=" " read -r filter_name filter_action filter_expr; do
        filter_expr=$(echo "$filter_expr" | tr -d '\r\n')
        filter_cmd="bcftools filter -m+ -s$filter_name -$filter_action '$filter_expr'"
        pipeline_cmd="$pipeline_cmd | $filter_cmd"
        debug_msg "Applied filter from file: $filter_cmd"
    done < "$filters_file"
fi

# Apply PASS filter if the --only-pass option is set
if $only_pass; then
    pass_filter_cmd="bcftools view -f PASS"
    pipeline_cmd="$pipeline_cmd | $pass_filter_cmd"
    debug_msg "Applied PASS filter: $pass_filter_cmd"
fi

# Determine the output type based on the output file extension
if [[ -n "$output_vcf" ]]; then
    if [[ "$output_vcf" == *.vcf.gz ]]; then
        output_type="z"  # Compressed VCF
    elif [[ "$output_vcf" == *.vcf ]]; then
        output_type="v"  # Uncompressed VCF
    else
        log_msg "Error: Unrecognized output file format for $output_vcf."
        exit 1
    fi
    pipeline_cmd="$pipeline_cmd | bcftools view -O$output_type -o $output_vcf"
else
    pipeline_cmd="$pipeline_cmd | bcftools view"
fi

# Print the final composed pipeline command in debug mode
debug_msg "Executing pipeline: $pipeline_cmd"

# Execute the composed pipeline command
eval "$pipeline_cmd"

# Check for errors
if [[ $? -ne 0 ]]; then
    log_msg "Error: Failed to filter the VCF."
    exit 1
fi

# Step 7: Generate stats file if the --generate-stats option is set and output_vcf is provided
if $generate_stats && [[ -n "$output_vcf" ]]; then
    stats_output="${output_vcf%.vcf.gz}.stats.txt"
    debug_msg "Generating stats file: $stats_output"
    bcftools stats "$output_vcf" > "$stats_output"
    if [[ $? -ne 0 ]]; then
        log_msg "Error: Failed to generate stats file."
        exit 1
    fi
    log_msg "Stats file saved to $stats_output"

    # If plotting is requested
    if $plot_stats; then
        log_msg "Plotting stats to $plot_output_dir"
        plot_output=$(mktemp)

        # Run the plot-vcfstats command and capture its output
        plot-vcfstats "$stats_output" -p "$plot_output_dir" > "$plot_output" 2>&1
        plot_exit_code=$?

        # Process and log the output from the plot-vcfstats command
        while IFS= read -r line; do
            log_msg "Plot-vcfstats output: $line"
        done < "$plot_output"

        if [[ $plot_exit_code -ne 0 ]]; then
            log_msg "Error: Failed to plot stats."
            rm -f "$plot_output"
            exit 1
        fi

        log_msg "Plots saved to $plot_output_dir"
        rm -f "$plot_output"
    fi
else
    debug_msg "Stats generation skipped (either --generate-stats was not set or no output file provided)."
fi

# Cleanup temporary files
if $cleanup; then
    debug_msg "Cleaning up temporary directory: $tmp_dir"
    rm -rf "$tmp_dir"
else
    debug_msg "Temporary directory not cleaned up: $tmp_dir"
fi

if [[ -n "$output_vcf" ]]; then
    log_msg "Filtered VCF saved to $output_vcf"
fi

# Disable debugging if it was enabled
if $debug; then
    set +x  # Disable command tracing
fi
