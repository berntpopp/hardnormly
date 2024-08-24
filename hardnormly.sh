#!/bin/bash

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
tmp_dir="tmp"  # Default temporary directory
cleanup=true  # Default to cleaning up the temporary directory
slop=20  # Default slop value (in base pairs)
only_pass=false  # Option to filter only PASS variants
generate_stats=false  # Option to generate stats file

# Function to display help message
show_help() {
    echo "Usage: $0 -v <vcf_file> -f <fasta_file> [-o <output_vcf>] [options]"
    echo ""
    echo "Options:"
    echo "  -v, --vcf            Input VCF file (required)"
    echo "  -f, --fasta          Reference FASTA file for normalization (required)"
    echo "  -b, --include-bed    BED file(s) for inclusion (can specify multiple)"
    echo "  -e, --exclude-bed    BED file(s) for exclusion (can specify multiple)"
    echo "  -g, --genome         Genome file for slop operation (will be auto-generated if not provided)"
    echo "  --genome-build       Genome build to use for UCSC MySQL query (default: hg19)"
    echo "  --slop               Slop size for region padding (default: 20bp)"
    echo "  -o, --output         Output VCF file (if not specified, output will be sent to stdout)"
    echo "  --filters            Inline bcftools filter expression (can specify multiple; format: filter_name action expression)"
    echo "  --filters-file       File containing bcftools filter expressions (format: filter_name action expression)"
    echo "  --only-pass          Filter to retain only PASS variants"
    echo "  --generate-stats     Generate stats file from the output VCF"
    echo "  --tmp-dir            Temporary directory to use (default: tmp)"
    echo "  --no-cleanup         Do not clean up the temporary directory after execution"
    echo "  --debug              Enable debug mode (prints all executed commands and detailed messages)"
    echo "  -h, --help           Display this help message"
    exit 1
}

# Function to print debug messages
debug_msg() {
    if $debug; then
        echo "[DEBUG] $1"
    fi
}

# Parsing command-line arguments
while [[ "$#" -gt 0 ]]; do
    case "$1" in
        -v|--vcf) vcf_file="$2"; shift ;;
        -f|--fasta) fasta_file="$2"; shift ;;
        -b|--include-bed) include_bed_files+=("$2"); shift ;;
        -e|--exclude-bed) exclude_bed_files+=("$2"); shift ;;
        -g|--genome) genome_file="$2"; shift ;;
        --genome-build) genome_build="$2"; shift ;;
        --slop) slop="$2"; shift ;;
        -o|--output) output_vcf="$2"; shift ;;
        --filters) filters+=("$2"); shift ;;
        --filters-file) filters_file="$2"; shift ;;
        --only-pass) only_pass=true ;;
        --generate-stats) generate_stats=true ;;
        --tmp-dir) tmp_dir="$2"; shift ;;
        --no-cleanup) cleanup=false ;;
        --debug) debug=true ;;
        -h|--help) show_help ;;
        *) echo "Unknown parameter: $1"; show_help ;;
    esac
    shift
done

# Enable debugging if --debug flag is set
if $debug; then
    set -x  # Enable command tracing (prints every command)
fi

# Check required parameters
if [[ -z "$vcf_file" || -z "$fasta_file" ]]; then
    echo "Error: Missing required parameters."
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
    awk -v annot="$annotation" '{OFS="\t"; print $1, $2, $3, annot}' "$bed_file" | sort -k1,1 -k2,2n > "$output_file"
}

# Step 1: Create the genome file (if not provided)
if [[ -z "$genome_file" ]]; then
    genome_file="$tmp_dir/${genome_build}.genome"
    echo "Creating genome file: $genome_file for genome build: $genome_build"
    mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -e "select chrom, size from ${genome_build}.chromInfo" | grep -v "^chrom" | sed 's/chr//g' > "$genome_file"
    if [[ $? -ne 0 ]]; then
        echo "Error: Failed to create genome file from UCSC MySQL database."
        exit 1
    fi
    debug_msg "Genome file created: $genome_file"
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
    echo "Combining normalized exclusion BED files..."
    bedtools multiinter -i "${normalized_exclude_bed_files[@]}" | \
    awk '{OFS="\t"; print $1, $2, $3, "1"}' > "$tmp_dir/merged_exclude_regions.bed"
elif [[ "${#normalized_exclude_bed_files[@]}" -eq 1 ]]; then
    cp "${normalized_exclude_bed_files[0]}" "$tmp_dir/merged_exclude_regions.bed"
fi

# Combine normalized include BED files using bedtools multiinter
if [[ "${#normalized_include_bed_files[@]}" -gt 1 ]]; then
    echo "Combining normalized inclusion BED files..."
    bedtools multiinter -i "${normalized_include_bed_files[@]}" | \
    awk '{OFS="\t"; print $1, $2, $3, "1"}' > "$tmp_dir/merged_include_regions.bed"
elif [[ "${#normalized_include_bed_files[@]}" -eq 1 ]]; then
    cp "${normalized_include_bed_files[0]}" "$tmp_dir/merged_include_regions.bed"
fi

# Compress and index the merged BED files for inclusion
if [[ -f "$tmp_dir/merged_include_regions.bed" ]]; then
    bgzip -f "$tmp_dir/merged_include_regions.bed"  # Force overwrite
    tabix -p bed "$tmp_dir/merged_include_regions.bed.gz"
    debug_msg "Inclusion BED files normalized, merged, and indexed: $tmp_dir/merged_include_regions.bed.gz"
else
    debug_msg "No inclusion BED files provided; skipping inclusion annotation."
fi

# Compress and index the merged BED files for exclusion
if [[ -f "$tmp_dir/merged_exclude_regions.bed" ]]; then
    bgzip -f "$tmp_dir/merged_exclude_regions.bed"  # Force overwrite
    tabix -p bed "$tmp_dir/merged_exclude_regions.bed.gz"
    debug_msg "Exclusion BED files normalized, merged, and indexed: $tmp_dir/merged_exclude_regions.bed.gz"
else
    debug_msg "No exclusion BED files provided; skipping exclusion annotation."
fi

# Step 3: Create header files for INFO fields with the correct format
if [[ -f "$tmp_dir/merged_include_regions.bed.gz" ]]; then
    echo '##INFO=<ID=INCLUDE_REGION,Number=1,Type=Integer,Description="Included region">' > "$tmp_dir/include_regions.hdr"
    debug_msg "Created include_regions.hdr file"
fi
if [[ -f "$tmp_dir/merged_exclude_regions.bed.gz" ]]; then
    echo '##INFO=<ID=EXCLUDE_REGION,Number=1,Type=Integer,Description="Excluded region">' > "$tmp_dir/exclude_regions.hdr"
    debug_msg "Created exclude_regions.hdr file"
fi

# Step 4: Annotate the VCF file with the BED regions
echo "Annotating VCF with BED regions..."

# Annotate with inclusion regions if the file exists
if [[ -f "$tmp_dir/merged_include_regions.bed.gz" ]]; then
    bcftools annotate -a "$tmp_dir/merged_include_regions.bed.gz" -h "$tmp_dir/include_regions.hdr" -c CHROM,FROM,TO,INCLUDE_REGION "$vcf_file" -Oz -o "$tmp_dir/temp_include_annotated.vcf.gz"
    if [[ $? -ne 0 ]]; then
        echo "Error: Failed to annotate VCF with inclusion regions."
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
        echo "Error: Failed to annotate VCF with exclusion regions."
        exit 1
    fi
    vcf_file="$tmp_dir/temp_exclude_annotated.vcf.gz"
    debug_msg "Annotated VCF with exclusion regions: $vcf_file"
else
    debug_msg "Skipping annotation with exclusion regions because the file does not exist."
fi

# Step 5: Normalize the VCF file and write to an intermediate file
normalized_vcf="$tmp_dir/normalized.vcf.gz"
bcftools norm -m-any --force -a --atom-overlaps . -W -f "$fasta_file" "$vcf_file" -Oz -o "$normalized_vcf"
if [[ $? -ne 0 ]]; then
    echo "Error: Failed to normalize the VCF."
    exit 1
fi
debug_msg "Normalized VCF written to: $normalized_vcf"

# Step 6: Apply filters to the normalized VCF
echo "Filtering the normalized VCF file..."

# Start building the pipeline command
pipeline_cmd="bcftools view $normalized_vcf | bcftools +fill-tags"

# Apply filters based on the inclusion and exclusion annotations
if [[ -f "$tmp_dir/merged_include_regions.bed.gz" ]]; then
    filter_cmd="bcftools filter -s NOT_IN_INCLUDE_REGION -m+ -e 'INFO/INCLUDE_REGION!=1'"
    pipeline_cmd="$pipeline_cmd | $filter_cmd"
fi

if [[ -f "$tmp_dir/merged_exclude_regions.bed.gz" ]]; then
    filter_cmd="bcftools filter -s IN_EXCLUDE_REGION -m+ -e 'INFO/EXCLUDE_REGION=1'"
    pipeline_cmd="$pipeline_cmd | $filter_cmd"
fi

# Apply inline filters
for filter in "${filters[@]}"; do
    IFS=" " read -r filter_name filter_action filter_expr <<< "$filter"
    filter_cmd="bcftools filter -m+ -s$filter_name -$filter_action '$filter_expr'"
    pipeline_cmd="$pipeline_cmd | $filter_cmd"
done

# Apply filters from file, if provided
if [[ -n "$filters_file" ]]; then
    while IFS=" " read -r filter_name filter_action filter_expr; do
        filter_cmd="bcftools filter -m+ -s$filter_name -$filter_action '$filter_expr'"
        pipeline_cmd="$pipeline_cmd | $filter_cmd"
    done < "$filters_file"
fi

# Apply PASS filter if the --only-pass option is set
if $only_pass; then
    pass_filter_cmd="bcftools view -f PASS"
    pipeline_cmd="$pipeline_cmd | $pass_filter_cmd"
fi

# Determine the output type based on the output file extension
if [[ -n "$output_vcf" ]]; then
    if [[ "$output_vcf" == *.vcf.gz ]]; then
        output_type="z"  # Compressed VCF
    elif [[ "$output_vcf" == *.vcf ]]; then
        output_type="v"  # Uncompressed VCF
    else
        echo "Error: Unrecognized output file format for $output_vcf."
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
    echo "Error: Failed to filter the VCF."
    exit 1
fi

# Step 7: Generate stats file if the --generate-stats option is set
if $generate_stats; then
    stats_output="${output_vcf%.vcf.gz}.stats.txt"
    debug_msg "Generating stats file: $stats_output"
    bcftools stats "$output_vcf" > "$stats_output"
    if [[ $? -ne 0 ]]; then
        echo "Error: Failed to generate stats file."
        exit 1
    fi
    echo "Stats file saved to $stats_output"
fi

# Cleanup temporary files
if $cleanup; then
    debug_msg "Cleaning up temporary directory: $tmp_dir"
    rm -rf "$tmp_dir"
else
    debug_msg "Temporary directory not cleaned up: $tmp_dir"
fi

if [[ -n "$output_vcf" ]]; then
    echo "Filtered VCF saved to $output_vcf"
fi

# Disable debugging if it was enabled
if $debug; then
    set +x  # Disable command tracing
fi