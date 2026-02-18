#!/bin/bash
# lib/bed.sh — BED file normalization, merging, compression, and header creation
# Provides: normalize_bed, merge_include_beds, merge_exclude_beds,
#           compress_index_bed, create_header_file
# Requires: lib/logging.sh (log_msg, debug_msg, error_msg, run_cmd)

[[ -n "${_LIB_BED_LOADED:-}" ]] && return 0
readonly _LIB_BED_LOADED=1

# normalize_bed — annotate and sort a BED file
# Usage: normalize_bed <bed_file> <annotation> <output_file>
normalize_bed() {
	local bed_file="$1"
	local annotation="$2"
	local output_file="$3"
	debug_msg "Normalizing BED file: $bed_file with annotation: $annotation"
	awk -v annot="$annotation" '{OFS="\t"; print $1, $2, $3, annot}' "$bed_file" \
		| bedtools sort -i - >"$output_file" \
		|| {
			error_msg "normalize_bed: bedtools sort pipeline failed for $bed_file"
			return 1
		}
	debug_msg "Normalized BED file written to: $output_file"
}

# merge_include_beds — intersect (and slop) one or more normalized include BED files
# Usage: merge_include_beds <output_file> <slop> <genome_file> <bed_file>...
# If zero BED files are provided, output_file is not created.
merge_include_beds() {
	local output_file="$1"
	local slop="$2"
	local genome_file="$3"
	shift 3
	local bed_files=("$@")

	if [[ "${#bed_files[@]}" -gt 1 ]]; then
		log_msg "Intersecting and padding normalized inclusion BED files..."
		bedtools intersect -a "${bed_files[0]}" -b "${bed_files[@]:1}" \
			| bedtools sort -i - \
			| bedtools slop -b "$slop" -g "$genome_file" >"$output_file" \
			|| {
				error_msg "merge_include_beds: bedtools intersect/sort/slop pipeline failed"
				return 1
			}
		debug_msg "Combined inclusion regions written to: $output_file"
	elif [[ "${#bed_files[@]}" -eq 1 ]]; then
		log_msg "Padding single normalized inclusion BED file..."
		run_cmd bedtools slop -b "$slop" -g "$genome_file" -i "${bed_files[0]}" >"$output_file"
		debug_msg "Single inclusion BED file padded and written to: $output_file"
	fi
}

# merge_exclude_beds — merge one or more normalized exclude BED files via multiinter
# Usage: merge_exclude_beds <output_file> <bed_file>...
# If zero BED files are provided, output_file is not created.
merge_exclude_beds() {
	local output_file="$1"
	shift
	local bed_files=("$@")

	if [[ "${#bed_files[@]}" -gt 1 ]]; then
		log_msg "Combining normalized exclusion BED files..."
		bedtools multiinter -i "${bed_files[@]}" \
			| bedtools sort -i - \
			| awk '{OFS="\t"; print $1, $2, $3, "1"}' >"$output_file" \
			|| {
				error_msg "merge_exclude_beds: bedtools multiinter/sort pipeline failed"
				return 1
			}
		debug_msg "Combined exclusion regions written to: $output_file"
	elif [[ "${#bed_files[@]}" -eq 1 ]]; then
		cp "${bed_files[0]}" "$output_file"
		debug_msg "Single exclusion BED file copied to: $output_file"
	fi
}

# compress_index_bed — bgzip and tabix-index a BED file in-place
# Usage: compress_index_bed <bed_file>
# Produces <bed_file>.gz and <bed_file>.gz.tbi
compress_index_bed() {
	local bed_file="$1"
	run_cmd bgzip -f "$bed_file"
	run_cmd tabix -p bed "${bed_file}.gz"
}

# create_header_file — write a VCF INFO header line for a BED annotation field
# Usage: create_header_file <field_name> <description> <output_hdr_file>
create_header_file() {
	local field_name="$1"
	local description="$2"
	local output_hdr_file="$3"
	echo "##INFO=<ID=${field_name},Number=1,Type=Integer,Description=\"${description}\">" >"$output_hdr_file"
	debug_msg "Created header file: $output_hdr_file"
}
