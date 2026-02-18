#!/bin/bash
# lib/annotate.sh — VCF annotation with BED region INFO fields
# Provides: annotate_vcf_with_regions
# Requires: lib/logging.sh (run_cmd)

[[ -n "${_LIB_ANNOTATE_LOADED:-}" ]] && return 0
readonly _LIB_ANNOTATE_LOADED=1

# annotate_vcf_with_regions — annotate a VCF file with a BED region INFO field
# Usage: annotate_vcf_with_regions <vcf_file> <bed_gz> <header_file> <field_name> <output_vcf>
# Runs bcftools annotate via run_cmd; returns non-zero on failure.
annotate_vcf_with_regions() {
	local vcf_file="$1"
	local bed_gz="$2"
	local header_file="$3"
	local field_name="$4"
	local output_vcf="$5"
	run_cmd bcftools annotate \
		-a "$bed_gz" \
		-h "$header_file" \
		-c CHROM,FROM,TO,"$field_name" \
		"$vcf_file" \
		-Oz -o "$output_vcf"
}
