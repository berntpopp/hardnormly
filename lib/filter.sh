#!/bin/bash
# lib/filter.sh — VCF filter pipeline: fill-tags, sequential filters, output writing
# Provides: init_filter_pipeline, apply_filter_stages, write_filtered_output
# Requires: lib/logging.sh (log_msg, debug_msg)

[[ -n "${_LIB_FILTER_LOADED:-}" ]] && return 0
readonly _LIB_FILTER_LOADED=1

# init_filter_pipeline — initialize filter pipeline BCF with fill-tags applied
# Usage: init_filter_pipeline <input_vcf> <tmp_dir>
# Creates $tmp_dir/filter_current.bcf; returns 1 on failure.
init_filter_pipeline() {
	local input_vcf="$1"
	local tmp_dir="$2"

	bcftools view "$input_vcf" \
		| bcftools +fill-tags -Ob -o "${tmp_dir}/filter_current.bcf" \
		|| {
			log_msg "Error: Failed to initialize filter pipeline."
			return 1
		}
}

# apply_filter_stages — apply filter stages sequentially via BCF temp files
# Usage: apply_filter_stages <tmp_dir> <stage>...
# Each stage is "name|action|expression" (pipe-delimited).
# Returns 1 if any filter stage fails.
apply_filter_stages() {
	local tmp_dir="$1"
	shift
	local stages=("$@")

	local stage stage_name stage_action stage_expr
	for stage in "${stages[@]}"; do
		stage_name=""
		stage_action=""
		stage_expr=""
		IFS="|" read -r stage_name stage_action stage_expr <<<"$stage"
		debug_msg "Applying filter: $stage_name ($stage_action) '$stage_expr'"
		bcftools filter -m+ -s "$stage_name" "-${stage_action}" "$stage_expr" \
			"${tmp_dir}/filter_current.bcf" -Ob -o "${tmp_dir}/filter_next.bcf" \
			|| {
				log_msg "Error: Filter '$stage_name' failed."
				return 1
			}
		mv "${tmp_dir}/filter_next.bcf" "${tmp_dir}/filter_current.bcf"
	done
}

# write_filtered_output — detect output format, optionally apply PASS filter, write final output
# Usage: write_filtered_output <tmp_dir> <output_vcf> <only_pass> <auto_index>
# Returns 1 on failure.
write_filtered_output() {
	local tmp_dir="$1"
	local output_vcf="$2"
	local only_pass="$3"
	local auto_index="$4"

	local output_args=()
	if [[ -n "$output_vcf" ]]; then
		if [[ "$output_vcf" == *.vcf.gz ]]; then
			output_args+=("-Oz")
			if [[ "$auto_index" == "true" ]]; then
				output_args+=("--write-index=tbi")
				debug_msg "Auto-index enabled for compressed output."
			fi
		elif [[ "$output_vcf" == *.vcf ]]; then
			output_args+=("-Ov")
		else
			log_msg "Error: Unrecognized output file format for $output_vcf."
			return 1
		fi
		output_args+=("-o" "$output_vcf")
	fi

	if [[ "$only_pass" == "true" ]]; then
		bcftools view -f PASS "${output_args[@]}" "${tmp_dir}/filter_current.bcf" \
			|| {
				log_msg "Error: PASS filter failed."
				return 1
			}
	else
		bcftools view "${output_args[@]}" "${tmp_dir}/filter_current.bcf" \
			|| {
				log_msg "Error: Failed to write output."
				return 1
			}
	fi
}
