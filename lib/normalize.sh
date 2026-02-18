#!/bin/bash
# lib/normalize.sh — VCF normalization via bcftools norm
# Provides: normalize_vcf
# Requires: lib/logging.sh (log_msg, debug_msg)

[[ -n "${_LIB_NORMALIZE_LOADED:-}" ]] && return 0
readonly _LIB_NORMALIZE_LOADED=1

# normalize_vcf — normalize a VCF file using bcftools norm
# Usage: normalize_vcf <vcf_file> <fasta_file> <output_vcf> <tmp_dir>
# Does NOT use run_cmd — bcftools norm emits warnings on stderr even on success.
# Returns 1 on failure.
normalize_vcf() {
	local vcf_file="$1"
	local fasta_file="$2"
	local output_vcf="$3"
	local tmp_dir="$4"

	local norm_stderr norm_stdout
	norm_stderr=$(mktemp "${tmp_dir}/norm-stderr-XXXXXX")
	norm_stdout=$(mktemp "${tmp_dir}/norm-stdout-XXXXXX")

	bcftools norm -m-any --force -a --atom-overlaps . --write-index=tbi \
		-f "$fasta_file" "$vcf_file" \
		-Oz -o "$output_vcf" \
		2>"$norm_stderr" 1>"$norm_stdout" \
		|| {
			log_msg "Error: Failed to normalize the VCF."
			log_msg "bcftools norm error details: $(cat "$norm_stderr")"
			rm -f "$norm_stderr" "$norm_stdout"
			return 1
		}

	# Log any warnings (even if the command succeeded)
	if grep -q "Warning" "$norm_stderr"; then
		log_msg "bcftools norm warnings: $(cat "$norm_stderr")"
	fi

	# Log the summary line from stdout (e.g., Lines total/split/joined/realigned/skipped)
	if grep -q "Lines" "$norm_stdout"; then
		log_msg "bcftools norm summary: $(grep 'Lines' "$norm_stdout")"
	fi

	rm -f "$norm_stderr" "$norm_stdout"
	debug_msg "Normalized VCF written to: $output_vcf"
}
