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

	# Drop records on contigs absent from the FASTA (e.g. hs38DH decoy/HLA/ALT
	# contigs vs a no-alt FASTA): bcftools norm aborts on them (issue #27).
	# Targets file (CHROM, 1, length) avoids region parsing of names with ':'.
	local contigs_file="${tmp_dir}/fasta_contigs.tsv"
	if [[ -f "${fasta_file}.fai" ]]; then
		awk -v OFS='\t' '{print $1, 1, $2}' "${fasta_file}.fai" >"$contigs_file"
	else
		# No index: take names from headers; END is irrelevant to whole-contig matching
		gzip -cdf "$fasta_file" | awk -v OFS='\t' '/^>/{print substr($1, 2), 1, 2147483647}' >"$contigs_file"
	fi
	local n_dropped
	n_dropped=$(bcftools view -H -T "^${contigs_file}" "$vcf_file" | wc -l) || {
		log_msg "Error: Failed to check VCF contigs against the FASTA."
		return 1
	}
	if [[ "$n_dropped" -gt 0 ]]; then
		log_msg "Dropping ${n_dropped} records on contigs not present in the FASTA."
		bcftools view -T "$contigs_file" "$vcf_file" -Oz -o "${tmp_dir}/fasta_contigs_only.vcf.gz" || {
			log_msg "Error: Failed to restrict the VCF to FASTA contigs."
			return 1
		}
		vcf_file="${tmp_dir}/fasta_contigs_only.vcf.gz"
	fi

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
