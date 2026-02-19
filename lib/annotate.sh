#!/bin/bash
# lib/annotate.sh — VCF annotation with BED region INFO fields
# Provides: preprocess_vcf_filters, annotate_vcf_with_regions, strip_vcf_annotations
# Requires: lib/logging.sh (run_cmd, log_msg)

[[ -n "${_LIB_ANNOTATE_LOADED:-}" ]] && return 0
readonly _LIB_ANNOTATE_LOADED=1

# preprocess_vcf_filters — detect and fix non-compliant FILTER tags in a VCF
# Handles: (1) FILTER values in data not defined in the ##FILTER header section
#           (2) FILTER IDs containing invalid characters (e.g., dots per VCF spec)
# Uses bgzip+awk to read VCFs without triggering bcftools/htslib FILTER validation errors.
# VCF spec (4.1+): FILTER IDs should be alphanumeric; defined in ##FILTER header lines.
# bcftools annotate (unlike bcftools view) converts FILTER IDs to integers internally
# (BCF representation) and exits 255 when a FILTER ID is not defined in the header.
# Usage: preprocess_vcf_filters <vcf_in> <vcf_out>
# Returns 0 always; <vcf_out> is only created if issues are detected.
preprocess_vcf_filters() {
	local vcf_in="$1"
	local vcf_out="$2"

	# Pass 1: detect issues via bgzip (bypasses bcftools htslib FILTER validation)
	# Collect defined FILTER IDs from header and used FILTER values from data.
	# Use 'flds' for the split() array and 'k' for for-in loops to avoid
	# gawk "attempt to use array in scalar context" error when the same variable
	# name is used as both an array (split target) and a scalar (for-in key).
	local detect_out
	detect_out=$(
		bgzip -c -d "$vcf_in" 2>/dev/null | awk '
			BEGIN { has_dots = 0 }
			/^##FILTER=</ {
				start = index($0, "ID=")
				if (start > 0) {
					rest = substr($0, start + 3)
					end1 = index(rest, ","); end2 = index(rest, ">")
					if (end1 == 0) end1 = length(rest) + 1
					if (end2 == 0) end2 = length(rest) + 1
					endc = (end1 < end2) ? end1 : end2
					id = substr(rest, 1, endc - 1)
					if (index(id, ".") > 0) has_dots = 1
					gsub(/\./, "_", id)
					defined[id] = 1
				}
				next
			}
			/^#/ { next }
			NF >= 7 && $7 != "." && $7 != "PASS" {
				n = split($7, flds, ";")
				for (i = 1; i <= n; i++) {
					if (flds[i] == "PASS" || flds[i] == ".") continue
					if (index(flds[i], ".") > 0) has_dots = 1
					fi = flds[i]; gsub(/\./, "_", fi)
					if (!(fi in defined)) missing[fi] = 1
				}
			}
			END {
				cnt = 0
				for (k in missing) cnt++
				if (cnt > 0) {
					printf "MISSING:"
					sep = ""
					for (k in missing) { printf "%s%s", sep, k; sep = "," }
					printf "\n"
				}
				if (has_dots) print "HAS_DOTS"
			}
		' 2>/dev/null
	)

	# No issues found — skip preprocessing
	if [[ -z "$detect_out" ]]; then
		return 0
	fi

	local missing_csv=""
	local has_dots=false
	while IFS= read -r line; do
		case "$line" in
			MISSING:*) missing_csv="${line#MISSING:}" ;;
			HAS_DOTS) has_dots=true ;;
		esac
	done <<<"$detect_out"

	[[ -n "$missing_csv" ]] && log_msg "Warning: FILTER(s) used in VCF data but undefined in header: ${missing_csv} — adding header definitions."
	[[ "$has_dots" == "true" ]] && log_msg "Warning: FILTER IDs with '.' detected (non-VCF-spec) — renaming '.' to '_' for bcftools compatibility."

	# Pass 2: write fixed VCF — add missing ##FILTER headers and rename invalid IDs.
	# Use 'flds' for split() array and 'k' for for-in loops (same variable naming
	# convention as Pass 1 to avoid gawk scalar/array context errors).
	bgzip -c -d "$vcf_in" | awk \
		-v missing_csv="$missing_csv" \
		'
		BEGIN {
			OFS = "\t"
			if (missing_csv != "") {
				n = split(missing_csv, mf, ",")
				for (i = 1; i <= n; i++) if (mf[i] != "") need_header[mf[i]] = 1
			}
		}
		/^##FILTER=</ {
			line = $0
			start = index(line, "ID=")
			if (start > 0) {
				rest = substr(line, start + 3)
				end1 = index(rest, ","); end2 = index(rest, ">")
				if (end1 == 0) end1 = length(rest) + 1
				if (end2 == 0) end2 = length(rest) + 1
				endc = (end1 < end2) ? end1 : end2
				id = substr(rest, 1, endc - 1)
				new_id = id; gsub(/\./, "_", new_id)
				if (id != new_id) {
					escaped = id; gsub(/\./, "\\.", escaped)
					sub("ID=" escaped, "ID=" new_id, line)
				}
			}
			print line; next
		}
		/^#CHROM/ {
			# Inject missing FILTER header lines immediately before the #CHROM column line
			for (k in need_header) {
				print "##FILTER=<ID=" k ",Description=\"Filter imported from source VCF\">"
			}
			print; next
		}
		/^#/ { print; next }
		{
			# Rename dots in FILTER column (column 7); use flds[] for split array
			if (NF >= 7 && $7 != "." && $7 != "PASS") {
				n = split($7, flds, ";")
				result = ""
				for (i = 1; i <= n; i++) {
					fi = flds[i]; gsub(/\./, "_", fi)
					result = (i == 1) ? fi : result ";" fi
				}
				$7 = result
			}
			print
		}
		' | bgzip -c >"$vcf_out" \
		|| {
			error_msg "preprocess_vcf_filters: pass 2 rewrite pipeline failed for $vcf_in"
			return 1
		}

	run_cmd tabix -p vcf "$vcf_out"
}

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

# strip_vcf_annotations — remove specified INFO fields from a VCF via bcftools annotate -x
# Usage: strip_vcf_annotations <vcf_file> <strip_list> <output_vcf>
# strip_list is comma-separated (e.g., "INFO/CSQ,INFO/ANN"), matching bcftools annotate -x syntax.
strip_vcf_annotations() {
	local vcf_file="$1"
	local strip_list="$2"
	local output_vcf="$3"
	run_cmd bcftools annotate -x "$strip_list" "$vcf_file" -Oz -o "$output_vcf"
}
