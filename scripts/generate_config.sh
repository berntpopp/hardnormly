#!/usr/bin/env bash
# scripts/generate_config.sh — Interactive config generator for the hardnormly pipeline family
#
# Generates config.yaml and samples.tsv templates for one or more of:
#   alignment   — sm-alignment  (FASTQ → analysis-ready BAM)
#   calling     — sm-vcf-calling (BAM → filtered VCF)
#   hardnormly  — hardnormly workflow (VCF → normalized/hard-filtered VCF)
#
# Usage:
#   scripts/generate_config.sh [OPTIONS]
#
# Options:
#   -p, --pipeline PIPELINE  Pipeline to configure (repeatable):
#                              alignment   — sm-alignment only
#                              calling     — sm-vcf-calling only
#                              hardnormly  — hardnormly workflow only
#                              all         — all three (default)
#   -o, --output-dir DIR     Write generated files here (default: current directory)
#   -n, --non-interactive    Accept all defaults without prompting
#   -y, --yes                Skip write confirmation (write files without asking)
#   -h, --help               Show this help and exit
#
# Examples:
#   # Interactive — configure all three pipelines, write to ./configs/
#   scripts/generate_config.sh -o configs/
#
#   # Configure only the hardnormly workflow config
#   scripts/generate_config.sh -p hardnormly -o config/
#
#   # Non-interactive — write templates with all defaults, no prompts
#   scripts/generate_config.sh -p hardnormly -n -y -o config/
#
#   # Configure alignment and calling, skip confirmation
#   scripts/generate_config.sh -p alignment -p calling -y -o project-configs/

set -euo pipefail
shopt -s inherit_errexit

# ---------------------------------------------------------------------------
# Globals
# ---------------------------------------------------------------------------

NON_INTERACTIVE="false"
SKIP_CONFIRM="false"
OUTPUT_DIR="."
declare -a SELECTED_PIPELINES=()

# Arrays declared globally for nameref to work across functions
declare -a AL_KNOWN_SITES=()
declare -a HN_INCLUDE_BEDS=()
declare -a HN_EXCLUDE_BEDS=()

# ---------------------------------------------------------------------------
# Terminal colors (disabled when not a TTY)
# ---------------------------------------------------------------------------

if [[ -t 1 ]]; then
	BOLD='\033[1m'
	DIM='\033[2m'
	BLUE='\033[0;34m'
	CYAN='\033[0;36m'
	GREEN='\033[0;32m'
	YELLOW='\033[1;33m'
	MAGENTA='\033[0;35m'
	RED='\033[0;31m'
	RESET='\033[0m'
else
	BOLD='' DIM='' BLUE='' CYAN='' GREEN='' YELLOW='' MAGENTA='' RED='' RESET=''
fi

# ---------------------------------------------------------------------------
# Output helpers
# ---------------------------------------------------------------------------

die() {
	printf '%bERROR: %s%b\n' "$RED" "$1" "$RESET" >&2
	exit 1
}
info() { printf '%b%s%b\n' "$CYAN" "$1" "$RESET"; }
success() { printf '%b  ✔  %s%b\n' "$GREEN" "$1" "$RESET"; }
warn() { printf '%b  ⚠  %s%b\n' "$YELLOW" "$1" "$RESET"; }
skipped() { printf '%b  –  %s%b\n' "$DIM" "$1" "$RESET"; }

# Pipeline banner — large section heading
header() {
	local width=60
	printf '\n'
	printf '%b' "$BOLD$BLUE"
	printf '%.0s─' $(seq 1 $width)
	printf '\n  %s\n' "$1"
	printf '%.0s─' $(seq 1 $width)
	printf '%b\n\n' "$RESET"
}

# Sub-section heading within a pipeline
section() { printf '\n%b▸ %s%b\n' "$BOLD$MAGENTA" "$1" "$RESET"; }

# Hint line below a prompt for extra context
hint() { printf '  %b%s%b\n' "$DIM" "$1" "$RESET"; }

# Warn if a given path doesn't exist (non-fatal)
validate_path() {
	local label="$1" path="$2"
	[[ -z "$path" ]] && return
	[[ -e "$path" ]] || warn "${label}: file not found — ${path}"
}

# ---------------------------------------------------------------------------
# Prompt helpers
# ---------------------------------------------------------------------------

# prompt NAMEREF QUESTION DEFAULT [HINT]
# Prompts for a string value; uses DEFAULT if empty or non-interactive.
prompt() {
	local -n _pv="$1"
	local question="$2"
	local default="${3:-}"
	local help_text="${4:-}"

	if [[ -n "$help_text" ]]; then hint "$help_text"; fi

	if [[ "$NON_INTERACTIVE" == "true" ]]; then
		_pv="$default"
		if [[ -n "$default" ]]; then
			skipped "${question}: ${default}"
		else
			skipped "${question}: (empty)"
		fi
		return
	fi

	local hint_str=""
	[[ -n "$default" ]] && hint_str=" ${DIM}[${default}]${RESET}"
	printf '%b%s%b%s: ' "$BOLD" "$question" "$RESET" "$hint_str"
	IFS= read -r _pv || true
	[[ -z "$_pv" ]] && _pv="$default"
}

# prompt_bool NAMEREF QUESTION DEFAULT(y|n) [HINT]
# Prompts for yes/no; stores "true" or "false".
prompt_bool() {
	local -n _pb="$1"
	local question="$2"
	local default="${3:-n}"
	local help_text="${4:-}"

	if [[ -n "$help_text" ]]; then hint "$help_text"; fi

	if [[ "$NON_INTERACTIVE" == "true" ]]; then
		case "$default" in
			y | yes | true | 1) _pb="true" ;;
			*) _pb="false" ;;
		esac
		if [[ "$_pb" == "true" ]]; then
			skipped "${question}: yes"
		else
			skipped "${question}: no"
		fi
		return
	fi

	local hint_str="[y/N]"
	case "$default" in y | yes | true | 1) hint_str="[Y/n]" ;; esac

	printf '%b%s%b %s: ' "$BOLD" "$question" "$RESET" "$hint_str"
	local answer=""
	IFS= read -r answer || true
	answer="${answer,,}"
	case "$answer" in
		y | yes) _pb="true" ;;
		n | no) _pb="false" ;;
		*)
			case "$default" in
				y | yes | true | 1) _pb="true" ;;
				*) _pb="false" ;;
			esac
			;;
	esac
}

# prompt_choice NAMEREF QUESTION "value:description" ...
# Numbered menu; stores the chosen value. First option is the default.
prompt_choice() {
	local -n _pc="$1"
	local question="$2"
	shift 2
	local -a options=("$@")

	if [[ "$NON_INTERACTIVE" == "true" ]]; then
		_pc="${options[0]%%:*}"
		skipped "${question}: ${_pc} (default)"
		return
	fi

	printf '\n%b%s%b\n' "$BOLD" "$question" "$RESET"
	local i=1
	for opt in "${options[@]}"; do
		printf '  %b%d)%b %-16s %b%s%b\n' \
			"$CYAN" "$i" "$RESET" \
			"${opt%%:*}" \
			"$DIM" "${opt#*:}" "$RESET"
		((i++))
	done

	local answer="" default_val="${options[0]%%:*}"
	printf '%bChoice%b [1-%d, default: %s]: ' "$BOLD" "$RESET" "${#options[@]}" "$default_val"
	IFS= read -r answer || true

	if [[ -z "$answer" ]]; then
		_pc="$default_val"
	elif [[ "$answer" =~ ^[0-9]+$ ]] && ((answer >= 1 && answer <= ${#options[@]})); then
		_pc="${options[$((answer - 1))]%%:*}"
	else
		_pc="$answer"
	fi
}

# prompt_array ARRAY_NAME QUESTION EXAMPLE [HINT]
# Collects multiple values; empty line finishes. Assigns to the named global array.
prompt_array() {
	local arr_name="$1"
	local question="$2"
	local example="${3:-}"
	local help_text="${4:-}"

	local -n _pa_ref="$arr_name"
	_pa_ref=()

	if [[ -n "$help_text" ]]; then hint "$help_text"; fi

	if [[ "$NON_INTERACTIVE" == "true" ]]; then
		[[ -n "$example" ]] && _pa_ref=("$example")
		skipped "${question}: (empty — fill in manually)"
		return
	fi

	printf '\n%b%s%b\n' "$BOLD" "$question" "$RESET"
	printf '  %bOne path per line. Empty line to finish.%b\n' "$DIM" "$RESET"
	[[ -n "$example" ]] && printf '  %bExample: %s%b\n' "$DIM" "$example" "$RESET"

	local item=""
	while true; do
		printf '  %b›%b ' "$CYAN" "$RESET"
		IFS= read -r item || true
		[[ -z "$item" ]] && break
		_pa_ref+=("$item")
	done
}

# prompt_optional NAMEREF QUESTION HINT
# Like prompt() but explicitly marks the field as optional.
prompt_optional() {
	local -n _po="$1"
	local question="$2"
	local help_text="${3:-}"
	prompt _po "${question} ${DIM}(optional)${RESET}" "" "$help_text"
}

# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

show_help() {
	sed -n '/^# Usage:/,/^[^#]/{ /^[^#]/d; s/^# \{0,1\}//; p }' "${BASH_SOURCE[0]}"
	exit 0
}

parse_args() {
	while [[ $# -gt 0 ]]; do
		case "$1" in
			-p | --pipeline)
				[[ -z "${2:-}" ]] && die "--pipeline requires a value"
				SELECTED_PIPELINES+=("$2")
				shift 2
				;;
			-o | --output-dir)
				[[ -z "${2:-}" ]] && die "--output-dir requires a value"
				OUTPUT_DIR="$2"
				shift 2
				;;
			-n | --non-interactive)
				NON_INTERACTIVE="true"
				shift
				;;
			-y | --yes)
				SKIP_CONFIRM="true"
				shift
				;;
			-h | --help)
				show_help
				;;
			*)
				die "Unknown option: $1 (use -h for help)"
				;;
		esac
	done

	[[ "${#SELECTED_PIPELINES[@]}" -eq 0 ]] && SELECTED_PIPELINES=(all)

	# Expand "all" and deduplicate while preserving order
	local -A _seen=()
	local -a _expanded=()
	for p in "${SELECTED_PIPELINES[@]}"; do
		if [[ "$p" == "all" ]]; then
			for sub in alignment calling hardnormly; do
				[[ -z "${_seen[$sub]+x}" ]] && {
					_expanded+=("$sub")
					_seen["$sub"]=1
				}
			done
		else
			case "$p" in
				alignment | calling | hardnormly)
					[[ -z "${_seen[$p]+x}" ]] && {
						_expanded+=("$p")
						_seen["$p"]=1
					}
					;;
				*)
					die "Unknown pipeline: '$p' (valid: alignment, calling, hardnormly, all)"
					;;
			esac
		fi
	done
	SELECTED_PIPELINES=("${_expanded[@]}")
}

# ---------------------------------------------------------------------------
# YAML helpers
# ---------------------------------------------------------------------------

# Always double-quote string values to avoid ambiguity
yaml_str() {
	local v="${1:-}"
	v="${v//\\/\\\\}"
	v="${v//\"/\\\"}"
	printf '"%s"' "$v"
}

yaml_bool() { [[ "${1:-}" == "true" ]] && printf 'true' || printf 'false'; }

# ===========================================================================
# sm-alignment
# ===========================================================================

configure_alignment() {
	header "sm-alignment  ·  FASTQ → Analysis-ready BAM"

	# --- Reference --------------------------------------------------------
	section "Reference"
	prompt AL_GENOME "Reference genome FASTA (uncompressed)" "" \
		"Absolute or relative path to the reference .fa / .fasta file"
	validate_path "ref.genome" "$AL_GENOME"

	prompt_choice AL_BUILD "Genome build" \
		"GRCh38:GRCh38 / hg38  (recommended for new projects)" \
		"GRCh37:GRCh37 / hg19"

	prompt_array AL_KNOWN_SITES \
		"Known variant sites VCFs for BQSR" \
		"/path/to/dbsnp.vcf.gz" \
		"dbSNP, Mills & 1000G gold standard indels, 1000G phase 1 snps — one path per line"

	# --- Paths ------------------------------------------------------------
	section "Paths"
	prompt AL_FASTQ_FOLDER "Input FASTQ directory" "" \
		"Root directory containing all FASTQ files (subfolders supported via samples.tsv)"
	validate_path "paths.fastq_folder" "$AL_FASTQ_FOLDER"

	prompt AL_OUTPUT_FOLDER "Output directory for BAMs" "results/alignment"
	prompt AL_SAMPLES_TSV "Samples metadata TSV path" "config/samples.tsv" \
		"Will be created from the template written by this script"
	prompt AL_LOG_SUBDIR "Log subdirectory name" "logs"

	# --- FASTQ naming & read group ----------------------------------------
	section "FASTQ naming & read groups"
	prompt AL_R1_SUFFIX "R1 filename suffix" "_R1_001.fastq.gz" \
		"Everything after the sample basename that identifies read 1"
	prompt AL_R2_SUFFIX "R2 filename suffix" "_R2_001.fastq.gz"
	prompt AL_PLATFORM "Sequencing platform for PL read group tag" "ILLUMINA" \
		"Standard values: ILLUMINA  ONT  PACBIO  ELEMENT"

	# --- Trimming ---------------------------------------------------------
	section "Adapter trimming  (BBDuk)"
	prompt_bool AL_TRIM_ENABLED "Enable BBDuk adapter trimming?" "n" \
		"Trims adapters and low-quality bases before alignment"

	# --- QC ---------------------------------------------------------------
	section "Quality control"
	prompt_bool AL_QC_ENABLED "Enable QC rules (master switch)?" "y"
	prompt_bool AL_QC_FASTQC "  FastQC on raw FASTQs?" "y"
	prompt_bool AL_QC_SAMTOOLS_STATS "  samtools stats on final BAM?" "y"
	prompt_bool AL_QC_SAMTOOLS_FLAGSTAT "  samtools flagstat on final BAM?" "y"
	prompt_bool AL_QC_PICARD "  Picard CollectMultipleMetrics?" "y"

	prompt_bool AL_QC_QUALIMAP "  Qualimap bamqc (high memory)?" "n" \
		"Detailed coverage analysis — requires extra memory (~16 GB)"
	if [[ "$AL_QC_QUALIMAP" == "true" ]]; then
		prompt_optional AL_QUALIMAP_BED "  Qualimap feature BED" \
			"BED file with target regions for exome-aware coverage stats"
		validate_path "qc.qualimap_feature_file" "$AL_QUALIMAP_BED"
	else
		AL_QUALIMAP_BED=""
	fi

	# --- Region subsetting ------------------------------------------------
	section "BAM subsetting  (optional)"
	prompt_optional AL_SUBSET_BED "Subset BAM to region BED" \
		"Extract only reads overlapping this BED after final BAM is produced; leave empty to disable"
	validate_path "subset.bed_file" "$AL_SUBSET_BED"

	# --- Processing -------------------------------------------------------
	section "Processing"
	prompt AL_COMPRESSION "BAM compression level (1–9)" "6" \
		"Higher = smaller files, slower write; 6 is a good default"
}

write_alignment_config() {
	local outfile="${OUTPUT_DIR}/sm-alignment-config.yaml"
	mkdir -p "$OUTPUT_DIR"

	{
		printf '# sm-alignment config.yaml\n'
		printf '# Generated by hardnormly scripts/generate_config.sh\n'
		printf '# Copy to:  sm-alignment/config/config.yaml\n\n'

		printf 'ref:\n'
		printf '  genome: %s\n' "$(yaml_str "$AL_GENOME")"
		printf '  build: %s\n' "$(yaml_str "$AL_BUILD")"
		printf '  known_sites:\n'
		if [[ "${#AL_KNOWN_SITES[@]}" -gt 0 && -n "${AL_KNOWN_SITES[0]}" ]]; then
			for s in "${AL_KNOWN_SITES[@]}"; do printf '    - %s\n' "$(yaml_str "$s")"; done
		else
			printf '    []  # TODO: add dbSNP / Mills / 1000G VCF paths\n'
		fi

		printf '\npaths:\n'
		printf '  samples: %s\n' "$(yaml_str "$AL_SAMPLES_TSV")"
		printf '  fastq_folder: %s\n' "$(yaml_str "$AL_FASTQ_FOLDER")"
		printf '  output_folder: %s\n' "$(yaml_str "$AL_OUTPUT_FOLDER")"
		printf '  log_subdir: %s\n' "$(yaml_str "$AL_LOG_SUBDIR")"

		printf '\nread_group:\n'
		printf '  platform: %s\n' "$(yaml_str "$AL_PLATFORM")"

		printf '\nfastq:\n'
		printf '  r1_suffix: %s\n' "$(yaml_str "$AL_R1_SUFFIX")"
		printf '  r2_suffix: %s\n' "$(yaml_str "$AL_R2_SUFFIX")"

		printf '\ntrimming:\n'
		printf '  enabled: %s\n' "$(yaml_bool "$AL_TRIM_ENABLED")"
		if [[ "$AL_TRIM_ENABLED" == "true" ]]; then
			printf '  # BBDuk adapter trimming defaults (adjust as needed)\n'
			printf '  bbduk_ref: "adapters,artifacts"\n'
			printf '  ktrim: "r"\n'
			printf '  k: 23\n'
			printf '  mink: 11\n'
			printf '  hdist: 1\n'
			printf '  tpe: "t"\n'
			printf '  tbo: "t"\n'
			printf '  ftl: 5\n'
			printf '  trimpolyg: 3\n'
			printf '  trimpolya: 3\n'
			printf '  qtrim: "t"\n'
			printf '  trimq: 10\n'
			printf '  ziplevel: 5\n'
		fi

		printf '\nqc:\n'
		printf '  enabled: %s\n' "$(yaml_bool "$AL_QC_ENABLED")"
		printf '  fastqc: %s\n' "$(yaml_bool "$AL_QC_FASTQC")"
		printf '  samtools_stats: %s\n' "$(yaml_bool "$AL_QC_SAMTOOLS_STATS")"
		printf '  samtools_flagstat: %s\n' "$(yaml_bool "$AL_QC_SAMTOOLS_FLAGSTAT")"
		printf '  picard_collect_metrics: %s\n' "$(yaml_bool "$AL_QC_PICARD")"
		printf '  qualimap: %s\n' "$(yaml_bool "$AL_QC_QUALIMAP")"
		if [[ -n "$AL_QUALIMAP_BED" ]]; then
			printf '  qualimap_feature_file: %s\n' "$(yaml_str "$AL_QUALIMAP_BED")"
		else
			printf '  qualimap_feature_file: ""  # BED for targeted coverage stats\n'
		fi

		printf '\nsubset:\n'
		if [[ -n "$AL_SUBSET_BED" ]]; then
			printf '  bed_file: %s\n' "$(yaml_str "$AL_SUBSET_BED")"
		else
			printf '  bed_file: ""  # Leave empty to disable BAM subsetting\n'
		fi

		printf '\nprocessing:\n'
		printf '  compression_level: %s\n' "$AL_COMPRESSION"

		printf '\n# Advanced params (edit as needed)\n'
		printf 'params:\n'
		printf '  bwa_mem:\n'
		printf '    extra: ""\n'
		printf '  samtools:\n'
		printf '    sort_extra: ""\n'
		printf '    merge_extra: ""\n'
		printf '  gatk:\n'
		printf '    MarkDuplicates: "--CREATE_INDEX true --VALIDATION_STRINGENCY SILENT"\n'
		printf '    BaseRecalibrator: ""\n'
		printf '    ApplyBQSR: ""\n'
	} >"$outfile"

	success "sm-alignment-config.yaml → ${outfile}"
}

write_alignment_samples() {
	local outfile="${OUTPUT_DIR}/sm-alignment-samples.tsv"

	{
		printf 'fastq_files_basename\tlane\tproject_sample\tmdc_project\tsubfolder\n'
		printf '# Each row = one FASTQ pair. Rows sharing project_sample are merged into one BAM.\n'
		printf '# fastq_files_basename: prefix before %s\n' "$AL_R1_SUFFIX"
		printf '# lane:                 e.g. L001 — used in read group @RG:ID\n'
		printf '# project_sample:       logical sample name — all rows merged into one BAM\n'
		printf '# mdc_project:          project ID for @RG:PU tag\n'
		printf '# subfolder:            subdirectory under fastq_folder (leave empty if flat)\n'
		printf '#\n'
		printf '# Example — two lanes of the same sample:\n'
		printf '# Sample1_S1_L001\tL001\tSample1\tProjectX\t\n'
		printf '# Sample1_S1_L002\tL002\tSample1\tProjectX\t\n'
	} >"$outfile"

	success "sm-alignment-samples.tsv   → ${outfile}"
}

# ===========================================================================
# sm-vcf-calling
# ===========================================================================

configure_calling() {
	header "sm-vcf-calling  ·  BAM → Filtered VCF"

	# --- Caller -----------------------------------------------------------
	section "Variant caller"
	prompt_choice CALL_CALLER "Which caller to run?" \
		"mutect2:GATK Mutect2  (somatic, tumor-normal, tumor-only, germline)" \
		"freebayes:FreeBayes  (germline)" \
		"all:Both Mutect2 and FreeBayes"

	# --- Reference --------------------------------------------------------
	section "Reference"
	prompt CALL_GENOME "Reference genome FASTA" "" \
		"Same reference used for alignment"
	validate_path "ref.genome" "$CALL_GENOME"

	prompt_choice CALL_BUILD "Genome build" \
		"GRCh38:GRCh38 / hg38  (recommended for new projects)" \
		"GRCh37:GRCh37 / hg19"

	# --- GATK resources ---------------------------------------------------
	if [[ "$CALL_CALLER" == "mutect2" || "$CALL_CALLER" == "all" ]]; then
		section "GATK resources  (required for Mutect2)"
		info "  These VCFs must match your reference genome build"
		prompt CALL_PON "Panel of Normals VCF" "" \
			"Pre-built PoN VCF from GATK best practices"
		validate_path "gatk_resources.panel_of_normals" "$CALL_PON"

		prompt CALL_GNOMAD "gnomAD allele frequency VCF" "" \
			"af-only-gnomad.*.vcf.gz from GATK resource bundle"
		validate_path "gatk_resources.af_only_gnomad" "$CALL_GNOMAD"

		prompt CALL_GNOMAD_CB "gnomAD common biallelic VCF" "" \
			"small_exac_common_3.*.vcf.gz from GATK resource bundle"
		validate_path "gatk_resources.common_biallelic_gnomad" "$CALL_GNOMAD_CB"
	else
		CALL_PON="" CALL_GNOMAD="" CALL_GNOMAD_CB=""
	fi

	# --- Paths ------------------------------------------------------------
	section "Paths"
	prompt CALL_BAM_FOLDER "Input BAM directory" "" \
		"Directory containing the analysis-ready BAMs from sm-alignment"
	validate_path "paths.bam_folder" "$CALL_BAM_FOLDER"

	prompt CALL_OUTPUT_FOLDER "Output directory for VCFs" "results/calling"
	prompt CALL_SAMPLES_TSV "Samples metadata TSV path" "config/samples.tsv"
	prompt CALL_LOG_SUBDIR "Log subdirectory name" "logs"
	prompt CALL_BAM_EXT "BAM file extension (suffix)" ".merged.dedup.bqsr.bam" \
		"Must match bam.final_suffix in sm-alignment config"

	# --- Scatter ----------------------------------------------------------
	section "Scatter / parallelisation"
	prompt_choice CALL_SCATTER "Scatter strategy" \
		"chromosome:Per-chromosome  (fast, recommended)" \
		"interval:Fixed-size intervals  (set count below)" \
		"none:No scatter  (single job per sample — slow for WGS)"

	if [[ "$CALL_SCATTER" == "interval" ]]; then
		prompt CALL_SCATTER_COUNT "Number of scatter intervals" "400" \
			"More intervals = more parallelism; 200–500 works well for WES"
	else
		CALL_SCATTER_COUNT="400"
	fi

	# --- Mutect2 advanced -------------------------------------------------
	if [[ "$CALL_CALLER" == "mutect2" || "$CALL_CALLER" == "all" ]]; then
		section "Mutect2 advanced params"
		hint "Defaults are optimised for somatic calling and PureCN compatibility"
		prompt_bool CALL_M2_GERMLINE "  genotype_germline_sites (required for PureCN)?" "y"
		prompt_bool CALL_M2_PON_SITES "  genotype_pon_sites (required for PureCN)?" "y"
		prompt_bool CALL_M2_NO_ADAPTIVE_PRUNING "  disable_adaptive_pruning (prevents low-AF dropout)?" "y"
		prompt_optional CALL_M2_EXTRA "  Additional Mutect2 flags" \
			"Passed verbatim to mutect2 — e.g. '--max-reads-per-alignment-start 0'"
	else
		CALL_M2_GERMLINE="true"
		CALL_M2_PON_SITES="true"
		CALL_M2_NO_ADAPTIVE_PRUNING="true"
		CALL_M2_EXTRA=""
	fi

	# --- FreeBayes advanced -----------------------------------------------
	if [[ "$CALL_CALLER" == "freebayes" || "$CALL_CALLER" == "all" ]]; then
		section "FreeBayes advanced params"
		prompt CALL_FB_EXTRA "FreeBayes extra flags" \
			"--min-coverage 20 --limit-coverage 500 --use-best-n-alleles 4 --standard-filters" \
			"Passed verbatim to freebayes"
	else
		CALL_FB_EXTRA="--min-coverage 20 --limit-coverage 500 --use-best-n-alleles 4 --standard-filters"
	fi

	# --- PureCN -----------------------------------------------------------
	section "PureCN copy number analysis  (optional)"
	prompt_bool CALL_PURECN "Enable PureCN?" "n" \
		"Tumour purity, ploidy, and copy number from Mutect2 output"

	if [[ "$CALL_PURECN" == "true" ]]; then
		prompt_choice CALL_PURECN_GENOME "PureCN genome identifier" \
			"hg38:GRCh38 / hg38" \
			"hg19:GRCh37 / hg19"
		prompt CALL_PURECN_BED "Capture bait BED file" "" \
			"Required: BED with target region coordinates"
		validate_path "purecn.intervals_bed" "$CALL_PURECN_BED"

		prompt_optional CALL_PURECN_NORMALDB "Pre-built normalDB.rds" \
			"Leave empty to build from matched normals in this run"
		validate_path "purecn.normaldb" "$CALL_PURECN_NORMALDB"

		prompt_optional CALL_PURECN_MAPPING_BIAS "Pre-built mapping_bias.rds" \
			"Optional — generated alongside normalDB"
		validate_path "purecn.mapping_bias" "$CALL_PURECN_MAPPING_BIAS"

		prompt_optional CALL_PURECN_SNP_BLACKLIST "SNP blacklist BED" \
			"Optional simple-repeat BED for SNP filtering"
		prompt_optional CALL_PURECN_EXTRA "Extra PureCN.R arguments"

		prompt CALL_PURECN_SEED "PureCN random seed" "123" \
			"Fix for reproducible results"
		prompt_bool CALL_PURECN_POSTOPTIMIZE "Run PureCN post-optimisation?" "y"
	else
		CALL_PURECN_GENOME="hg38"
		CALL_PURECN_BED=""
		CALL_PURECN_NORMALDB=""
		CALL_PURECN_MAPPING_BIAS=""
		CALL_PURECN_SNP_BLACKLIST=""
		CALL_PURECN_EXTRA=""
		CALL_PURECN_SEED="123"
		CALL_PURECN_POSTOPTIMIZE="true"
	fi
}

write_calling_config() {
	local outfile="${OUTPUT_DIR}/sm-calling-config.yaml"
	mkdir -p "$OUTPUT_DIR"

	{
		printf '# sm-vcf-calling config.yaml\n'
		printf '# Generated by hardnormly scripts/generate_config.sh\n'
		printf '# Copy to:  sm-vcf-calling/config/config.yaml\n\n'

		printf 'caller: %s\n' "$(yaml_str "$CALL_CALLER")"

		printf '\nref:\n'
		printf '  genome: %s\n' "$(yaml_str "$CALL_GENOME")"
		printf '  build: %s\n' "$(yaml_str "$CALL_BUILD")"

		printf '\ngatk_resources:\n'
		printf '  panel_of_normals: %s\n' "$(yaml_str "$CALL_PON")"
		printf '  af_only_gnomad: %s\n' "$(yaml_str "$CALL_GNOMAD")"
		printf '  common_biallelic_gnomad: %s\n' "$(yaml_str "$CALL_GNOMAD_CB")"

		printf '\npaths:\n'
		printf '  samples: %s\n' "$(yaml_str "$CALL_SAMPLES_TSV")"
		printf '  bam_folder: %s\n' "$(yaml_str "$CALL_BAM_FOLDER")"
		printf '  output_folder: %s\n' "$(yaml_str "$CALL_OUTPUT_FOLDER")"
		printf '  log_subdir: %s\n' "$(yaml_str "$CALL_LOG_SUBDIR")"

		printf '\nbam:\n'
		printf '  file_extension: %s\n' "$(yaml_str "$CALL_BAM_EXT")"

		printf '\nscatter:\n'
		printf '  mode: %s\n' "$(yaml_str "$CALL_SCATTER")"
		printf '  count: %s\n' "$CALL_SCATTER_COUNT"

		printf '\nparams:\n'
		printf '  mutect2:\n'
		printf '    genotype_germline_sites: %s\n' "$(yaml_bool "$CALL_M2_GERMLINE")"
		printf '    genotype_pon_sites: %s\n' "$(yaml_bool "$CALL_M2_PON_SITES")"
		printf '    disable_adaptive_pruning: %s\n' "$(yaml_bool "$CALL_M2_NO_ADAPTIVE_PRUNING")"
		printf '    extra: %s\n' "$(yaml_str "$CALL_M2_EXTRA")"
		printf '  freebayes:\n'
		printf '    extra: %s\n' "$(yaml_str "$CALL_FB_EXTRA")"
		printf '  bcftools_norm:\n'
		printf '    extra: "-m-any --force -a --atom-overlaps ."\n'
		printf '  bcftools_stats:\n'
		printf '    extra: ""\n'

		printf '\npurecn:\n'
		printf '  enabled: %s\n' "$(yaml_bool "$CALL_PURECN")"
		printf '  genome: %s\n' "$(yaml_str "$CALL_PURECN_GENOME")"
		printf '  intervals_bed: %s\n' "$(yaml_str "$CALL_PURECN_BED")"
		printf '  normaldb: %s\n' "$(yaml_str "$CALL_PURECN_NORMALDB")"
		printf '  mapping_bias: %s\n' "$(yaml_str "$CALL_PURECN_MAPPING_BIAS")"
		printf '  snp_blacklist: %s\n' "$(yaml_str "$CALL_PURECN_SNP_BLACKLIST")"
		printf '  extra: %s\n' "$(yaml_str "$CALL_PURECN_EXTRA")"
		printf '  seed: %s\n' "$CALL_PURECN_SEED"
		printf '  postoptimize: %s\n' "$(yaml_bool "$CALL_PURECN_POSTOPTIMIZE")"
	} >"$outfile"

	success "sm-calling-config.yaml     → ${outfile}"
}

write_calling_samples() {
	local outfile="${OUTPUT_DIR}/sm-calling-samples.tsv"

	{
		printf 'sample\ttumor_bam\tnormal_bam\tanalysis_type\n'
		printf '# sample:        unique identifier — used in output filenames\n'
		printf '# tumor_bam:     BAM basename without extension (not full path)\n'
		printf '# normal_bam:    matched normal BAM basename, or "." if none\n'
		printf '# analysis_type: tumor_only | tumor_normal | germline\n'
		printf '#\n'
		printf '# Examples:\n'
		printf '# IND001_To\tIND001.tumor\t.\ttumor_only\n'
		printf '# IND002_TN\tIND002.tumor\tIND002.normal\ttumor_normal\n'
		printf '# IND003_G\tIND003\t.\tgermline\n'
	} >"$outfile"

	success "sm-calling-samples.tsv     → ${outfile}"
}

# ===========================================================================
# hardnormly workflow
# ===========================================================================

configure_hardnormly() {
	header "hardnormly  ·  VCF → Normalized + Hard-filtered VCF"

	# --- Reference --------------------------------------------------------
	section "Reference"
	prompt HN_FASTA "Reference FASTA" "" \
		"Same reference used for alignment and calling"
	validate_path "ref.fasta" "$HN_FASTA"

	prompt HN_GENOME_FILE "Chromosome sizes (genome) file" "defaults/hg19.genome" \
		"Used by bedtools slop to clamp regions at chromosome boundaries"
	validate_path "ref.genome_file" "$HN_GENOME_FILE"

	prompt_choice HN_BUILD "Genome build" \
		"GRCh37:GRCh37 / hg19" \
		"GRCh38:GRCh38 / hg38"

	# --- Paths ------------------------------------------------------------
	section "Paths"
	prompt HN_VCF_LIST "Input VCF list file" "input/vcfs.txt" \
		"One VCF path per line; run: ls results/calling/**/*.vcf.gz > input/vcfs.txt"
	prompt HN_OUTPUT_FOLDER "Output directory" "results/hardnormly"
	prompt HN_LOG_SUBDIR "Log subdirectory name" "logs"

	# --- Region BED files -------------------------------------------------
	section "Region BED files"
	prompt_array HN_INCLUDE_BEDS \
		"Include region BED files  (target capture / exome regions)" \
		"" \
		"Variants outside these regions get EXCLUDE_REGION annotation; leave empty to process all variants"
	for bed in "${HN_INCLUDE_BEDS[@]}"; do validate_path "include_beds" "$bed"; done

	prompt_array HN_EXCLUDE_BEDS \
		"Exclude region BED files  (blacklist / low-complexity / segmental duplications)" \
		"" \
		"Run 'scripts/generate_exclusion_bed.sh' or 'hardnormly.sh generate-exclusion-bed' to create these"
	for bed in "${HN_EXCLUDE_BEDS[@]}"; do validate_path "exclude_beds" "$bed"; done

	prompt HN_SLOP "Region padding in base pairs (slop)" "100" \
		"Extends include regions by this many bp on each side before annotation"

	# --- Filtering --------------------------------------------------------
	section "Hard filtering"
	prompt_choice HN_FILTER_MODE "Filter source" \
		"caller:Built-in caller preset  (gatk or freebayes)" \
		"file:Custom filters file  (3-column TSV: name action expression)" \
		"none:No hard filters  (normalise only)"

	case "$HN_FILTER_MODE" in
		caller)
			prompt_choice HN_CALLER "Caller preset" \
				"gatk:GATK HaplotypeCaller  (7 filters: DP, VAF, SNP hard, INDEL hard)" \
				"freebayes:FreeBayes  (9 filters: DP, VAF, QUAL, strand bias, read-pos bias)"
			HN_FILTERS_FILE=""
			;;
		file)
			HN_CALLER=""
			prompt HN_FILTERS_FILE "Path to filters file" "" \
				"Format per line: <name> <e|i> <bcftools expression>"
			validate_path "filtering.filters_file" "$HN_FILTERS_FILE"
			;;
		none)
			HN_CALLER=""
			HN_FILTERS_FILE=""
			;;
	esac

	prompt_bool HN_ONLY_PASS "Keep only PASS variants in output?" "n" \
		"Removes soft-filtered variants; set to 'no' to retain all variants with FILTER tags"

	# --- Strip annotations ------------------------------------------------
	section "Annotation stripping  (optional)"
	hint "Removes INFO fields before normalisation and filtering — useful to drop"
	hint "large annotation fields (CSQ, ANN, LOF) that slow down downstream tools"
	prompt_optional HN_STRIP_ANNOTATIONS "INFO fields to strip" \
		"Comma-separated, e.g.: INFO/CSQ,INFO/ANN,INFO/LOF"

	# --- Stats & plots ----------------------------------------------------
	section "Stats & QC"
	prompt_bool HN_GENERATE_STATS "Generate bcftools stats per sample?" "y"
	prompt_bool HN_AUTO_INDEX "Auto-index output VCF with tabix?" "y"
	prompt_bool HN_PLOT_STATS "Generate stats plots?" "n" \
		"Requires plot-vcfstats (matplotlib); only useful when generate_stats=true"

	if [[ "$HN_PLOT_STATS" == "true" ]]; then
		prompt HN_PLOT_OUTPUT_DIR "Stats plots base directory" "${HN_OUTPUT_FOLDER}/plots" \
			"Per-sample plots go in a subdirectory named after each VCF"
	else
		HN_PLOT_OUTPUT_DIR="${HN_OUTPUT_FOLDER}/plots"
	fi
}

write_hardnormly_config() {
	local outfile="${OUTPUT_DIR}/hardnormly-config.yaml"
	mkdir -p "$OUTPUT_DIR"

	{
		printf '# hardnormly config.yaml\n'
		printf '# Generated by hardnormly scripts/generate_config.sh\n'
		printf '# Copy to:  config/config.yaml in this repo\n\n'

		printf 'ref:\n'
		printf '  fasta: %s\n' "$(yaml_str "$HN_FASTA")"
		printf '  genome_file: %s\n' "$(yaml_str "$HN_GENOME_FILE")"
		printf '  build: %s\n' "$(yaml_str "$HN_BUILD")"

		printf '\npaths:\n'
		printf '  vcf_list: %s\n' "$(yaml_str "$HN_VCF_LIST")"
		printf '  output_folder: %s\n' "$(yaml_str "$HN_OUTPUT_FOLDER")"
		printf '  log_subdir: %s\n' "$(yaml_str "$HN_LOG_SUBDIR")"

		printf '\nregions:\n'
		printf '  include_beds:\n'
		if [[ "${#HN_INCLUDE_BEDS[@]}" -gt 0 && -n "${HN_INCLUDE_BEDS[0]}" ]]; then
			for bed in "${HN_INCLUDE_BEDS[@]}"; do printf '    - %s\n' "$(yaml_str "$bed")"; done
		else
			printf '    []  # TODO: add target capture BED file(s)\n'
		fi
		printf '  exclude_beds:\n'
		if [[ "${#HN_EXCLUDE_BEDS[@]}" -gt 0 && -n "${HN_EXCLUDE_BEDS[0]}" ]]; then
			for bed in "${HN_EXCLUDE_BEDS[@]}"; do printf '    - %s\n' "$(yaml_str "$bed")"; done
		else
			printf '    []  # Optional: run generate-exclusion-bed to create blacklist\n'
		fi
		printf '  slop: %s\n' "$HN_SLOP"

		printf '\nfiltering:\n'
		if [[ -n "$HN_CALLER" ]]; then
			printf '  caller: %s\n' "$(yaml_str "$HN_CALLER")"
		else
			printf '  # caller: gatk  # uncomment to use built-in preset\n'
		fi
		if [[ -n "$HN_FILTERS_FILE" ]]; then
			printf '  filters_file: %s\n' "$(yaml_str "$HN_FILTERS_FILE")"
		else
			printf '  # filters_file: ""  # uncomment to use custom filters\n'
		fi
		if [[ -n "$HN_STRIP_ANNOTATIONS" ]]; then
			printf '  strip_annotations: %s\n' "$(yaml_str "$HN_STRIP_ANNOTATIONS")"
		else
			printf '  strip_annotations: ""  # e.g. "INFO/CSQ,INFO/ANN" to drop before filtering\n'
		fi
		printf '  only_pass: %s\n' "$(yaml_bool "$HN_ONLY_PASS")"

		printf '\nprocessing:\n'
		printf '  generate_stats: %s\n' "$(yaml_bool "$HN_GENERATE_STATS")"
		printf '  auto_index: %s\n' "$(yaml_bool "$HN_AUTO_INDEX")"
		printf '  plot_stats: %s\n' "$(yaml_bool "$HN_PLOT_STATS")"
		if [[ "$HN_PLOT_STATS" == "true" ]]; then
			printf '  plot_output_dir: %s\n' "$(yaml_str "$HN_PLOT_OUTPUT_DIR")"
		else
			printf '  # plot_output_dir: ""  # defaults to {output_folder}/plots\n'
		fi
	} >"$outfile"

	success "hardnormly-config.yaml     → ${outfile}"
}

# ===========================================================================
# Confirmation + Summary
# ===========================================================================

confirm_write() {
	[[ "$SKIP_CONFIRM" == "true" || "$NON_INTERACTIVE" == "true" ]] && return

	printf '\n%b─────────────────────────────────────────────────────────────%b\n' "$DIM" "$RESET"
	printf '%bFiles that will be written to:  %s%b\n' "$BOLD" "$OUTPUT_DIR" "$RESET"

	for pipeline in "${SELECTED_PIPELINES[@]}"; do
		case "$pipeline" in
			alignment)
				printf '  sm-alignment-config.yaml\n'
				printf '  sm-alignment-samples.tsv\n'
				;;
			calling)
				printf '  sm-calling-config.yaml\n'
				printf '  sm-calling-samples.tsv\n'
				;;
			hardnormly)
				printf '  hardnormly-config.yaml\n'
				;;
		esac
	done
	printf '%b─────────────────────────────────────────────────────────────%b\n' "$DIM" "$RESET"

	printf '\n%bWrite these files?%b [Y/n]: ' "$BOLD" "$RESET"
	local answer=""
	IFS= read -r answer || true
	answer="${answer,,}"
	case "$answer" in
		n | no)
			info "Aborted — no files written."
			exit 0
			;;
	esac
}

print_summary() {
	printf '\n%b─────────────────────────────────────────────────────────────%b\n' "$DIM" "$RESET"
	printf '%b  Done!%b\n' "$GREEN$BOLD" "$RESET"
	printf '\n'

	for pipeline in "${SELECTED_PIPELINES[@]}"; do
		case "$pipeline" in
			alignment)
				printf '%bsm-alignment next steps:%b\n' "$BOLD" "$RESET"
				printf '  1. Fill in %ssm-alignment-samples.tsv%s with your sample metadata\n' "$CYAN" "$RESET"
				printf '  2. cp %s/sm-alignment-config.yaml  sm-alignment/config/config.yaml\n' "$OUTPUT_DIR"
				printf '  3. cd sm-alignment && snakemake --snakefile workflow/Snakefile --workflow-profile profiles/default\n'
				printf '\n'
				;;
			calling)
				printf '%bsm-vcf-calling next steps:%b\n' "$BOLD" "$RESET"
				printf '  1. Fill in %ssm-calling-samples.tsv%s with your sample metadata\n' "$CYAN" "$RESET"
				printf '  2. cp %s/sm-calling-config.yaml  sm-vcf-calling/config/config.yaml\n' "$OUTPUT_DIR"
				printf '  3. cd sm-vcf-calling && snakemake --snakefile workflow/Snakefile --workflow-profile profiles/default\n'
				printf '\n'
				;;
			hardnormly)
				printf '%bhardnormly next steps:%b\n' "$BOLD" "$RESET"
				printf '  1. Create VCF list:  ls results/calling/**/*.vcf.gz > %s\n' "$HN_VCF_LIST"
				printf '  2. cp %s/hardnormly-config.yaml  config/config.yaml\n' "$OUTPUT_DIR"
				printf '  3. snakemake --snakefile workflow/Snakefile --configfile config/config.yaml\n'
				printf '\n'
				;;
		esac
	done
}

# ===========================================================================
# Main
# ===========================================================================

main() {
	parse_args "$@"
	[[ -n "$OUTPUT_DIR" ]] || die "--output-dir cannot be empty"

	printf '\n%bhardnormly · pipeline config generator%b\n' "$BOLD$BLUE" "$RESET"
	info "  Pipelines : ${SELECTED_PIPELINES[*]}"
	info "  Output dir: ${OUTPUT_DIR}"
	[[ "$NON_INTERACTIVE" == "true" ]] \
		&& warn "  Mode: non-interactive — all prompts skipped, using defaults"

	for pipeline in "${SELECTED_PIPELINES[@]}"; do
		case "$pipeline" in
			alignment) configure_alignment ;;
			calling) configure_calling ;;
			hardnormly) configure_hardnormly ;;
		esac
	done

	confirm_write

	printf '\n'
	for pipeline in "${SELECTED_PIPELINES[@]}"; do
		case "$pipeline" in
			alignment)
				write_alignment_config
				write_alignment_samples
				;;
			calling)
				write_calling_config
				write_calling_samples
				;;
			hardnormly)
				write_hardnormly_config
				;;
		esac
	done

	print_summary
}

main "$@"
