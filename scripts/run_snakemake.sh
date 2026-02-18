#!/bin/bash
#
#SBATCH --job-name=sm_hardnormly
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --time=168:00:00
#SBATCH --mem-per-cpu=2200M
#SBATCH --output=slurm_logs/%x-%j.log

set -euo pipefail
# inherit_errexit requires Bash >= 4.4; guard for portability on older HPC systems
if [[ "${BASH_VERSINFO[0]}" -gt 4 || ("${BASH_VERSINFO[0]}" -eq 4 && "${BASH_VERSINFO[1]}" -ge 4) ]]; then
	shopt -s inherit_errexit
fi

# ── Cluster auto-detection ──────────────────────────────────────────────────

detect_cluster() {
	local hostname
	hostname=$(hostname -f 2>/dev/null || hostname)

	if [[ "$hostname" == *".internal.bih"* ]] || [[ "$hostname" == *"hpc-login"* ]]; then
		echo "bih"
	elif [[ "$hostname" == *"charite"* ]] || [[ "$hostname" == *"hpc2"* ]]; then
		echo "charite"
	else
		echo "local"
	fi
}

CLUSTER=$(detect_cluster)
echo "Detected cluster: $CLUSTER"

# ── Conda activation ───────────────────────────────────────────────────────

if [[ "$CLUSTER" == "charite" ]]; then
	# shellcheck source=/dev/null
	source "$HOME/.bashrc"
fi

eval "$(conda shell.bash hook)"
conda activate snakemake

# ── TMPDIR setup ────────────────────────────────────────────────────────────

if [[ "$CLUSTER" == "bih" ]]; then
	export TMPDIR="${HOME}/scratch/tmp"
elif [[ "$CLUSTER" == "charite" ]]; then
	export TMPDIR="/scratch/${USER}/tmp"
else
	export TMPDIR="${TMPDIR:-/tmp}"
fi

mkdir -p "$TMPDIR"
TMPDIR=$(mktemp -d "${TMPDIR}/hardnormly.XXXXXX")
export TMPDIR
trap 'rm -rf "$TMPDIR"' EXIT

# ── SLURM log directory ────────────────────────────────────────────────────

mkdir -p slurm_logs
export SBATCH_DEFAULTS="--output=slurm_logs/%x-%j.log"

# ── Arguments ───────────────────────────────────────────────────────────────

SNAKEFILE="${1:-workflow/Snakefile}"
CONFIGFILE="${2:-config/config.yaml}"
shift 2 2>/dev/null || true

# ── Profile selection ───────────────────────────────────────────────────────

CLUSTER_PROFILE_ARGS=()
if [[ "$CLUSTER" == "charite" ]]; then
	CLUSTER_PROFILE_ARGS=(--profile profiles/charite)
elif [[ "$CLUSTER" == "bih" ]]; then
	CLUSTER_PROFILE_ARGS=(--profile profiles/bih)
fi

# ── Run Snakemake ───────────────────────────────────────────────────────────

echo "Starting hardnormly Snakemake workflow at $(date)"

snakemake \
	--snakefile "$SNAKEFILE" \
	--configfile "$CONFIGFILE" \
	--workflow-profile profiles/default \
	"${CLUSTER_PROFILE_ARGS[@]}" \
	"$@"

echo "Finished at $(date)"
