#!/bin/bash
# Verify Task 2 artifacts

set -euo pipefail
export PATH=/home/bernt/miniconda3/envs/hardnormly/bin:/usr/bin:/bin

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$REPO_ROOT"

echo "=== VCF list line counts ==="
echo "gatk_vcfs.txt: $(wc -l < tests/data/gatk_vcfs.txt)"
echo "freebayes_vcfs.txt: $(wc -l < tests/data/freebayes_vcfs.txt)"
echo "real_vcfs.txt: $(wc -l < tests/data/real_vcfs.txt)"

echo ""
echo "=== VCF path existence checks ==="
for list in tests/data/gatk_vcfs.txt tests/data/freebayes_vcfs.txt tests/data/real_vcfs.txt; do
    echo "--- $list ---"
    while IFS= read -r f || [[ -n "$f" ]]; do
        [[ -z "$f" ]] && continue
        if test -f "$f"; then
            echo "  OK: $f"
        else
            echo "  MISSING: $f"
        fi
    done < "$list"
done

echo ""
echo "=== YAML config field checks ==="
for config in tests/data/test_config_gatk.yaml tests/data/test_config_freebayes.yaml tests/data/test_config_real.yaml; do
    echo "--- $config ---"
    grep -E 'output_folder|plot_stats|build|filters_file|vcf_list' "$config" || echo "  (grep failed)"
done

echo ""
echo "=== All sections present ==="
for config in tests/data/test_config_gatk.yaml tests/data/test_config_freebayes.yaml tests/data/test_config_real.yaml; do
    echo "--- $config ---"
    for section in ref paths regions filtering processing; do
        grep -q "^${section}:" "$config" && echo "  ${section}: present" || echo "  ${section}: MISSING"
    done
done
