#!/usr/bin/env bash
# tests/setup_bats.sh — Install BATS test libraries (pinned versions)
# Run once after cloning, or in CI before running tests.
# Clones bats-core, bats-assert, bats-support, bats-file into tests/.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Pinned versions
BATS_CORE_VERSION="v1.11.1"
BATS_ASSERT_VERSION="v2.2.0"
BATS_SUPPORT_VERSION="v0.3.0"
BATS_FILE_VERSION="v0.4.0"

clone_if_missing() {
	local repo="$1"
	local dest="$2"
	local version="$3"

	if [[ -d "$dest" ]]; then
		echo "Already installed: $dest"
		return
	fi

	echo "Installing $repo@$version → $dest"
	git clone -c core.autocrlf=false --depth 1 --branch "$version" "https://github.com/bats-core/$repo.git" "$dest"
}

clone_if_missing bats-core "$SCRIPT_DIR/bats" "$BATS_CORE_VERSION"
clone_if_missing bats-support "$SCRIPT_DIR/test_helper/bats-support" "$BATS_SUPPORT_VERSION"
clone_if_missing bats-assert "$SCRIPT_DIR/test_helper/bats-assert" "$BATS_ASSERT_VERSION"
clone_if_missing bats-file "$SCRIPT_DIR/test_helper/bats-file" "$BATS_FILE_VERSION"

echo "BATS setup complete."
