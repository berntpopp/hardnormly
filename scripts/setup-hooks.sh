#!/bin/bash
# scripts/setup-hooks.sh
# Configures git to use .githooks/ directory for this repository.
# Run once after cloning: bash scripts/setup-hooks.sh

set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd "${script_dir}/.." && pwd)"

git -C "$repo_root" config core.hooksPath .githooks
chmod +x "$repo_root"/.githooks/*
echo "Git hooks configured: .githooks/"
echo "ShellCheck and shfmt will run on staged .sh files before each commit."
