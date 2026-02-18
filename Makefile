# Makefile for hardnormly development tasks
# Usage: make lint, make format, make setup-hooks

SHFMT_FLAGS := -i 0 -bn -ci

# Shell scripts to check (excludes Snakemake launcher which is out of scope)
SH_FILES := hardnormly.sh .githooks/pre-commit scripts/setup-hooks.sh

.PHONY: lint format check-format check-shellcheck setup-hooks help

## Run all lint checks (shellcheck + shfmt)
lint: check-shellcheck check-format

## Run ShellCheck on all shell scripts
check-shellcheck:
	shellcheck $(SH_FILES)

## Check shfmt formatting (no changes, exit 1 if diff)
check-format:
	shfmt -d $(SHFMT_FLAGS) $(SH_FILES)

## Apply shfmt formatting (writes changes in place)
format:
	shfmt -w $(SHFMT_FLAGS) $(SH_FILES)

## Configure git hooks for this repository
setup-hooks:
	bash scripts/setup-hooks.sh

## Show available targets
help:
	@echo "Available targets:"
	@echo "  make lint          - Run ShellCheck + shfmt check"
	@echo "  make format        - Apply shfmt formatting"
	@echo "  make setup-hooks   - Configure git pre-commit hooks"
	@echo "  make help          - Show this help"
