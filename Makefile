# Makefile for hardnormly development tasks
# Usage: make lint, make format, make setup-hooks

SHFMT_FLAGS := -i 0 -bn -ci

# Shell scripts to check (excludes Snakemake launcher which is out of scope)
SH_FILES := hardnormly.sh lib/logging.sh lib/cli.sh lib/genome.sh lib/bed.sh lib/annotate.sh \
	.githooks/pre-commit scripts/setup-hooks.sh tests/setup_bats.sh

.PHONY: lint format check-format check-shellcheck setup-hooks test test-debug help

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

## Run BATS tests
test:
	tests/bats/bin/bats --print-output-on-failure tests/*.bats

## Run BATS tests with verbose output and preserved temp dirs (for debugging)
test-debug:
	tests/bats/bin/bats --no-tempdir-cleanup --verbose-run tests/*.bats

## Show available targets
help:
	@echo "Available targets:"
	@echo "  make lint          - Run ShellCheck + shfmt check"
	@echo "  make format        - Apply shfmt formatting"
	@echo "  make setup-hooks   - Configure git pre-commit hooks"
	@echo "  make test          - Run BATS tests"
	@echo "  make test-debug    - Run BATS tests (verbose, keep temp dirs)"
	@echo "  make help          - Show this help"
