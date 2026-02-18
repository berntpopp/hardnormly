# Architecture

hardnormly is a modular bash application. The main script (`hardnormly.sh`) orchestrates the pipeline and subcommands by sourcing eight focused library modules.

## Module Structure

| Module | Lines | Purpose | Key Functions |
|--------|-------|---------|---------------|
| `lib/logging.sh` | ~100 | Logging and command execution | `log_msg`, `debug_msg`, `error_msg`, `run_cmd`, `run_cmd_with_retry` |
| `lib/cli.sh` | ~310 | Argument parsing and validation | `parse_args`, `validate_args`, `parse_filter_args`, `show_help` |
| `lib/genome.sh` | ~50 | Genome file creation | `create_genome_file` |
| `lib/bed.sh` | ~95 | BED file processing | `normalize_bed`, `merge_include_beds`, `merge_exclude_beds`, `compress_index_bed`, `create_header_file` |
| `lib/annotate.sh` | ~35 | VCF region annotation | `annotate_vcf_with_regions`, `strip_vcf_annotations` |
| `lib/normalize.sh` | ~45 | VCF normalization | `normalize_vcf` |
| `lib/filter.sh` | ~90 | Filter pipeline construction | `init_filter_pipeline`, `apply_filter_stages`, `write_filtered_output` |
| `lib/stats.sh` | ~45 | Stats generation and plotting | `generate_stats`, `plot_stats_output` |

## Source Order

Modules are sourced in dependency order — `logging.sh` first because other modules depend on `log_msg`, `error_msg`, and `run_cmd`:

```bash
source "${_SCRIPT_DIR}/lib/logging.sh"   # Must be first
source "${_SCRIPT_DIR}/lib/cli.sh"
source "${_SCRIPT_DIR}/lib/genome.sh"
source "${_SCRIPT_DIR}/lib/bed.sh"
source "${_SCRIPT_DIR}/lib/annotate.sh"
source "${_SCRIPT_DIR}/lib/normalize.sh"
source "${_SCRIPT_DIR}/lib/filter.sh"
source "${_SCRIPT_DIR}/lib/stats.sh"
```

Each module uses an include guard to prevent double-sourcing:

```bash
[[ -n "${_LIB_X_LOADED:-}" ]] && return 0
readonly _LIB_X_LOADED=1
```

## Error Handling

The script runs with `set -Eeuo pipefail`:

- `-E` — ERR traps are inherited by functions
- `-e` — Exit on error
- `-u` — Error on undefined variables
- `-o pipefail` — Pipeline fails if any command in the pipe fails

### Patterns

**Single external commands** use `run_cmd` (captures stderr, reports command name on failure):

```bash
run_cmd bgzip -f "$bed_file"
run_cmd tabix -p bed "${bed_file}.gz"
```

**Pipelines** use `|| { error_msg; return 1 }` guards (pipes can't be wrapped with `run_cmd`):

```bash
awk ... "$bed_file" | bedtools sort -i - >"$output_file" \
    || {
        error_msg "normalize_bed: pipeline failed for $bed_file"
        return 1
    }
```

**Cleanup** is handled by a single EXIT trap:

```bash
trap cleanup_handler EXIT
```

The handler preserves the original exit code and removes the temp directory (unless `--no-cleanup` is set). In debug mode, files are listed before deletion.

## Development

### Lint and Format

```bash
make lint           # ShellCheck + shfmt check
make format         # Apply shfmt formatting
make setup-hooks    # Install pre-commit hooks
```

Formatting uses tabs, binary ops at line start, case body indent (`shfmt -i 0 -bn -ci`).

### Testing

```bash
make test           # Run full BATS test suite
make test-debug     # Verbose output, keep temp dirs
```

Test structure:

| File | Tests | What it covers |
|------|-------|----------------|
| `01-smoke.bats` | Smoke | `--help`, `--version`, missing args |
| `02-gatk-filters.bats` | Filter | Each GATK filter tags correctly |
| `03-freebayes-filters.bats` | Filter | Each Freebayes filter tags correctly |
| `04-integration.bats` | Integration | Full pipeline on real data |
| `05-lib-logging.bats` | Unit | `log_msg`, `debug_msg`, `error_msg`, `run_cmd` |
| `06-lib-cli.bats` | Unit | `parse_args`, `validate_args`, `parse_filter_args` |
| `07-lib-modules.bats` | Unit | All other lib/ module functions |

### CI

GitHub Actions runs on every push and PR:

1. ShellCheck on all shell scripts
2. shfmt format check
3. BATS test suite

### Adding a New Module

1. Create `lib/newmodule.sh` with include guard
2. Source it in `hardnormly.sh` (after logging.sh)
3. Add tests in `tests/` (unit tests for functions, integration for pipeline changes)
4. The Makefile auto-discovers `lib/*.sh` — no manual edits needed
