#!/usr/bin/env bats
setup() {
	load 'test_helper/common'
	_common_setup
}

@test "debug: exit code after EXIT trap with rm -rf" {
	# Test: does EXIT trap rm -rf affect exit code?
	run bash -c '
set -Eeuo pipefail
_tmp=$(mktemp -d)
trap "rm -rf '"'"'$_tmp'"'"'" EXIT
echo "before exit"
exit 0'
	echo "status=$status output=$output"
	assert_success
}

@test "debug: function that sets trap and exits" {
	run bash -c '
set -Eeuo pipefail
cleanup_handler() {
    local ec=$?
    trap - EXIT
    exit "$ec"
}
trap cleanup_handler EXIT
foo() {
    local _tmp
    _tmp=$(mktemp -d)
    # shellcheck disable=SC2064
    trap "rm -rf '"'"'$_tmp'"'"'" EXIT
    echo hello
}
foo
exit 0'
	echo "status=$status output=$output"
	assert_success
}
