#!/bin/bash
set -e

# 3. append the elementary matrix to the output file of last (see test file)
# 5. add generators for the symmetric group
# 6. add rank n - 1 identity

bold() {
    printf "\033[1m%s\033[0m\n" "$*"
}

if [[ $# -ne 1 ]]; then
  bold "error: expected 1 argument, got $#!"
  exit 1
elif ! [[ "$1" =~ ^[1-9]$ ]] ; then
  bold "error: expected an integer, $1 is not an integer!"
  exit 1
fi

sh/build.sh >> /dev/null

if [[ "$1" =~ ^6$ ]] ; then
  build/libsemigroups/test_bmat8_enum "[004]" | grep "=>" 
fi
