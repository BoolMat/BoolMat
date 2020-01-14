#!/bin/bash
set -e

bold() {
    printf "\033[1m%s\033[0m\n" "$*"
}

run() {
  build/libsemigroups/test_bmat8_enum "[$1]" 
  bold "converting C++ output to GAP ..."
  echo "LoadPackage(\"semigroups\", false);; Read(\"src/int_to_bmat.g\");; FileOfIntsToListOfBMats($1); QUIT;" |
    ~/gap/bin/gap.sh -b -A -m 768m -o 1g -T 2>&1
}

if [[ $# -ne 1 ]]; then
  bold "error: expected 1 argument, got $#!"
  exit 1
elif ! [[ "$1" =~ ^[1-9]$ ]] ; then
  bold "error: expected an integer, $1 is not an integer!"
  exit 1
fi



bold "rebuilding ..."
sh/build.sh 

if [[ "$1" =~ ^5$ || "$1" =~ ^6$ || "$1" =~ ^7$ ]] ; then
  run $1
fi
