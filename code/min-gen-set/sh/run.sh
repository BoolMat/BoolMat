#!/bin/bash
set -e

bold() {
    printf "\033[1m%s\033[0m\n" "$*"
}

xbold() {
    printf "\033[1m* %07d %s \033[0m\n" "$1" "$2"
}

GAPROOT="$HOME/gap/"
GAP_SH="$GAPROOT/bin/gap.sh -b -A -m 512m -T -q"

run() {
  bold "* running C++ code ************************************************************"
  build/libsemigroups/test_bmat8_enum "[$1]"
  bold "* using GAP to convert row space numbers to digraphs **************************"
  $GAP_SH << GAP_COMMANDS
  LoadPackage("semigroups", false);;
  Read("src/gap.g");;
  WriteAugmentedDigraphs("build/output/row_space_numbers_$1.txt",
                         "build/output/digraphs_$1.d6.gz",
                         $1);;
  QUIT;
GAP_COMMANDS
  nr=$(cat build/output/digraphs_$1.d6.gz | gzip -d | wc -l)
  xbold $nr "augmented digraphs written to: build/output/digraphs_$1.d6.gz ********"
  bold "* using GAP to find maximal row spaces ****************************************"
  for ((rem=1;i<=$2;i++)); do
      ($GAP_SH << GAP_COMMANDS
      LoadPackage("semigroups");;
      Read("src/gap.g");;
      FilterByHomomorphisms($1, 
                            "build/output/digraphs_$1.d6.gz", 
                            "build/output/max_digraphs_$1_$rem.d6.gz", 
                            7, 
                            $rem);;
      QUIT;
GAP_COMMANDS
) &
  done
  wait
  cat build/output/max_digraphs_$1_*.d6.gz > build/output/max_digraphs_$1.d6.gz
  nr=$(cat build/output/max_digraphs_$1.d6.gz | gzip -d | wc -l)
  xbold $nr "maximal digraphs written to: build/output/max_digraphs_$1.d6.gz ******"
  cat build/output/bmat-gens-$1-*.gz > build/output/bmat-gens-$1.gz
  xbold $((nr+3)) "min. gen. set written to:    build/output/bmat-gens-$1.gz ************"
}

if [[ $# -gt 2 || $# -lt 1 ]]; then
  bold "error: expected 1 or 2 argument, got $#!"
  exit 1
elif ! [[ "$1" =~ ^[1-9]$ ]] ; then
  bold "error: expected an integer as first arg, $1 is not an integer!"
  exit 1
elif [[ $# -eq 2 && "$2" =~ ^[1-9]\+$ ]] ; then
  bold "error: expected an integer as second arg, $2 is not an integer!"
  exit 1
fi
if [[ $# -eq 1 ]] ; then
  nr_cores=1
else
  nr_cores=$2
fi

bold "* rebuilding C++ code *********************************************************"
sh/build.sh

if [[ "$1" =~ ^[3-8]$ ]]; then
  time run $1 $nr_cores
else
  bold "error: expected a value between 4 and 7 (inclusive), got $1"
fi
