#!/bin/bash
set -e

bold() {
    printf "\033[1m%s\033[0m\n" "$*"
}

xbold() {
    printf "\033[1m* %07d %s \033[0m\n" "$1" "$2"
}

GAP_SH="gap -b -A -m 512m -T -q"

run() {
  if [ $1 -eq 8 ]; then
    build/libsemigroups/test_bmat8_enum "[generate8]"
    build/libsemigroups/test_bmat8_enum "[filter8]"
  else
    build/libsemigroups/test_bmat8_enum "[$1]"
  fi
  bold "* using GAP to convert row space numbers to digraphs **************************"
  $GAP_SH << GAP_COMMANDS
  LoadPackage("semigroups", false);;
  Read("src/gap.g");;
  WriteAugmentedDigraphs("build/output/row_space_numbers_$1_prefiltered.txt",
                         "build/output/digraphs_$1_prefiltered.d6.gz",
                         $1);;
  QUIT;
GAP_COMMANDS
  nr=$(cat build/output/digraphs_$1.d6.gz | gzip -d | wc -l)
  xbold $nr "augmented digraphs written to: build/output/digraphs_$1.d6.gz ********"
  bold "* using GAP to find maximal row spaces ****************************************"
  for ((rem=1;rem<=$2;rem++)); do
      ($GAP_SH << GAP_COMMANDS
      LoadPackage("semigroups");;
      Read("src/gap.g");;
      FilterByHomomorphisms($1,
                            "build/output/digraphs_$1.d6.gz",
                            "build/output/max_digraphs_$1_$rem.d6.gz",
                            $2,
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
  xbold $((nr+4)) "min. gen. set written to:    build/output/bmat-gens-$1.gz ************"
  $GAP_SH << GAP_COMMANDS
    LoadPackage("semigroups");;
    Read("src/gap.g");;
    FileOfBoolMatsToFileOfInts($1,
                               "build/output/bmat-gens-$1.gz",
                               "build/output/bmat-int-gens-$1.txt");;
GAP_COMMANDS
  xbold $((nr+4)) "min. gen. set written to:    build/output/bmat-int-gens-$1.txt *******"

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

sh/build.sh

if [[ "$1" =~ ^[3-8]$ ]]; then
  time run $1 $nr_cores
else
  bold "error: expected a value between 4 and 7 (inclusive), got $1"
fi
