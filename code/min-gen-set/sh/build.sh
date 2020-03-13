#!/bin/bash
set -e

# This script is used to build the things needed to run the run.sh script.

bold() {
    printf "\033[1m%s\033[0m\n" "$*"
}

if [ ! -d  build ]; then
  bold "* creating build directory ****************************************************"
  mkdir -p build/output
  cd build

  bold "* cloning libsemigroups *******************************************************"
  git clone --branch enum --depth=1 https://github.com/james-d-mitchell/libsemigroups
  cd libsemigroups/extern

  bold "* downloading fmt *************************************************************"
  curl -L -O https://github.com/fmtlib/fmt/archive/5.3.0.tar.gz
  tar -xzf 5.3.0.tar.gz && rm -f 5.3.0.tar.gz

  # bold "* downloading HPCombi *********************************************************"
  # curl -L -O https://github.com/hivert/HPCombi/archive/v0.0.6.tar.gz
  # tar -xzf v0.0.6.tar.gz && rm -f v0.0.6.tar.gz
  # mv HPCombi-0.0.6 HPCombi

  bold "* cloning HPCombi *************************************************************"
  git clone --branch bmat8 https://github.com/james-d-mitchell/HPCombi.git
  cd HPCombi && git checkout 4ba8fef && cd ..

  bold "* cloning bliss ***************************************************************"
  git clone --branch master https://github.com/james-d-mitchell/bliss
  cd bliss && git checkout df8f056 && cd ..
  mv bliss bliss-0.73 && cd ..

  bold "* copying files from src/ to build/ *******************************************"
  cp ../../src/Makefile.am .
  cp ../../src/test-bmat8-enum.cpp tests

  bold "* running autoconf on libsemigroups *******************************************"
  ./autogen.sh && ./configure CXXFLAGS="-fopenmp"

  bold "* configuring libsemigroups ***************************************************"
  cd ../..
fi

bold "* making C++ code *************************************************************"
rsync --progress -r -u src/test-bmat8-enum.cpp build/libsemigroups/tests
cd build/libsemigroups && make test_bmat8_enum -j8 && cd ../..
