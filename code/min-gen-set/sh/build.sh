#!/bin/bash
set -e

# This script is used to build the things needed to run the run.sh script.

if [ ! -d  build ]; then
  mkdir -p build/output
  cd build
  git clone --branch enum --depth=1 https://github.com/james-d-mitchell/libsemigroups
  cd libsemigroups/extern
  curl -L -O https://github.com/fmtlib/fmt/archive/5.3.0.tar.gz
  tar -xzf 5.3.0.tar.gz && rm -f 5.3.0.tar.gz
  curl -L -O https://github.com/hivert/HPCombi/archive/v0.0.6.tar.gz
  tar -xzf v0.0.6.tar.gz && rm -f v0.0.6.tar.gz
  mv extern/HPCombi-0.0.6 extern/HPCombi

  git clone --depth=1  https://github.com/james-d-mitchell/bliss/tree/05bf42aa825ab04834fc44e0537618e118a06bc0
  mv bliss bliss-0.73

  cp ../../src/Makefile.am .
  cp ../../src/test-bmat8-enum.cpp tests
  ./autogen.sh && ./configure
  cd ../..
fi

cp src/test-bmat8-enum.cpp build/libsemigroups/tests
cd build/libsemigroups && make test_bmat8_enum -j8 && cd ../..
