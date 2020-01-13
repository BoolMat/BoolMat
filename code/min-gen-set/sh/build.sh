#!/bin/bash
set -e

# This script is used to build the things needed to run the run.sh script.

if [ ! -d  build ]; then
  mkdir -p build/output build/output
  cd build
  git clone --branch v1.0.3 --depth=1 https://github.com/libsemigroups/libsemigroups
  cd libsemigroups/extern
  curl -L -O https://github.com/fmtlib/fmt/archive/5.3.0.tar.gz
  tar -xzf 5.3.0.tar.gz && rm -f 5.3.0.tar.gz
  curl -L -O https://github.com/hivert/HPCombi/archive/v0.0.5.tar.gz
  tar -xzf v0.0.5.tar.gz && rm -f v0.0.5.tar.gz
  curl -L -O http://www.tcs.hut.fi/Software/bliss/bliss-0.73.zip
  unzip bliss-0.73.zip && rm -f bliss-0.73.zip && cd ..
  sed -i '' -e '62d' extern/bliss-0.73/bliss.cc

  mv extern/HPCombi-0.0.5 extern/HPCombi

  cp ../../src/Makefile.am .
  cp ../../src/test-bmat8-enum.cpp tests
  ./autogen.sh && ./configure
  cd ../..
fi

cp src/test-bmat8-enum.cpp build/libsemigroups/tests
cd build/libsemigroups && make test_bmat8_enum -j8 && cd ../..
