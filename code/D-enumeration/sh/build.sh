#!/bin/bash
set -e

if [ ! -d  build ]; then
  mkdir -p bin
  mkdir -p build && cd build
  git clone --branch v1.0.3 --depth=1 https://github.com/libsemigroups/libsemigroups
  cd libsemigroups/extern
  curl -L -O https://github.com/fmtlib/fmt/archive/5.3.0.tar.gz
  tar -xzf 5.3.0.tar.gz && rm -f 5.3.0.tar.gz
  curl -L -O https://github.com/hivert/HPCombi/archive/v0.0.5.tar.gz
  tar -xzf v0.0.5.tar.gz && rm -f v0.0.5.tar.gz && cd ..
  mv extern/HPCombi-0.0.5 extern/HPCombi
  cp ../../src/Makefile.am .
  ./autogen.sh && ./configure
  cd ../..
fi

cp src/test-d-enum.cpp build/libsemigroups/tests
cd build/libsemigroups && make test_d_enum -j8 && cd ../..
if [ ! -L bin/test_d_enum ]; then
  ln -s build/libsemigroups/test_d_enum bin/test_d_enum
fi
