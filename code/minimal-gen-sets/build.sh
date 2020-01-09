#!/bin/bash
set -e

git clone --branch v1.0.3 --depth=1 https://github.com/libsemigroups/libsemigroups
cd libsemigroups/extern
curl -L -O https://github.com/fmtlib/fmt/archive/5.3.0.tar.gz
tar -xzf 5.3.0.tar.gz && rm -f 5.3.0.tar.gz
curl -L -O https://github.com/hivert/HPCombi/archive/v0.0.5.tar.gz
tar -xzf v0.0.5.tar.gz && rm -f v0.0.5.tar.gz
curl -L -O http://www.tcs.hut.fi/Software/bliss/bliss-0.73.zip
unzip bliss-0.73.zip && rm -f bliss-0.73.zip && cd ..

mv extern/HPCombi-0.0.5 extern/HPCombi

cp ../Makefile.am .
cp ../test-bmat8-enum.cpp tests
./autogen.sh && ./configure
cd ..
./rebuild.sh
