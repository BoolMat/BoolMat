#!/bin/bash
set -e

cp test-d-enum.cpp libsemigroups/tests
cd libsemigroups && make test_d_enum -j8 
if [ ! -L d_enum ]; then
  ln -s libsemigroups/test_d_enum d_enum
fi
