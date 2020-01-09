#!/bin/bash
set -e

cp test-d-enum.cpp libsemigroups/tests
cd libsemigroups && make test_d_enum -j8 
ln -s libsemigrous/test_d_enum d_enum
