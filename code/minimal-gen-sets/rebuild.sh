#!/bin/bash
set -e

cp test-bmat8-enum.cpp libsemigroups/tests
cd libsemigroups && make test_bmat8_enum -j8
ln -s libsemigroups/test_bmat8_enum test_bmat8_enum
