// Copyright (C) 2020 J. D. Mitchell and Finn Smith
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#include <cstddef>        // for size_t
#include <cstdint>        // for uint64_t
#include <iosfwd>         // for ostream, ostringstream, stringbuf
#include <set>            // for set
#include <unordered_set>  // for unordered_set
#include <vector>         // for vector

#include "hpcombi.hpp"         // for BMat8, operator<<
#include "catch.hpp"         // for REQUIRE, REQUIRE_THROWS_AS, REQUIRE_NOTHROW
#include "froidure-pin.hpp"  // for FroidurePin
#include "test-main.hpp"     // for LIBSEMIGROUPS_TEST_CASE
#include "timer.hpp"         // for Timer

namespace {
  bool is_row_reduced(HPCombi::BMat8 const& x) {
    return x.nr_rows() == x.row_space_basis().nr_rows();
  }

  bool is_col_reduced(HPCombi::BMat8 const& x) {
    return is_row_reduced(x.transpose());
  }
}

namespace libsemigroups {
  LIBSEMIGROUPS_TEST_CASE("BMat8", "001", "is_row_reduced/is_col_reduced",
                          "[quick]") {
    HPCombi::BMat8 bm1(0);
    REQUIRE(bm1.transpose() == bm1);
    REQUIRE(::is_row_reduced(bm1));
    REQUIRE(::is_col_reduced(bm1));
    REQUIRE(::is_row_reduced(bm1.one()));
    REQUIRE(::is_col_reduced(bm1.one()));
    REQUIRE(!::is_row_reduced(HPCombi::BMat8(0xFFFFFFFFFF)));
    REQUIRE(!::is_col_reduced(HPCombi::BMat8(0xFFFFFFFFFF)));
  }
}
