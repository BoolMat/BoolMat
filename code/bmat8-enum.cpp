#include <bitset>
#include <fstream>
#include <iostream>
#include <unordered_set>  // for unordered_set
#include <vector>         // for vector

#include "catch.hpp"
#include "test-main.hpp"

#include "bmat8.hpp"
#include "element.hpp"
#include "froidure-pin.hpp"
#include "runner.hpp"
#include "stl.hpp"
#include "timer.hpp"

#include "../extern/bliss-0.73/graph.hh"

namespace bmat8_enum {
  using namespace libsemigroups;
  typedef bliss_digraphs::Digraph bliss_digraph;

  bliss_digraphs::Stats stats;

  bool is_regular_element(const BMat8 bm) {
    size_t data = bm.to_int();
    return bm
               * BMat8(~(bm * BMat8(~data).transpose() * bm).to_int())
                     .transpose()
               * bm
           == bm;
  }

  std::vector<BMat8> const BMAT8_ONES       = {BMat8(0x0000000000000000),
                                         BMat8(0x8000000000000000),
                                         BMat8(0x8040000000000000),
                                         BMat8(0x8040200000000000),
                                         BMat8(0x8040201000000000),
                                         BMat8(0x8040201008000000),
                                         BMat8(0x8040201008040000),
                                         BMat8(0x8040201008040200),
                                         BMat8(0x8040201008040201)};
  std::vector<BMat8> const BMAT8_ELEMENTARY = {BMat8(0x0000000000000000),
                                               BMat8(0x0000000000000001),
                                               BMat8(0xc040000000000000),
                                               BMat8(0xc040200000000000),
                                               BMat8(0xc040201000000000),
                                               BMat8(0xc040201008000000),
                                               BMat8(0xc040201008040000),
                                               BMat8(0xc040201008040200),
                                               BMat8(0xc040201008040201)};

  BMat8 BMat8_from_rows(std::vector<uint8_t> const& rows) {
    LIBSEMIGROUPS_ASSERT(rows.size() <= 8);
    LIBSEMIGROUPS_ASSERT(0 < rows.size());
    size_t out = 0;
    for (size_t i = 0; i < rows.size(); ++i) {
      out = (out << 8) | rows[i];
    }
    out = out << 8 * (8 - rows.size());
    return BMat8(out);
  }

  template <size_t N>
  BMat8 BMat8_from_perm(typename Perm<N>::type x) {
    LIBSEMIGROUPS_ASSERT(N <= 8);
    std::vector<uint8_t> rows;
    for (size_t i = 0; i < N; ++i) {
      rows.push_back(1 << (7 - x[i]));
    }
    return BMat8_from_rows(rows);
  }

  bliss_digraph bliss_digraph_from_BMat8(BMat8 bm, size_t dim = 8) {
    bliss_digraph out = bliss_digraph(2 * dim);
    size_t        x   = bm.to_int();
    for (size_t i = 0; i < dim; ++i) {
      out.change_color(i, 0);
      out.change_color(dim + i, 1);
      for (size_t j = 0; j < dim; ++j) {
        if ((x >> (63 - 8 * i - j)) & 1) {
          out.add_edge(i, dim + j);
        }
      }
    }
    return out;
  }

  void bliss_hook_function(void*, const unsigned int, const unsigned int*) {}

  BMat8 permuted_BMat8(BMat8 bm, size_t dim, const unsigned int* perm) {
    HPCombi::epu8 row_perm;
    HPCombi::epu8 col_perm;
    for (size_t i = 0; i < dim; ++i) {
      row_perm[perm[i]] = i;
      col_perm[perm[i+dim]-dim] = i;
    }
    for (size_t i = dim; i < 8; ++i) {
      row_perm[i] = i;
      col_perm[i] = i;
    }

    return bm.row_permuted(row_perm).col_permuted(col_perm);
  }

  BMat8 canonical_BMat8(BMat8 bm, size_t dim) {
    bliss_digraph dg = bliss_digraph_from_BMat8(bm, dim);
    return permuted_BMat8(
        bm, dim, dg.canonical_form(stats, &bliss_hook_function, nullptr));
  }

  bool is_row_reduced(BMat8 bm) {
    return bm.nr_rows() == bm.row_space_basis().nr_rows();
  }

  bool is_col_reduced(BMat8 bm) {
    return is_row_reduced(bm.transpose());
  }

  std::bitset<256> row_space_bitset(BMat8 bm) {
    std::bitset<256>     lookup;
    std::vector<uint8_t> row_vec = bm.row_space_basis().rows();
    auto                 last = std::remove(row_vec.begin(), row_vec.end(), 0);
    row_vec.erase(last, row_vec.end());
    for (uint8_t x : row_vec) {
      lookup.set(x);
    }
    lookup.set(0);
    std::vector<uint8_t> row_space(row_vec.begin(), row_vec.end());
    for (size_t i = 0; i < row_space.size(); ++i) {
      for (uint8_t row : row_vec) {
        uint8_t x = row_space[i] | row;
        if (!lookup[x]) {
          row_space.push_back(x);
          lookup.set(x);
        }
      }
    }
    return lookup;
  }

  std::bitset<256> col_space_bitset(BMat8 bm) {
    return row_space_bitset(bm.transpose());
  }
  
  std::vector<uint8_t> row_space_vector(BMat8 bm) {
    std::bitset<256>     lookup;
    std::vector<uint8_t> row_vec = bm.row_space_basis().rows();
    auto                 last = std::remove(row_vec.begin(), row_vec.end(), 0);
    row_vec.erase(last, row_vec.end());
    for (uint8_t x : row_vec) {
      lookup.set(x);
    }
    lookup.set(0);
    std::vector<uint8_t> row_space(row_vec.begin(), row_vec.end());
    for (size_t i = 0; i < row_space.size(); ++i) {
      for (uint8_t row : row_vec) {
        uint8_t x = row_space[i] | row;
        if (!lookup[x]) {
          row_space.push_back(x);
          lookup.set(x);
        }
      }
    }
    return row_space;
  }

  bool is_row_trim(BMat8 bm, size_t dim = 8) {
    std::vector<uint8_t> rows = bm.rows();
    for (size_t i = 0; i < dim; ++i) {
      for (size_t j = 0; j < dim; ++j) {
        if (rows[i] && (i != j) && ((rows[i] | rows[j]) == rows[j])) {
          return false;
        }
      }
    }
    return true;
  }

  bool is_col_trim(BMat8 bm, size_t dim = 8) {
    return is_row_trim(bm.transpose(), dim);
  }

  bool is_trim(BMat8 bm, size_t dim = 8) {
    return is_row_trim(bm, dim) && is_col_trim(bm, dim);
  }

  int nr_ones(uint8_t x) {
    return __builtin_popcount(x);
  }

  // given a BMat8 bm of dimension dim, find all the BMat8s of
  // dimension one higher which are extensions of bm with 0s in the new column
  // and maximal row spaces
  // TODO: obvious optimisations
  std::vector<BMat8> simple_prime_extensions(BMat8 bm, size_t dim) {
    std::vector<uint8_t> rows       = row_space_vector(bm);
    std::vector<BMat8>   extensions(0);

    for (uint8_t new_row = 0; new_row < (1 << dim); ++new_row) {
      bool contains = false;
      for (uint8_t row : rows) {
        if ((row | (new_row << 1)) == (new_row << 1)) {
          contains = true;
          break;
        }
      }
      if (!contains) {
        extensions.push_back(
            BMat8(bm.to_int() | (((new_row << 1) | 1) << (7 - dim))));
      }
    }
    return extensions;
  }

  bool is_fully_indecomposable(BMat8 bm, size_t dim) {
    std::vector<uint8_t> rows = bm.rows();
    for (size_t x = 1; x < (1 << dim - 1); ++x) {
      uint8_t res = 0;
      size_t count = nr_ones(x); 
      for (size_t i = 0; i < dim; ++i) {
        res = ((x >> i) & 1) ? res | rows[i] : res;
      }
      if (nr_ones(res) < count + 1) {
        return false;
      }
    }
    return true;
  }

  class BMatEnumerator : public internal::Runner {
   public:
    BMatEnumerator(size_t dim, bool trim)
        : _n(dim),
          _max((1 << _n) - 1),
          _rows(8, 0),
          _set(),
          _out(),
          _row_seen(),
          _row_orb_by_row(8, std::vector<uint8_t>()),
          _first_row(),
          _min_ones(),
          _trim(trim) {}

   private:
    void dive(size_t k) {
      LIBSEMIGROUPS_ASSERT(k > 0);

      // by permuting columns we can assume if a new 1 appears
      // it appears in the next column
      size_t next_one = 1;
      while ((_rows[k - 1] >> next_one) != 0) {
        next_one += 1;
      }
      next_one = 1 << next_one;
      if (k < _n - 1) {
        for (uint8_t row = _rows[k - 1]; row < _max; ++row) {
          if (!_row_seen[row] && nr_ones(row) >= _min_ones) {
            if (_trim) {
              if ((row > next_one) && ((row & next_one) == 0)) {
                continue;
              }
              for (size_t i = _first_row; i < k; ++i) {
                if (_rows[i] && ((_rows[i] | row) == row)) {
                  goto next_loop;
                }
              }
            }
            _rows[k] = row;

            // UPDATE ROW ORB
            for (size_t i = _first_row; i < k; ++i) {
              for (uint8_t old_row : _row_orb_by_row[i]) {
                uint8_t new_row = old_row | row;
                if (!_row_seen[new_row]) {
                  _row_orb_by_row[k].push_back(new_row);
                  _row_seen[new_row] = true;
                }
              }
            }

            // DIVE
            dive(k + 1);

            // RESET ROW ORB
            for (uint8_t row : _row_orb_by_row[k]) {
              _row_seen[row] = false;
            }
            _row_orb_by_row[k].clear();
          }
        next_loop:;
        }
      } else if (k == _n - 1) {
        for (uint8_t row = _rows[k - 1]; row <= (next_one << 1); ++row) {
          if (!_row_seen[row] && nr_ones(row) >= _min_ones) {
            if (_trim) {
              for (size_t i = _first_row; i < k; ++i) {
                if (_rows[i] && ((_rows[i] | row) == row)) {
                  goto next_loop_last_row;
                }
              }
            }
            _rows[k] = row;
            BMat8 bm = BMat8_from_rows(_rows);
            // move the matrix to the right place
            bm = BMat8(bm.to_int() << (8 - _n));
            if ((_trim && is_col_trim(bm, _n))
                || (!_trim && is_col_reduced(bm))) {
              BMat8 canon = canonical_BMat8(bm, _n);
              if (_set.find(canon) == _set.end()) {
                _set.insert(canon);
                _out.push_back(canon);
              }
            }
          }
        next_loop_last_row:;
          // uint8_t overflow!!!
          if (row == _max) {
            break;
          }
        }
        if (report()) {
          REPORT("found ",
                 _out.size(),
                 " reps so far, currently on \n",
                 internal::to_string(BMat8_from_rows(_rows)));
        }
        _rows[k] = 0;
      }
    }

   public:
    void run() {
      for (bool& x : _row_seen) {
        x = false;
      }
      _row_seen[0] = true;

      set_started(true);

      for (size_t i = 0; i < _n - 1; ++i) {
        _first_row = i;
        for (size_t j = 1; j < _n; ++j) {
          uint8_t row = (1 << j) - 1;
          _rows[i]    = row;
          _min_ones   = nr_ones(row);
          _row_orb_by_row[i].clear();
          _row_orb_by_row[i].push_back(0);
          _row_orb_by_row[i].push_back(row);
          _row_seen[row] = true;

          dive(i + 1);

          _row_seen[row] = false;
        }
        _rows[i] = 0;
      }
      _out.push_back(canonical_BMat8(BMat8(static_cast<size_t>(1) << 63), _n));
      _out.push_back(BMat8(0));
      set_finished(true);
      report_why_we_stopped();
    }

    const std::vector<BMat8>& reps() {
      if (!started()) {
        run();
      }
      return _out;
    }

   private:
    size_t                            _n;
    size_t                            _max;
    std::vector<uint8_t>              _rows;
    std::unordered_set<BMat8>         _set;
    std::vector<BMat8>                _out;
    std::array<bool, 256>             _row_seen;
    std::vector<std::vector<uint8_t>> _row_orb_by_row;
    size_t                            _first_row;
    int                               _min_ones;
    bool                              _trim;
  };

  template <size_t _dim>
  class Filterer : public internal::Runner {
   public:
    Filterer(std::string        in,
             std::string        out,
             std::vector<BMat8> filterers,
             bool               permute)
        : _in(in),
          _out(out),
          _filtered(0),
          _filterers(std::move(filterers)),
          _permute(permute) {}

    void run() {
      set_started(true);

      std::vector<BMat8> bmat_enum;
      std::ifstream      f(_in);
      std::string        line;
      while (std::getline(f, line)) {
        bmat_enum.push_back(BMat8(std::stoul(line)));
      }
      f.close();

      using Perm = typename Perm<_dim>::type;
      HPCombi::epu8 cycle;
      HPCombi::epu8 transposition;
      // TODO: fix for tiny perms
      transposition[0] = 1;
      transposition[1] = 0;
      cycle[0] = 1;
      for (uint8_t i = 2; i < _dim; ++i) {
        cycle[i- 1] = i;
        transposition[i] = i;
      }
      cycle[_dim - 1] = 0;

      std::vector<Perm> S_gens = {Perm(cycle), Perm(transposition)};

      FroidurePin<Perm> S(S_gens);
      S.enumerate();

      std::vector<BMat8> S_bmats;
      for (Perm p : S) {
        S_bmats.push_back(BMat8_from_perm<_dim>(p));
      }

      //TODO evil copying
      if (_filterers.size() == 0) {
        _filterers = bmat_enum;
      }

      std::vector<std::vector<std::bitset<256>>> row_spaces(
          (1 << _dim) + 1, std::vector<std::bitset<256>>());
      for (BMat8 x : _filterers) {
        std::bitset<256> bitset = x.row_space_bitset_ref();
        size_t count = bitset.count();
        row_spaces[count].push_back(bitset);
        if (_permute) {
          for (BMat8 y : S_bmats) {
            bitset = row_space_bitset(x * y);
            row_spaces[count].push_back(bitset);
          }
        }
      }

      for (size_t i = 0; i < bmat_enum.size(); ++i) {
        BMat8            bm     = bmat_enum[i];
        std::bitset<256> bitset = row_space_bitset(bm);
        bool             found  = false;
        if (bitset.count() == 1 << _dim) {
          continue;
        }
        for (size_t i = bitset.count() + 1; i < (1 << _dim); ++i) {
          for (std::bitset<256> bs : row_spaces[i]) {
            if ((bitset | bs) == bs) {
              // TODO: not this
              found = true;
              break;
            }
          }
        }
        if (!found) {
          _filtered.push_back(bm);
        }
        if (report()) {
          REPORT("On ",
                 i + 1,
                 " out of ",
                 bmat_enum.size(),
                 ", keeping ",
                 _filtered.size());
        }
      }

      /*
      std::vector<std::vector<bool>> lookup((1 << _dim) + 1,
                                            std::vector<bool>(0));

      std::vector<std::vector<BMat8>> bmats_by_card((1 << _dim) + 1,
                                                    std::vector<BMat8>(0));
      for (BMat8 bm : bmat_enum) {
        if (bm.row_space_included(BMAT8_ELEMENTARY[_dim])) {
          continue;
        }
        size_t card = bm.row_space_size();
        if (card == (1 << _dim)) {
          continue;
        }
        bmats_by_card[card].push_back(bm);
        lookup[card].push_back(true);
      }

      for (size_t i = 1; i < (1 << _dim) + 1; ++i) {
        size_t m = (1 << _dim) - i;
        for (size_t j = 0; j < bmats_by_card[m].size(); ++j) {
          BMat8 x = bmats_by_card[m][j];
          for (size_t k = m + 1; k < (1 << _dim); ++k) {
            for (size_t l = 0; l < bmats_by_card[k].size(); ++l) {
              if (!lookup[k][l]) {
                continue;
              }
              BMat8 y = bmats_by_card[k][l];
              size_t A = 0;
              for (; A < S_bmats.size() - 1; ++A) {
                std::pair<bool, bool> incl = BMat8::row_space_included2(
                    x, y * S_bmats[A], x, y * (S_bmats[++A]));
                if (incl.first || incl.second) {
                  lookup[m][j] = false;
                  goto next_x;
                }
              }
              if (A == S_bmats.size() - 1) {
                if (x.row_space_included(y * S_bmats[A])) {
                  lookup[m][j] = false;
                  goto next_x;
                }
              }
            }
          }
          if (report()) {
            REPORT("i = ", i);
          }
          _filtered.push_back(x);
          next_x:;
        }
      }
      */

      /*
      std::vector<bool> lookup(bmat_enum.size(), true);

      for (size_t i = 0; i < bmat_enum.size(); ++i) {
        if (bmat_enum[i].row_space_size() == (1 << _dim)) {
          lookup[i] = false;
        }
      }

      for (size_t i = 0; i < bmat_enum.size(); ++i) {
        BMat8 x = bmat_enum[i];
        bool found = false;
        if (x.row_space_size() == (1 << _dim)) {
          goto next_x;
        }
        for (size_t j = 0; j < bmat_enum.size(); ++j) {
          if (i != j && lookup[j]) {
            BMat8 y = bmat_enum[j];
            for (BMat8 const &a : S_bmats) {
              if (x.row_space_included(y * a)) {
                found = true;
                lookup[i] = false;
                goto next_x;
              }
            }
          }
        }
        for (BMat8 const& a : S_bmats) {
          if (x.row_space_included(BMAT8_ELEMENTARY[_dim] * a)) {
            found = true;
            break;
          }
        }
        if (!found) {
          _filtered.push_back(x);
        }
      next_x:
        if (report()) {
          REPORT("On ",
                 i + 1,
                 " out of ",
                 bmat_enum.size(),
                 ", keeping ",
                 _filtered.size());
        }
      }
    */

      std::ofstream o;
      o.open(_out, std::ios::out | std::ios::trunc);
      for (BMat8 i : _filtered) {
        o << i.to_int() << "\n";
      }
      o.close();

      set_finished(true);
      report_why_we_stopped();
    }

    const std::vector<BMat8>& reps() {
      if (!started()) {
        run();
      }
      return _filtered;
    }

   private:
    std::string        _in;
    std::string        _out;
    std::vector<BMat8> _filtered;
    std::vector<BMat8> _filterers;
    bool                _permute;
  };

}  // namespace bmat8_enum

namespace libsemigroups {
  using namespace bmat8_enum;
  LIBSEMIGROUPS_TEST_CASE("BMat8 enum",
                          "001",
                          "test bliss canonicalisation",
                          "[quick][bliss]") {
    BMat8         x   = BMat8({{0, 0, 0}, {0, 0, 1}, {0, 1, 1}});
    BMat8         y   = BMat8({{0, 0, 0}, {0, 0, 1}, {1, 0, 1}});
    bliss_digraph dgx = bliss_digraph_from_BMat8(x, 3);
    bliss_digraph dgy = bliss_digraph_from_BMat8(y, 3);

    dgx.write_dot("dgx.dot");

    bliss_digraph dgx2 = bliss_digraph(6);
    dgx2.change_color(3, 1);
    dgx2.change_color(4, 1);
    dgx2.change_color(5, 1);
    dgx2.add_edge(1, 5);
    dgx2.add_edge(2, 4);
    dgx2.add_edge(2, 5);

    REQUIRE(dgx.cmp(dgx2) == 0);

    BMat8 canonx = permuted_BMat8(
        x, 3, dgx.canonical_form(stats, &bliss_hook_function, nullptr));
    BMat8 canony = permuted_BMat8(
        y, 3, dgy.canonical_form(stats, &bliss_hook_function, nullptr));

    REQUIRE(canonx == canony);
  }

  LIBSEMIGROUPS_TEST_CASE("BMat8 enum",
                          "002",
                          "enumerate B4",
                          "[quick][enumerate]") {
    std::ofstream o;

    auto           rg = ReportGuard();
    BMatEnumerator enumerator(4, false);
    enumerator.report_every(std::chrono::nanoseconds(1));
    REQUIRE(enumerator.reps().size() == 60);
    o.open("bmat_enum_4.txt", std::ios::out | std::ios::trunc);
    for (BMat8 i : enumerator.reps()) {
      o << i.to_int() << "\n";
    }
    o.close();
  }

  LIBSEMIGROUPS_TEST_CASE("BMat8 enum",
                          "003",
                          "enumerate trim B5",
                          "[quick][enumerate]") {
    std::ofstream  o;
    auto           rg = ReportGuard();
    BMatEnumerator enumerator(5, true);
    REQUIRE(enumerator.reps().size() == 873);
    o.open("bmat_enum_5.txt", std::ios::out | std::ios::trunc);
    for (BMat8 i : enumerator.reps()) {
      o << i.to_int() << "\n";
    }
    o.close();
  }

  LIBSEMIGROUPS_TEST_CASE("BMat8 enum",
                          "004",
                          "enumerate trim B6",
                          "[standard][enumerate]") {
    std::ofstream  o;
    auto           rg = ReportGuard();
    BMatEnumerator enumerator(6, true);
    REQUIRE(enumerator.reps().size() == 394);
    o.open("bmat_enum_trim_6.txt", std::ios::out | std::ios::trunc);
    for (BMat8 i : enumerator.reps()) {
      o << i.to_int() << "\n";
    }
    o.close();
  }

  LIBSEMIGROUPS_TEST_CASE("BMat8 enum",
                          "800",
                          "enumerate trim B8",
                          "[extreme][enumerate]") {
    std::ofstream o;
    auto           rg = ReportGuard();
    BMatEnumerator enumerator_8_trim(8, true);
    o.open("bmat_trim_enum_8.txt", std::ios::out | std::ios::trunc);
    for (BMat8 i : enumerator_8_trim.reps()) {
      o << i.to_int() << "\n";
    }
    o.close();
    REQUIRE(enumerator_8_trim.reps().size() == 17120845);
  }

  LIBSEMIGROUPS_TEST_CASE("BMat8 enum", "006", "filter 5", "[quick]") {
    Filterer<5> f("bmat_trim_enum_5.txt", "bmat_gens_test_5.txt", {}, true);
    f.run();
    std::vector<BMat8> filtered = f.reps();
    REQUIRE(filtered.size() == 9);
  }

  LIBSEMIGROUPS_TEST_CASE("BMat8 enum", "007", "filter 6", "[quick]") {
    Filterer<6> f("bmat_enum_trim_6.txt", "bmat_gens_test_6.txt", {}, true);
    f.run();
    std::vector<BMat8> filtered = f.reps();
    REQUIRE(filtered.size() == 65);
  }

  LIBSEMIGROUPS_TEST_CASE("BMat8 enum", "008", "filter 7", "[extreme]") {
    auto        rg = ReportGuard();
    Filterer<7> f("bmat_trim_enum_7.txt", "bmat_gens_test_7.txt", {}, true);
    f.run();
    std::vector<BMat8> filtered = f.reps();
    REQUIRE(filtered.size() == 2143);
  }

  ////////////////////////////////////////////////////////////////////
  //
  // More experimental tests
  //
  // /////////////////////////////////////////////////////////////////
  //
  LIBSEMIGROUPS_TEST_CASE("BMat8 enum",
                          "110",
                          "filter 6 faster?",
                          "[extreme]") {
    std::vector<BMat8> bmat6_enum;
    std::ifstream      f("bmat_trim_enum_6.txt");
    std::string        line;
    while (std::getline(f, line)) {
      bmat6_enum.push_back(BMat8(std::stoul(line)));
    }
    f.close();
    REQUIRE(bmat6_enum.size() == 394);

    std::cout << "finished reading!" << std::endl;

    BMat8 E({{1, 1, 0, 0, 0, 0},
             {0, 1, 0, 0, 0, 0},
             {0, 0, 1, 0, 0, 0},
             {0, 0, 0, 1, 0, 0},
             {0, 0, 0, 0, 1, 0},
             {1, 0, 0, 0, 0, 1}});

    using Perm = typename Perm<6>::type;
    const std::vector<Perm> S6_gens
        = {Perm({1, 2, 3, 4, 5, 0}), Perm({1, 0, 2, 3, 4, 5})};

    FroidurePin<Perm> S6(S6_gens);
    REQUIRE(S6.size() == 720);
    std::cout << "finished computing S6!" << std::endl;

    std::vector<BMat8> S6_bmats;
    for (Perm p : S6) {
      S6_bmats.push_back(BMat8_from_perm<6>(p));
    }

    std::vector<std::vector<std::bitset<256>>> row_spaces(
        65, std::vector<std::bitset<256>>(0));
    size_t             count = 0;
    size_t             kept  = 0;
    std::vector<BMat8> filtered(0);
    for (BMat8 x : bmat6_enum) {
      /*
      if (is_regular_element(x)) {
        continue;
      }
      */
      size_t                        size = row_space_bitset(x).count();
      std::vector<std::bitset<256>> new_bitsets(0);
      bool                          add_new_bitsets = true;
      for (BMat8 y : S6_bmats) {
        /*
        if (row_space_bitset(x * y * x).count() == size) {
          add_new_bitsets = false;
          break;
        }
        */
        std::bitset<256> bitset = row_space_bitset(x * y);
        new_bitsets.push_back(bitset);
      }
      if (add_new_bitsets) {
        for (std::bitset<256> bs : new_bitsets) {
          row_spaces[bs.count()].push_back(bs);
        }
        filtered.push_back(x);
        kept++;
      }
    }

    std::cout << "kept " << kept << " out of 394" << std::endl;

    for (BMat8 y : S6_bmats) {
      std::bitset<256> bitset = row_space_bitset(E * y);
      row_spaces[bitset.count()].push_back(bitset);
    }

    std::cout << "finished producing bitsets!" << std::endl;
    size_t i = 0;

    for (size_t i = 1; i < 65; ++i) {
      std::unordered_set<std::bitset<256>> set(row_spaces[i].begin(),
                                               row_spaces[i].end());
      row_spaces[i].assign(set.begin(), set.end());
    }

    std::cout << "removed duplicates!" << std::endl;

    std::vector<BMat8> gens(0);

    for (BMat8 bm : filtered) {
      std::bitset<256> bitset = row_space_bitset(bm);
      bool             found  = false;
      // TODO: not this
      if (bitset.count() == 64) {
        found = true;  // get rid of any permutation matrices
      }
      for (size_t i = bitset.count() + 1; i < 64 && !found; ++i) {
        for (size_t j = 0; j < row_spaces[i].size(); ++j) {
          if ((bitset | row_spaces[i][j]) == row_spaces[i][j]) {
            found = true;
            break;
          }
        }
      }
      if (!found) {
        gens.push_back(bm);
      }
    }
    for (Perm p : S6_gens) {
      gens.push_back(BMat8_from_perm<6>(p));
    }

    gens.push_back(E);
    gens.push_back(BMAT8_ONES[5]);
    std::ofstream o;
    o.open("bmat_gens_6_wrong?.txt", std::ios::out | std::ios::trunc);
    for (BMat8 i : gens) {
      o << i.to_int() << "\n";
    }
    o.close();
    REQUIRE(gens.size() == 69);
  }

  LIBSEMIGROUPS_TEST_CASE("BMat8 enum",
                          "111",
                          "filter 8 by elementary containment",
                          "[extreme]") {
    std::vector<BMat8> bmat8_enum;
    std::ifstream      f("bmat_trim_enum_8.txt");
    std::string        line;
    while (std::getline(f, line)) {
      bmat8_enum.push_back(BMat8(std::stoul(line)));
    }
    f.close();
    REQUIRE(bmat8_enum.size() == 17120845);

    std::cout << "finished reading!" << std::endl;

    BMat8 E({{1, 1, 0, 0, 0, 0, 0, 0},
             {0, 1, 0, 0, 0, 0, 0, 0},
             {0, 0, 1, 0, 0, 0, 0, 0},
             {0, 0, 0, 1, 0, 0, 0, 0},
             {0, 0, 0, 0, 1, 0, 0, 0},
             {0, 0, 0, 0, 0, 1, 0, 0},
             {0, 0, 0, 0, 0, 0, 1, 0},
             {1, 0, 0, 0, 0, 0, 0, 1}});

    using Perm = typename Perm<8>::type;
    const std::vector<Perm> S8_gens
        = {Perm({1, 2, 3, 4, 5, 6, 7, 0}), Perm({1, 0, 2, 3, 4, 5, 6, 7})};

    FroidurePin<Perm> S8(S8_gens);
    REQUIRE(S8.size() == 40320);
    std::cout << "finished computing S8!" << std::endl;

    std::vector<BMat8> S8_bmats;
    for (Perm p : S8) {
      S8_bmats.push_back(BMat8_from_perm<8>(p));
    }

    std::vector<std::bitset<256>> row_spaces(0);
    for (BMat8 y : S8_bmats) {
      std::bitset<256> bitset = row_space_bitset(E * y);
      row_spaces.push_back(bitset);
    }

    std::cout << "finished producing bitsets!" << std::endl;

    std::vector<BMat8> gens(0);

    size_t             count = 0;
    std::vector<BMat8> filtered(0);
    for (BMat8 bm : bmat8_enum) {
      std::bitset<256> bitset = row_space_bitset(bm);
      bool             found  = false;
      // TODO: not this
      if (bitset.count() == 256) {
        found = true;  // get rid of any permutation matrices
      }
      for (std::bitset<256> bs : row_spaces) {
        if ((bitset | bs) == bs) {
          found = true;
          break;
        }
      }

      if (!found) {
        filtered.push_back(bm);
      }
      if (count++ % 10000 == 0) {
        std::cout << "on " << count << ", keeping " << filtered.size()
                  << std::endl;
      }
    }

    std::ofstream o;
    o.open("bmat_filtered_8.txt", std::ios::out | std::ios::trunc);
    for (BMat8 i : filtered) {
      o << i.to_int() << "\n";
    }
    o.close();
    REQUIRE(filtered.size() == 17117115);
  }

  LIBSEMIGROUPS_TEST_CASE("BMat8 enum",
                          "112",
                          "filter 8 by extended 7 gens",
                          "[extreme]") {
    std::vector<BMat8> bmat7_gens(0);
    std::vector<BMat8> bmat8_enum(0);
    std::ifstream      f("bmat_gens_7.txt");
    std::string        line;
    while (std::getline(f, line)) {
      bmat7_gens.push_back(BMat8(std::stoul(line) | 1));
    }
    f.close();
    f = std::ifstream("bmat_filtered_8.txt");
    while (std::getline(f, line)) {
      bmat8_enum.push_back(BMat8(std::stoul(line)));
    }
    f.close();
    REQUIRE(bmat8_enum.size() == 17117115);

    std::cout << "finished reading!" << std::endl;

    size_t                                     i = 0;
    std::vector<std::vector<std::bitset<256>>> row_spaces(
        257, std::vector<std::bitset<256>>());

    std::vector<std::vector<size_t>> containment_count(257,
                                                       std::vector<size_t>(0));

    for (BMat8 x : bmat7_gens) {
      std::bitset<256> bitset = row_space_bitset(x);
      size_t           count  = bitset.count();
      if (count == 256) {
        continue;
      }
      row_spaces[count].push_back(bitset);
      containment_count[count].push_back(0);
      if (++i % 100 == 0) {
        std::cout << i << std::endl;
      }
    }

    std::unordered_set<BMat8> canon_gens;
    for (BMat8 bm : bmat7_gens) {
      canon_gens.insert(canonical_BMat8(bm, 8));
    }

    std::cout << "finished producing bitsets!" << std::endl;

    std::vector<BMat8> filtered(0);
    size_t             count = 0;
    for (BMat8 bm : bmat8_enum) {
      std::bitset<256> bitset = row_space_bitset(bm);
      bool             found  = false;
      size_t           size   = bitset.count();
      // TODO: not this
      if (size == 256) {
        found = true;  // get rid of any permutation matrices
      }
      if (canon_gens.find(bm) != canon_gens.end()) {
        found = true;
      }
      for (size_t i = size + 1; i < 257 && !found; ++i) {
        for (size_t j = 0; j < row_spaces[i].size(); ++j) {
          std::bitset<256> bs = row_spaces[i][j];
          if ((bitset | bs) == bs) {
            found = true;
            containment_count[i][j]++;
            break;
          }
        }
      }
      if (!found) {
        filtered.push_back(bm);
      }
      if (++count % 100000 == 0) {
        std::cout << "on " << count << ", keeping " << filtered.size()
                  << std::endl;
      }
    }

    for (size_t i = 0; i < 256; ++i) {
      for (size_t j = 0; j < containment_count[i].size(); ++j) {
        std::cout << "row space " << i << "," << j
                  << " killed: " << containment_count[i][j] << std::endl;
      }
    }

    std::ofstream o;
    o.open("bmat_filtered_8_2.txt", std::ios::out | std::ios::trunc);
    for (BMat8 i : filtered) {
      o << i.to_int() << "\n";
    }
    o.close();
    REQUIRE(filtered.size() == 15334021);
  }

  LIBSEMIGROUPS_TEST_CASE("BMat8 enum",
                          "114",
                          "filter 8 by special 7 gens",
                          "[extreme]") {
    std::vector<BMat8> bmat7_gens(0);
    std::vector<BMat8> bmat8_enum(0);
    std::ifstream      f("bmat_gens_7.txt");
    std::string        line;
    while (std::getline(f, line)) {
      bmat7_gens.push_back(BMat8(std::stoul(line) | 1));
    }
    f.close();
    f = std::ifstream("bmat_filtered_8_2.txt");
    while (std::getline(f, line)) {
      bmat8_enum.push_back(BMat8(std::stoul(line)));
    }
    f.close();
    REQUIRE(bmat8_enum.size() == 15334021);

    std::cout << "finished reading!" << std::endl;

    using Perm = typename Perm<8>::type;
    const std::vector<Perm> S8_gens
        = {Perm({1, 2, 3, 4, 5, 6, 7, 0}), Perm({1, 0, 2, 3, 4, 5, 6, 7})};

    FroidurePin<Perm> S8(S8_gens);
    REQUIRE(S8.size() == 40320);
    std::cout << "finished computing S8!" << std::endl;

    std::vector<BMat8> S8_bmats;
    for (Perm p : S8) {
      S8_bmats.push_back(BMat8_from_perm<8>(p));
    }

    size_t                                     i = 0;
    std::vector<std::vector<std::bitset<256>>> row_spaces(
        257, std::vector<std::bitset<256>>());

    for (BMat8 x : bmat7_gens) {
      std::bitset<256> bitset = row_space_bitset(x);
      size_t           count  = bitset.count();
      if (count != 160) {
        continue;
      }
      row_spaces[count].push_back(bitset);
      for (BMat8 const& A : S8_bmats) {
        row_spaces[count].push_back(row_space_bitset(x * A));
      }
      if (++i % 100 == 0) {
        std::cout << i << std::endl;
      }
    }

    std::unordered_set<BMat8> canon_gens;
    for (BMat8 bm : bmat7_gens) {
      canon_gens.insert(canonical_BMat8(bm, 8));
    }

    std::cout << "finished producing bitsets!" << std::endl;

    std::vector<BMat8> filtered(0);
    size_t             count = 0;
    for (BMat8 bm : bmat8_enum) {
      std::bitset<256> bitset = row_space_bitset(bm);
      bool             found  = false;
      size_t           size   = bitset.count();
      // TODO: not this
      if (size == 256) {
        found = true;  // get rid of any permutation matrices
      }
      if (canon_gens.find(bm) != canon_gens.end()) {
        found = true;
      }
      for (size_t i = size + 1; i < 257 && !found; ++i) {
        for (size_t j = 0; j < row_spaces[i].size(); ++j) {
          std::bitset<256> bs = row_spaces[i][j];
          if ((bitset | bs) == bs) {
            found = true;
            break;
          }
        }
      }
      if (!found) {
        filtered.push_back(bm);
      }
      if (++count % 100000 == 0) {
        std::cout << "on " << count << ", keeping " << filtered.size()
                  << std::endl;
      }
    }

    std::ofstream o;
    o.open("bmat_filtered_8_3.txt", std::ios::out | std::ios::trunc);
    for (BMat8 i : filtered) {
      o << i.to_int() << "\n";
    }
    o.close();
    REQUIRE(filtered.size() == 2923382);
  }
  
  LIBSEMIGROUPS_TEST_CASE("BMat8 enum",
                          "115",
                          "filter 8 by special 7 gens 2",
                          "[extreme]") {
    std::vector<BMat8> bmat7_gens(0);
    std::vector<BMat8> bmat8_enum(0);
    std::ifstream      f("bmat_gens_7.txt");
    std::string        line;
    while (std::getline(f, line)) {
      bmat7_gens.push_back(BMat8(std::stoul(line) | 1));
    }
    f.close();
    f = std::ifstream("bmat_filtered_8_3.txt");
    while (std::getline(f, line)) {
      bmat8_enum.push_back(BMat8(std::stoul(line)));
    }
    f.close();
    REQUIRE(bmat8_enum.size() == 2923382);

    std::cout << "finished reading!" << std::endl;

    using Perm = typename Perm<8>::type;
    const std::vector<Perm> S8_gens
        = {Perm({1, 2, 3, 4, 5, 6, 7, 0}), Perm({1, 0, 2, 3, 4, 5, 6, 7})};

    FroidurePin<Perm> S8(S8_gens);
    REQUIRE(S8.size() == 40320);
    std::cout << "finished computing S8!" << std::endl;
    std::vector<BMat8> S8_bmats;
    for (Perm p : S8) {
      S8_bmats.push_back(BMat8_from_perm<8>(p));
    }

    size_t                                     i = 0;
    std::vector<std::vector<std::bitset<256>>> row_spaces(
        257, std::vector<std::bitset<256>>());

    for (BMat8 x : bmat7_gens) {
      std::bitset<256> bitset = row_space_bitset(x);
      size_t           count  = bitset.count();
      if (!(count == 128 || count == 136 || count == 144)) {
        continue;
      }
      row_spaces[count].push_back(bitset);
      for (BMat8 const& A : S8_bmats) {
        row_spaces[count].push_back(row_space_bitset(x * A));
      }
      if (++i % 100 == 0) {
        std::cout << i << std::endl;
      }
    }

    std::unordered_set<BMat8> canon_gens;
    for (BMat8 bm : bmat7_gens) {
      canon_gens.insert(canonical_BMat8(bm, 8));
    }

    std::cout << "finished producing bitsets!" << std::endl;

    std::vector<BMat8> filtered(0);
    size_t             count = 0;
    for (BMat8 bm : bmat8_enum) {
      std::bitset<256> bitset = row_space_bitset(bm);
      bool             found  = false;
      size_t           size   = bitset.count();
      // TODO: not this
      if (size == 256) {
        found = true;  // get rid of any permutation matrices
      }
      if (canon_gens.find(bm) != canon_gens.end()) {
        found = true;
      }
      for (size_t i = size + 1; i < 257 && !found; ++i) {
        for (size_t j = 0; j < row_spaces[i].size(); ++j) {
          std::bitset<256> bs = row_spaces[i][j];
          if ((bitset | bs) == bs) {
            found = true;
            break;
          }
        }
      }
      if (!found) {
        filtered.push_back(bm);
      }
      if (++count % 1000 == 0) {
        std::cout << "on " << count << ", keeping " << filtered.size()
                  << std::endl;
      }
    }

    std::ofstream o;
    o.open("bmat_filtered_8_4.txt", std::ios::out | std::ios::trunc);
    for (BMat8 i : filtered) {
      o << i.to_int() << "\n";
    }
    o.close();
    REQUIRE(filtered.size() == 1661260);
  }
  
  LIBSEMIGROUPS_TEST_CASE("BMat8 enum",
                          "116",
                          "filter 8 by special 7 gens 3",
                          "[extreme]") {
    std::vector<BMat8> bmat7_gens(0);
    std::vector<BMat8> bmat8_enum(0);
    std::ifstream      f("bmat_gens_7.txt");
    std::string        line;
    while (std::getline(f, line)) {
      bmat7_gens.push_back(BMat8(std::stoul(line) | 1));
    }
    f.close();
    f = std::ifstream("bmat_filtered_8_4.txt");
    while (std::getline(f, line)) {
      bmat8_enum.push_back(BMat8(std::stoul(line)));
    }
    f.close();
    REQUIRE(bmat8_enum.size() == 1661260);

    std::cout << "finished reading!" << std::endl;

    using Perm = typename Perm<8>::type;
    const std::vector<Perm> S8_gens
        = {Perm({1, 2, 3, 4, 5, 6, 7, 0}), Perm({1, 0, 2, 3, 4, 5, 6, 7})};

    FroidurePin<Perm> S8(S8_gens);
    REQUIRE(S8.size() == 40320);
    std::cout << "finished computing S8!" << std::endl;

    std::vector<BMat8> S8_bmats;
    for (Perm p : S8) {
      S8_bmats.push_back(BMat8_from_perm<8>(p));
    }

    size_t                                     i = 0;
    std::vector<std::vector<std::bitset<256>>> row_spaces(
        257, std::vector<std::bitset<256>>());

    for (BMat8 x : bmat7_gens) {
      std::bitset<256> bitset = row_space_bitset(x);
      size_t           count  = bitset.count();
      if (!(count == 112 || count == 116 || count == 120)) {
        continue;
      }
      row_spaces[count].push_back(bitset);
      for (BMat8 const& A : S8_bmats) {
        row_spaces[count].push_back(row_space_bitset(x * A));
      }
      if (++i % 10000 == 0) {
        std::cout << i << std::endl;
      }
    }

    std::unordered_set<BMat8> canon_gens;
    for (BMat8 bm : bmat7_gens) {
      canon_gens.insert(canonical_BMat8(bm, 8));
    }

    std::cout << "finished producing bitsets!" << std::endl;

    std::vector<BMat8> filtered(0);
    size_t             count = 0;
    for (BMat8 bm : bmat8_enum) {
      std::bitset<256> bitset = row_space_bitset(bm);
      bool             found  = false;
      size_t           size   = bitset.count();
      // TODO: not this
      if (size == 256) {
        found = true;  // get rid of any permutation matrices
      }
      if (canon_gens.find(bm) != canon_gens.end()) {
        found = true;
      }
      for (size_t i = size + 1; i < 257 && !found; ++i) {
        for (size_t j = 0; j < row_spaces[i].size(); ++j) {
          std::bitset<256> bs = row_spaces[i][j];
          if ((bitset | bs) == bs) {
            found = true;
            break;
          }
        }
      }
      if (!found) {
        filtered.push_back(bm);
      }
      if (++count % 10000 == 0) {
        std::cout << "on " << count << ", keeping " << filtered.size()
                  << std::endl;
      }
    }

    std::ofstream o;
    o.open("bmat_filtered_8_5.txt", std::ios::out | std::ios::trunc);
    for (BMat8 i : filtered) {
      o << i.to_int() << "\n";
    }
    o.close();
    REQUIRE(filtered.size() == 1206208);
  }
  
  LIBSEMIGROUPS_TEST_CASE("BMat8 enum",
                          "117",
                          "filter 8 by special 7 gens 4",
                          "[extreme]") {
    std::vector<BMat8> bmat7_gens(0);
    std::vector<BMat8> bmat8_enum(0);
    std::ifstream      f("bmat_gens_7.txt");
    std::string        line;
    while (std::getline(f, line)) {
      bmat7_gens.push_back(BMat8(std::stoul(line) | 1));
    }
    f.close();
    f = std::ifstream("bmat_filtered_8_5.txt");
    while (std::getline(f, line)) {
      bmat8_enum.push_back(BMat8(std::stoul(line)));
    }
    f.close();
    REQUIRE(bmat8_enum.size() == 1206208);

    std::cout << "finished reading!" << std::endl;

    using Perm = typename Perm<8>::type;
    const std::vector<Perm> S8_gens
        = {Perm({1, 2, 3, 4, 5, 6, 7, 0}), Perm({1, 0, 2, 3, 4, 5, 6, 7})};

    FroidurePin<Perm> S8(S8_gens);
    REQUIRE(S8.size() == 40320);
    std::cout << "finished computing S8!" << std::endl;

    std::vector<BMat8> S8_bmats;
    for (Perm p : S8) {
      S8_bmats.push_back(BMat8_from_perm<8>(p));
    }

    size_t                                     i = 0;
    std::vector<std::vector<std::bitset<256>>> row_spaces(
        257, std::vector<std::bitset<256>>());

    for (BMat8 x : bmat7_gens) {
      std::bitset<256> bitset = row_space_bitset(x);
      size_t           count  = bitset.count();
      if (!(count == 104 || count == 108)) {
        continue;
      }
      row_spaces[count].push_back(bitset);
      for (BMat8 const& A : S8_bmats) {
        row_spaces[count].push_back(row_space_bitset(x * A));
      }
      if (++i % 10000 == 0) {
        std::cout << i << std::endl;
      }
    }

    std::unordered_set<BMat8> canon_gens;
    for (BMat8 bm : bmat7_gens) {
      canon_gens.insert(canonical_BMat8(bm, 8));
    }

    std::cout << "finished producing bitsets!" << std::endl;

    std::vector<BMat8> filtered(0);
    size_t             count = 0;
    for (BMat8 bm : bmat8_enum) {
      std::bitset<256> bitset = row_space_bitset(bm);
      bool             found  = false;
      size_t           size   = bitset.count();
      // TODO: not this
      if (size == 256) {
        found = true;  // get rid of any permutation matrices
      }
      if (canon_gens.find(bm) != canon_gens.end()) {
        found = true;
      }
      for (size_t i = size + 1; i < 257 && !found; ++i) {
        for (size_t j = 0; j < row_spaces[i].size(); ++j) {
          std::bitset<256> bs = row_spaces[i][j];
          if ((bitset | bs) == bs) {
            found = true;
            break;
          }
        }
      }
      if (!found) {
        filtered.push_back(bm);
      }
      if (++count % 10000 == 0) {
        std::cout << "on " << count << ", keeping " << filtered.size()
                  << std::endl;
      }
    }

    std::ofstream o;
    o.open("bmat_filtered_8_6.txt", std::ios::out | std::ios::trunc);
    for (BMat8 i : filtered) {
      o << i.to_int() << "\n";
    }
    o.close();
    REQUIRE(filtered.size() == 974417);
  }
  
  LIBSEMIGROUPS_TEST_CASE("BMat8 enum",
                          "118",
                          "filter 8 by special 7 gens 5",
                          "[extreme]") {
    std::vector<BMat8> bmat7_gens(0);
    std::vector<BMat8> bmat8_enum(0);
    std::ifstream      f("bmat_gens_7.txt");
    std::string        line;
    while (std::getline(f, line)) {
      bmat7_gens.push_back(BMat8(std::stoul(line) | 1));
    }
    f.close();
    f = std::ifstream("bmat_filtered_8_6.txt");
    while (std::getline(f, line)) {
      bmat8_enum.push_back(BMat8(std::stoul(line)));
    }
    f.close();
    REQUIRE(bmat8_enum.size() == 974417);

    std::cout << "finished reading!" << std::endl;

    using Perm = typename Perm<8>::type;
    const std::vector<Perm> S8_gens
        = {Perm({1, 2, 3, 4, 5, 6, 7, 0}), Perm({1, 0, 2, 3, 4, 5, 6, 7})};

    FroidurePin<Perm> S8(S8_gens);
    REQUIRE(S8.size() == 40320);
    std::cout << "finished computing S8!" << std::endl;

    std::vector<BMat8> S8_bmats;
    for (Perm p : S8) {
      S8_bmats.push_back(BMat8_from_perm<8>(p));
    }

    size_t                                     i = 0;
    std::vector<std::vector<std::bitset<256>>> row_spaces(
        257, std::vector<std::bitset<256>>());

    for (BMat8 x : bmat7_gens) {
      std::bitset<256> bitset = row_space_bitset(x);
      size_t           count  = bitset.count();
      if (!(count == 100 || count == 102)) {
        continue;
      }
      row_spaces[count].push_back(bitset);
      for (BMat8 const& A : S8_bmats) {
        row_spaces[count].push_back(row_space_bitset(x * A));
      }
      if (++i % 10000 == 0) {
        std::cout << i << std::endl;
      }
    }

    std::unordered_set<BMat8> canon_gens;
    for (BMat8 bm : bmat7_gens) {
      canon_gens.insert(canonical_BMat8(bm, 8));
    }

    std::cout << "finished producing bitsets!" << std::endl;

    std::vector<BMat8> filtered(0);
    size_t             count = 0;
    for (BMat8 bm : bmat8_enum) {
      std::bitset<256> bitset = row_space_bitset(bm);
      bool             found  = false;
      size_t           size   = bitset.count();
      // TODO: not this
      if (size == 256) {
        found = true;  // get rid of any permutation matrices
      }
      if (canon_gens.find(bm) != canon_gens.end()) {
        found = true;
      }
      for (size_t i = size + 1; i < 257 && !found; ++i) {
        for (size_t j = 0; j < row_spaces[i].size(); ++j) {
          std::bitset<256> bs = row_spaces[i][j];
          if ((bitset | bs) == bs) {
            found = true;
            break;
          }
        }
      }
      if (!found) {
        filtered.push_back(bm);
      }
      if (++count % 10000 == 0) {
        std::cout << "on " << count << ", keeping " << filtered.size()
                  << std::endl;
      }
    }

    std::ofstream o;
    o.open("bmat_filtered_8_7.txt", std::ios::out | std::ios::trunc);
    for (BMat8 i : filtered) {
      o << i.to_int() << "\n";
    }
    o.close();
    REQUIRE(filtered.size() == 851123);
  }
  
  LIBSEMIGROUPS_TEST_CASE("BMat8 enum",
                          "119",
                          "filter 8 by special 7 gens 6",
                          "[extreme]") {
    std::vector<BMat8> bmat7_gens(0);
    std::vector<BMat8> bmat8_enum(0);
    std::ifstream      f("bmat_gens_7.txt");
    std::string        line;
    while (std::getline(f, line)) {
      bmat7_gens.push_back(BMat8(std::stoul(line) | 1));
    }
    f.close();
    f = std::ifstream("bmat_filtered_8_7.txt");
    while (std::getline(f, line)) {
      bmat8_enum.push_back(BMat8(std::stoul(line)));
    }
    f.close();
    REQUIRE(bmat8_enum.size() == 851123);

    std::cout << "finished reading!" << std::endl;

    using Perm = typename Perm<8>::type;
    const std::vector<Perm> S8_gens
        = {Perm({1, 2, 3, 4, 5, 6, 7, 0}), Perm({1, 0, 2, 3, 4, 5, 6, 7})};

    FroidurePin<Perm> S8(S8_gens);
    REQUIRE(S8.size() == 40320);
    std::cout << "finished computing S8!" << std::endl;

    std::vector<BMat8> S8_bmats;
    for (Perm p : S8) {
      S8_bmats.push_back(BMat8_from_perm<8>(p));
    }

    size_t                                     i = 0;
    std::vector<std::vector<std::bitset<256>>> row_spaces(
        257, std::vector<std::bitset<256>>());

    for (BMat8 x : bmat7_gens) {
      std::bitset<256> bitset = row_space_bitset(x);
      size_t           count  = bitset.count();
      if (count != 96) {
        continue;
      }
      row_spaces[count].push_back(bitset);
      for (BMat8 const& A : S8_bmats) {
        row_spaces[count].push_back(row_space_bitset(x * A));
      }
      if (++i % 10000 == 0) {
        std::cout << i << std::endl;
      }
    }

    std::unordered_set<BMat8> canon_gens;
    for (BMat8 bm : bmat7_gens) {
      canon_gens.insert(canonical_BMat8(bm, 8));
    }

    std::cout << "finished producing bitsets!" << std::endl;

    std::vector<BMat8> filtered(0);
    size_t             count = 0;
    for (BMat8 bm : bmat8_enum) {
      std::bitset<256> bitset = row_space_bitset(bm);
      bool             found  = false;
      size_t           size   = bitset.count();
      // TODO: not this
      if (size == 256) {
        found = true;  // get rid of any permutation matrices
      }
      if (canon_gens.find(bm) != canon_gens.end()) {
        found = true;
      }
      for (size_t i = size + 1; i < 257 && !found; ++i) {
        for (size_t j = 0; j < row_spaces[i].size(); ++j) {
          std::bitset<256> bs = row_spaces[i][j];
          if ((bitset | bs) == bs) {
            found = true;
            break;
          }
        }
      }
      if (!found) {
        filtered.push_back(bm);
      }
      if (++count % 10000 == 0) {
        std::cout << "on " << count << ", keeping " << filtered.size()
                  << std::endl;
      }
    }

    std::ofstream o;
    o.open("bmat_filtered_8_8.txt", std::ios::out | std::ios::trunc);
    for (BMat8 i : filtered) {
      o << i.to_int() << "\n";
    }
    o.close();
    REQUIRE(filtered.size() == 793961);
  }
  
  LIBSEMIGROUPS_TEST_CASE("BMat8 enum",
                          "120",
                          "filter 8 by known primes",
                          "[extreme]") {
    std::vector<BMat8> known_primes
        = {
           BMat8({{1, 0, 0, 0, 0, 0, 0, 0},
                  {0, 1, 1, 0, 0, 0, 0, 0},
                  {0, 1, 0, 1, 0, 0, 0, 0},
                  {0, 0, 1, 0, 1, 0, 0, 0},
                  {0, 0, 0, 1, 0, 1, 0, 0},
                  {0, 0, 0, 0, 1, 0, 1, 0},
                  {0, 0, 0, 0, 0, 1, 0, 1},
                  {0, 0, 0, 0, 0, 0, 1, 1}})};

    std::vector<BMat8> bmat8_enum(0);
    std::string        line;
    std::ifstream f("bmat_filtered_8_8.txt");
    while (std::getline(f, line)) {
      bmat8_enum.push_back(BMat8(std::stoul(line)));
    }
    f.close();
    REQUIRE(bmat8_enum.size() == 793961);

    std::cout << "finished reading!" << std::endl;

    using Perm = typename Perm<8>::type;
    const std::vector<Perm> S8_gens
        = {Perm({1, 2, 3, 4, 5, 6, 7, 0}), Perm({1, 0, 2, 3, 4, 5, 6, 7})};

    FroidurePin<Perm> S8(S8_gens);
    REQUIRE(S8.size() == 40320);
    std::cout << "finished computing S8!" << std::endl;

    std::vector<BMat8> S8_bmats;
    for (Perm p : S8) {
      S8_bmats.push_back(BMat8_from_perm<8>(p));
    }

    size_t                                     i = 0;
    std::vector<std::vector<std::bitset<256>>> row_spaces(
        257, std::vector<std::bitset<256>>());

    for (BMat8 x : known_primes) {
      std::bitset<256> bitset = row_space_bitset(x);
      size_t           count  = bitset.count();
      row_spaces[count].push_back(bitset);
      for (BMat8 const& A : S8_bmats) {
        row_spaces[count].push_back(row_space_bitset(x * A));
      }
      if (++i % 10000 == 0) {
        std::cout << i << std::endl;
      }
    }

    std::cout << "finished producing bitsets!" << std::endl;

    std::vector<BMat8> filtered(0);
    size_t             count = 0;
    for (BMat8 bm : bmat8_enum) {
      std::bitset<256> bitset = row_space_bitset(bm);
      bool             found  = false;
      size_t           size   = bitset.count();
      // TODO: not this
      if (size == 256) {
        found = true;  // get rid of any permutation matrices
      }
      for (size_t i = size + 1; i < 257 && !found; ++i) {
        for (size_t j = 0; j < row_spaces[i].size(); ++j) {
          std::bitset<256> bs = row_spaces[i][j];
          if ((bitset | bs) == bs) {
            found = true;
            break;
          }
        }
      }
      if (!found) {
        filtered.push_back(bm);
      }
      if (++count % 100 == 0) {
        std::cout << "on " << count << ", keeping " << filtered.size()
                  << std::endl;
      }
    }
    REQUIRE(filtered.size() == 0);
  }
  
  LIBSEMIGROUPS_TEST_CASE("BMat8 enum",
                          "121",
                          "filter 8 without permuting",
                          "[extreme]") {

    std::vector<BMat8> bmat8_enum(0);
    std::string        line;
    std::ifstream f("bmat_filtered_8_8.txt");
    while (std::getline(f, line)) {
      bmat8_enum.push_back(BMat8(std::stoul(line)));
    }
    f.close();
    REQUIRE(bmat8_enum.size() == 793961);

    std::cout << "finished reading!" << std::endl;

    size_t                                     i = 0;
    std::vector<std::vector<std::bitset<256>>> row_spaces(
        257, std::vector<std::bitset<256>>());

    for (BMat8 x : bmat8_enum) {
      std::bitset<256> bitset = row_space_bitset(x);
      size_t           count  = bitset.count();
      row_spaces[count].push_back(bitset);
      if (++i % 10000 == 0) {
        std::cout << i << std::endl;
      }
    }

    std::cout << "finished producing bitsets!" << std::endl;

    std::vector<BMat8> filtered(0);
    size_t             count = 0;
    for (BMat8 bm : bmat8_enum) {
      std::bitset<256> bitset = row_space_bitset(bm);
      bool             found  = false;
      size_t           size   = bitset.count();
      // TODO: not this
      if (size == 256) {
        found = true;  // get rid of any permutation matrices
      }
      for (size_t i = size + 1; i < 257 && !found; ++i) {
        for (size_t j = 0; j < row_spaces[i].size(); ++j) {
          std::bitset<256> bs = row_spaces[i][j];
          if ((bitset | bs) == bs) {
            found = true;
            break;
          }
        }
      }
      if (!found) {
        filtered.push_back(bm);
      }
      if (++count % 10000 == 0) {
        std::cout << "on " << count << ", keeping " << filtered.size()
                  << std::endl;
      }
    }
    REQUIRE(filtered.size() == 793961);
  }
  
  LIBSEMIGROUPS_TEST_CASE("BMat8 enum",
                          "122",
                          "filter 8 by extended 7 gen columns",
                          "[extreme]") {
    std::vector<BMat8> bmat7_gens(0);
    std::vector<BMat8> bmat8_enum(0);
    std::ifstream      f("bmat_gens_7.txt");
    std::string        line;
    while (std::getline(f, line)) {
      bmat7_gens.push_back(BMat8(std::stoul(line) | 1));
    }
    f.close();
    f = std::ifstream("bmat_filtered_8_8.txt");
    while (std::getline(f, line)) {
      bmat8_enum.push_back(BMat8(std::stoul(line)));
    }
    f.close();
    REQUIRE(bmat8_enum.size() == 793961);

    std::cout << "finished reading!" << std::endl;

    size_t                                     i = 0;
    std::vector<std::vector<std::bitset<256>>> col_spaces(
        257, std::vector<std::bitset<256>>());

    std::vector<std::vector<size_t>> containment_count(257,
                                                       std::vector<size_t>(0));

    for (BMat8 x : bmat7_gens) {
      std::bitset<256> bitset = col_space_bitset(x);
      size_t           count  = bitset.count();
      if (count == 256) {
        continue;
      }
      col_spaces[count].push_back(bitset);
      containment_count[count].push_back(0);
      if (++i % 100 == 0) {
        std::cout << i << std::endl;
      }
    }

    std::cout << "finished producing bitsets!" << std::endl;

    std::vector<BMat8> filtered(0);
    size_t             count = 0;
    for (BMat8 bm : bmat8_enum) {
      std::bitset<256> bitset = col_space_bitset(bm);
      bool             found  = false;
      size_t           size   = bitset.count();
      for (size_t i = size + 1; i < 257 && !found; ++i) {
        for (size_t j = 0; j < col_spaces[i].size(); ++j) {
          std::bitset<256> bs = col_spaces[i][j];
          if ((bitset | bs) == bs) {
            found = true;
            containment_count[i][j]++;
            break;
          }
        }
      }
      if (!found) {
        filtered.push_back(bm);
      }
      if (++count % 100000 == 0) {
        std::cout << "on " << count << ", keeping " << filtered.size()
                  << std::endl;
      }
    }

    for (size_t i = 0; i < 256; ++i) {
      for (size_t j = 0; j < containment_count[i].size(); ++j) {
        std::cout << "col space " << i << "," << j
                  << " killed: " << containment_count[i][j] << std::endl;
      }
    }

    std::ofstream o;
    o.open("bmat_filtered_8_9.txt", std::ios::out | std::ios::trunc);
    for (BMat8 i : filtered) {
      o << i.to_int() << "\n";
    }
    o.close();
    REQUIRE(filtered.size() == 786570);
  }
  
  LIBSEMIGROUPS_TEST_CASE("BMat8 enum",
                          "123",
                          "filter 8 by extended 7 gen columns",
                          "[extreme]") {
    std::vector<BMat8> bmat7_gens(0);
    std::vector<BMat8> bmat8_enum(0);
    std::ifstream      f("bmat_gens_7.txt");
    std::string        line;
    while (std::getline(f, line)) {
      bmat7_gens.push_back(BMat8(std::stoul(line) | 1));
    }
    f.close();
    f = std::ifstream("bmat_filtered_8_9.txt");
    while (std::getline(f, line)) {
      bmat8_enum.push_back(BMat8(std::stoul(line)));
    }
    f.close();
    REQUIRE(bmat8_enum.size() == 786570);

    std::cout << "finished reading!" << std::endl;

    using Perm = typename Perm<8>::type;
    const std::vector<Perm> S8_gens
        = {Perm({1, 2, 3, 4, 5, 6, 7, 0}), Perm({1, 0, 2, 3, 4, 5, 6, 7})};

    FroidurePin<Perm> S8(S8_gens);
    REQUIRE(S8.size() == 40320);
    std::cout << "finished computing S8!" << std::endl;

    std::vector<BMat8> S8_bmats;
    for (Perm p : S8) {
      S8_bmats.push_back(BMat8_from_perm<8>(p));
    }
    
    std::vector<std::vector<std::bitset<256>>> col_spaces(
        257, std::vector<std::bitset<256>>());

    for (BMat8 x : bmat7_gens) {
      std::bitset<256> bitset = col_space_bitset(x);
      size_t           count  = bitset.count();
      if (!(count == 144 || count == 160)) {
        continue;
      }
      col_spaces[count].push_back(bitset);
      for (BMat8 const& A : S8_bmats) {
        col_spaces[count].push_back(col_space_bitset(A * x));
      }
    }

    std::cout << "finished producing bitsets!" << std::endl;

    std::vector<BMat8> filtered(0);
    size_t             count = 0;
    for (BMat8 bm : bmat8_enum) {
      std::bitset<256> bitset = col_space_bitset(bm);
      bool             found  = false;
      size_t           size   = bitset.count();
      for (size_t i = size + 1; i < 257 && !found; ++i) {
        for (size_t j = 0; j < col_spaces[i].size(); ++j) {
          std::bitset<256> bs = col_spaces[i][j];
          if ((bitset | bs) == bs) {
            found = true;
            break;
          }
        }
      }
      if (!found) {
        filtered.push_back(bm);
      }
      if (++count % 1000 == 0) {
        std::cout << "on " << count << ", keeping " << filtered.size()
                  << std::endl;
      }
    }

    std::ofstream o;
    o.open("bmat_filtered_8_10.txt", std::ios::out | std::ios::trunc);
    for (BMat8 i : filtered) {
      o << i.to_int() << "\n";
    }
    o.close();
    REQUIRE(filtered.size() == 520113);
  }
  
  LIBSEMIGROUPS_TEST_CASE("BMat8 enum",
                          "124",
                          "filter 8 by extended 7 gen columns 2",
                          "[extreme]") {
    std::vector<BMat8> bmat7_gens(0);
    std::vector<BMat8> bmat8_enum(0);
    std::ifstream      f("bmat_gens_7.txt");
    std::string        line;
    while (std::getline(f, line)) {
      bmat7_gens.push_back(BMat8(std::stoul(line) | 1));
    }
    f.close();
    f = std::ifstream("bmat_filtered_8_10.txt");
    while (std::getline(f, line)) {
      bmat8_enum.push_back(BMat8(std::stoul(line)));
    }
    f.close();
    REQUIRE(bmat8_enum.size() == 520113);

    std::cout << "finished reading!" << std::endl;

    using Perm = typename Perm<8>::type;
    const std::vector<Perm> S8_gens
        = {Perm({1, 2, 3, 4, 5, 6, 7, 0}), Perm({1, 0, 2, 3, 4, 5, 6, 7})};

    FroidurePin<Perm> S8(S8_gens);
    REQUIRE(S8.size() == 40320);
    std::cout << "finished computing S8!" << std::endl;

    std::vector<BMat8> S8_bmats;
    for (Perm p : S8) {
      S8_bmats.push_back(BMat8_from_perm<8>(p));
    }
    
    std::vector<std::vector<std::bitset<256>>> col_spaces(
        257, std::vector<std::bitset<256>>());

    for (BMat8 x : bmat7_gens) {
      std::bitset<256> bitset = col_space_bitset(x);
      size_t           count  = bitset.count();
      if (!(count == 136 || count == 128)) {
        continue;
      }
      col_spaces[count].push_back(bitset);
      for (BMat8 const& A : S8_bmats) {
        col_spaces[count].push_back(col_space_bitset(A * x));
      }
    }

    std::cout << "finished producing bitsets!" << std::endl;

    std::vector<BMat8> filtered(0);
    size_t             count = 0;
    for (BMat8 bm : bmat8_enum) {
      std::bitset<256> bitset = col_space_bitset(bm);
      bool             found  = false;
      size_t           size   = bitset.count();
      for (size_t i = size + 1; i < 257 && !found; ++i) {
        for (size_t j = 0; j < col_spaces[i].size(); ++j) {
          std::bitset<256> bs = col_spaces[i][j];
          if ((bitset | bs) == bs) {
            found = true;
            break;
          }
        }
      }
      if (!found) {
        filtered.push_back(bm);
      }
      if (++count % 1000 == 0) {
        std::cout << "on " << count << ", keeping " << filtered.size()
                  << std::endl;
      }
    }

    std::ofstream o;
    o.open("bmat_filtered_8_11.txt", std::ios::out | std::ios::trunc);
    for (BMat8 i : filtered) {
      o << i.to_int() << "\n";
    }
    o.close();
    REQUIRE(filtered.size() == 520113);
  }

  LIBSEMIGROUPS_TEST_CASE("BMat8 enum",
                          "130",
                          "bmat7 gens row space sizes",
                          "[extreme]") {
    auto               rg = ReportGuard();
    std::vector<BMat8> bmat_enum;
    std::ifstream      f("bmat_gens_7.txt");
    std::string        line;
    while (std::getline(f, line)) {
      bmat_enum.push_back(BMat8(std::stoul(line)));
    }
    f.close();
    std::vector<size_t> row_space_count(129, 0);
    for (BMat8 x : bmat_enum) {
      row_space_count[x.row_space_size()]++;
    }
    for (size_t i = 0; i < 129; ++i) {
      std::cout << "generator row spaces of size " << i << ": "
                << row_space_count[i] << std::endl;
    }
  }

  LIBSEMIGROUPS_TEST_CASE("BMat8 enum",
                          "131",
                          "bmat8 row spaces",
                          "[extreme]") {
    std::vector<BMat8> bmat_enum;
    std::ifstream      f("bmat_trim_enum_6.txt");
    std::string        line;
    while (std::getline(f, line)) {
      bmat_enum.push_back(BMat8(std::stoul(line)));
    }
    f.close();

    std::ofstream o;
    o.open("row_spaces_6.txt", std::ios::out | std::ios::trunc);
    for (BMat8 bm : bmat_enum) {
      std::vector<uint8_t> row_space = row_space_vector(bm);
      o << "{";
      for (uint8_t x: row_space) {
        o << +x << ", ";
      }
      o << "0}" << std::endl;;
    }
    o.close();
  }
  
  LIBSEMIGROUPS_TEST_CASE("BMat8 enum",
                          "132",
                          "bmat8 row space digraphs 5",
                          "[quick]") {
    std::vector<BMat8> bmat_enum;
    std::ifstream      f("bmat_trim_enum_5.txt");
    std::string        line;
    while (std::getline(f, line)) {
      bmat_enum.push_back(BMat8(std::stoul(line)));
    }
    f.close();

    std::ofstream o;
    o.open("row_spaces_digraphs_5.txt", std::ios::out | std::ios::trunc);
    o << "row_digraphs_5 := [" << std::endl;
    for (auto it = bmat_enum.cbegin(); it < bmat_enum.cend(); it++) {
      std::vector<uint8_t> row_space = row_space_vector(*it);
      o << "Digraph([";
      for (size_t i = 0; i < row_space.size(); ++i) {
        o << "[";
        for (size_t j = 0; j < row_space.size(); ++j) {
          if ((row_space[i] | row_space[j]) == row_space[i]) {
            o << j + 1 << ", ";
          }
        }
        o << row_space.size() + 1 << "], ";
      }
      o << "[" << row_space.size() + 1 << "]])";
      if (it != bmat_enum.cend() - 1) {
        o << ",";
      }
      o << std::endl;
    }
    o << "];;" << std::endl;
    o.close();
  }
  
  LIBSEMIGROUPS_TEST_CASE("BMat8 enum",
                          "133",
                          "bmat8 row space digraphs 6",
                          "[quick]") {
    std::vector<BMat8> bmat_enum;
    std::ifstream      f("bmat_trim_enum_6.txt");
    std::string        line;
    while (std::getline(f, line)) {
      bmat_enum.push_back(BMat8(std::stoul(line)));
    }
    f.close();

    std::ofstream o;
    o.open("row_spaces_digraphs_6.txt", std::ios::out | std::ios::trunc);
    o << "row_digraphs_6 := [" << std::endl;
    for (auto it = bmat_enum.cbegin(); it < bmat_enum.cend(); it++) {
      std::vector<uint8_t> row_space = row_space_vector(*it);
      o << "Digraph([";
      for (size_t i = 0; i < row_space.size(); ++i) {
        o << "[";
        for (size_t j = 0; j < row_space.size(); ++j) {
          if ((row_space[i] | row_space[j]) == row_space[i]) {
            o << j + 1 << ", ";
          }
        }
        o << row_space.size() + 1 << "], ";
      }
      o << "[" << row_space.size() + 1 << "]])";
      if (it != bmat_enum.cend() - 1) {
        o << ",";
      }
      o << std::endl;
    }
    o << "];;" << std::endl;
    o.close();
  }
  
  LIBSEMIGROUPS_TEST_CASE("BMat8 enum",
                          "134",
                          "bmat8 row space digraphs 7",
                          "[quick]") {
    std::vector<BMat8> bmat_enum;
    std::ifstream      f("bmat_trim_enum_7.txt");
    std::string        line;
    while (std::getline(f, line)) {
      bmat_enum.push_back(BMat8(std::stoul(line)));
    }
    f.close();

    std::ofstream o;
    o.open("row_spaces_digraphs_7.txt", std::ios::out | std::ios::trunc);
    o << "row_digraphs_7 := [" << std::endl;
    for (auto it = bmat_enum.cbegin(); it < bmat_enum.cend(); it++) {
      std::vector<uint8_t> row_space = row_space_vector(*it);
      o << "Digraph([";
      for (size_t i = 0; i < row_space.size(); ++i) {
        o << "[";
        for (size_t j = 0; j < row_space.size(); ++j) {
          if ((row_space[i] | row_space[j]) == row_space[i]) {
            o << j + 1 << ", ";
          }
        }
        o << row_space.size() + 1 << "], ";
      }
      o << "[" << row_space.size() + 1 << "]])";
      if (it != bmat_enum.cend() - 1) {
        o << ",";
      }
      o << std::endl;
    }
    o << "];;" << std::endl;
    o.close();
  }
  
  LIBSEMIGROUPS_TEST_CASE("BMat8 enum",
                          "135",
                          "bmat8 primes row space digraphs 5",
                          "[quick]") {
    std::vector<BMat8> bmat_enum;
    std::ifstream      f("bmat_gens_test_5.txt");
    std::string        line;
    while (std::getline(f, line)) {
      bmat_enum.push_back(BMat8(std::stoul(line)));
    }
    f.close();

    std::ofstream o;
    o.open("row_spaces_digraphs_primes_5.txt", std::ios::out | std::ios::trunc);
    o << "row_digraphs_primes_5 := [" << std::endl;
    for (auto it = bmat_enum.cbegin(); it < bmat_enum.cend(); it++) {
      std::vector<uint8_t> row_space = row_space_vector(*it);
      o << "Digraph([";
      for (size_t i = 0; i < row_space.size(); ++i) {
        o << "[";
        for (size_t j = 0; j < row_space.size(); ++j) {
          if ((row_space[i] | row_space[j]) == row_space[i]) {
            o << j + 1 << ", ";
          }
        }
        o << row_space.size() + 1 << "], ";
      }
      o << "[" << row_space.size() + 1 << "]])";
      if (it != bmat_enum.cend() - 1) {
        o << ",";
      }
      o << std::endl;
    }
    o << "];;" << std::endl;
    o.close();
  }
  
  LIBSEMIGROUPS_TEST_CASE("BMat8 enum",
                          "136",
                          "bmat8 primes row space digraphs 6",
                          "[quick]") {
    std::vector<BMat8> bmat_enum;
    std::ifstream      f("bmat_gens_test_6.txt");
    std::string        line;
    while (std::getline(f, line)) {
      bmat_enum.push_back(BMat8(std::stoul(line)));
    }
    f.close();

    std::ofstream o;
    o.open("row_spaces_digraphs_primes_6.txt", std::ios::out | std::ios::trunc);
    o << "row_digraphs_primes_6 := [" << std::endl;
    for (auto it = bmat_enum.cbegin(); it < bmat_enum.cend(); it++) {
      std::vector<uint8_t> row_space = row_space_vector(*it);
      o << "Digraph([";
      for (size_t i = 0; i < row_space.size(); ++i) {
        o << "[";
        for (size_t j = 0; j < row_space.size(); ++j) {
          if ((row_space[i] | row_space[j]) == row_space[i]) {
            o << j + 1 << ", ";
          }
        }
        o << row_space.size() + 1 << "], ";
      }
      o << "[" << row_space.size() + 1 << "]])";
      if (it != bmat_enum.cend() - 1) {
        o << ",";
      }
      o << std::endl;
    }
    o << "];;" << std::endl;
    o.close();
  }
  
  LIBSEMIGROUPS_TEST_CASE("BMat8 enum",
                          "137",
                          "bmat8 primes row space digraphs 7",
                          "[quick]") {
    std::vector<BMat8> bmat_enum;
    std::ifstream      f("bmat_gens_test_7.txt");
    std::string        line;
    while (std::getline(f, line)) {
      bmat_enum.push_back(BMat8(std::stoul(line)));
    }
    f.close();

    std::ofstream o;
    o.open("row_spaces_digraphs_primes_7.txt", std::ios::out | std::ios::trunc);
    o << "row_digraphs_primes_7 := [" << std::endl;
    for (auto it = bmat_enum.cbegin(); it < bmat_enum.cend(); it++) {
      if (it - bmat_enum.cbegin() < 10) {
        std::cout << *it << std::endl; 
      }
      std::vector<uint8_t> row_space = row_space_vector(*it);
      o << "Digraph([";
      for (size_t i = 0; i < row_space.size(); ++i) {
        o << "[";
        for (size_t j = 0; j < row_space.size(); ++j) {
          if ((row_space[i] | row_space[j]) == row_space[i]) {
            o << j + 1 << ", ";
          }
        }
        o << row_space.size() + 1 << "], ";
      }
      o << "[" << row_space.size() + 1 << "]])";
      if (it != bmat_enum.cend() - 1) {
        o << ",";
      }
      o << std::endl;
    }
    o << "];;" << std::endl;
    o.close();
  }
  
  LIBSEMIGROUPS_TEST_CASE("BMat8 enum", "160", "print matrices", "[quick]") {
    std::cout << BMat8(13853107707017691136) << std::endl;
    std::cout << BMat8(4620710844303409152) << std::endl;
    std::cout << BMat8(4647750068672397312) << std::endl;
    std::cout << BMat8(9241421688590303232) << std::endl;

    std::cout << BMAT8_ELEMENTARY[6].to_int() << std::endl;
    std::cout << BMAT8_ONES[5].to_int() << std::endl;
  }

  LIBSEMIGROUPS_TEST_CASE("BMat8 enum", "170", "prime extensions", "[quick]") {

    std::vector<BMat8> bmat_enum;
    std::ifstream      f("bmat_gens_7.txt");
    std::string        line;
    while (std::getline(f, line)) {
      bmat_enum.push_back(BMat8(std::stoul(line)));
    }
    f.close();
    BMat8 x = bmat_enum[1021];
    for (BMat8 &bm : simple_prime_extensions(x, 7)) {
      std::cout << bm << std::endl;
    }
  }
  
  LIBSEMIGROUPS_TEST_CASE("BMat8 enum",
                          "171",
                          "filter 8 by special 7 gens simple extensions",
                          "[extreme]") {
    std::vector<BMat8> bmat7_gens(0);
    std::vector<BMat8> bmat8_enum(0);
    std::ifstream      f("bmat_gens_7.txt");
    std::string        line;
    while (std::getline(f, line)) {
      bmat7_gens.push_back(BMat8(std::stoul(line)));
    }
    f.close();
    f = std::ifstream("bmat_filtered_8_11.txt");
    while (std::getline(f, line)) {
      bmat8_enum.push_back(BMat8(std::stoul(line)));
    }
    f.close();
    REQUIRE(bmat8_enum.size() == 484949);

    std::cout << "finished reading!" << std::endl;

    using Perm = typename Perm<8>::type;
    const std::vector<Perm> S8_gens
        = {Perm({1, 2, 3, 4, 5, 6, 7, 0}), Perm({1, 0, 2, 3, 4, 5, 6, 7})};

    FroidurePin<Perm> S8(S8_gens);
    REQUIRE(S8.size() == 40320);
    std::cout << "finished computing S8!" << std::endl;

    std::vector<BMat8> S8_bmats;
    for (Perm p : S8) {
      S8_bmats.push_back(BMat8_from_perm<8>(p));
    }

    size_t                                     i = 0;
    std::vector<std::vector<std::bitset<256>>> col_spaces(
        257, std::vector<std::bitset<256>>(0));
    std::vector<std::vector<size_t>> containment_count(257,
                                                       std::vector<size_t>(0));

    for (BMat8 gen : bmat7_gens) {
      for (BMat8 x : simple_prime_extensions(gen, 7)) {
        std::bitset<256> bitset = col_space_bitset(x);
        size_t           count  = bitset.count();
        if (!(count == 120 || count == 132 || count == 144 || count == 160)) {
          continue;
        }
        col_spaces[count].push_back(bitset);
        containment_count[count].push_back(0);
        for (BMat8 const& A : S8_bmats) {
          col_spaces[count].push_back(col_space_bitset(A * x));
          containment_count[count].push_back(0);
        }
        if (++i % 10000 == 0) {
          std::cout << i << std::endl;
        }
      }
    }
    
    std::cout << "finished producing bitsets!" << std::endl;

    std::vector<BMat8> filtered(0);
    size_t             count = 0;
    for (BMat8 bm : bmat8_enum) {
      std::bitset<256> bitset = col_space_bitset(bm);
      bool             found  = false;
      size_t           size   = bitset.count();
      for (size_t i = size + 1; i < 257 && !found; ++i) {
        for (size_t j = 0; j < col_spaces[i].size(); ++j) {
          std::bitset<256> bs = col_spaces[i][j];
          if ((bitset | bs) == bs) {
            found = true;
            std::cout << i <<", " << j << std::endl;
            containment_count[i][j]++;
            break;
          }
        }
      }
      if (!found) {
        filtered.push_back(bm);
      }
      if (++count % 10000 == 0) {
        std::cout << "on " << count << ", keeping " << filtered.size()
                  << std::endl;
      }
    }

    for (size_t i = 0; i < 256; ++i) {
      for (size_t j = 0; j < containment_count[i].size(); ++j) {
        if (containment_count[i][j]) {
          std::cout << "col space " << i << "," << j
                    << " killed: " << containment_count[i][j] << std::endl;
        }
      }
    }

    std::ofstream o;
    o.open("bmat_filtered_8_12.txt", std::ios::out | std::ios::trunc);
    for (BMat8 i : filtered) {
      o << i.to_int() << "\n";
    }
    o.close();
    REQUIRE(filtered.size() == 472207);
  }
  
  LIBSEMIGROUPS_TEST_CASE("BMat8 enum",
                          "172",
                          "filter 8 by special 7 gens simple extensions 2",
                          "[extreme]") {
    std::vector<BMat8> bmat7_gens(0);
    std::vector<BMat8> bmat8_enum(0);
    std::ifstream      f("bmat_gens_7.txt");
    std::string        line;
    while (std::getline(f, line)) {
      bmat7_gens.push_back(BMat8(std::stoul(line)));
    }
    f.close();
    f = std::ifstream("bmat_filtered_8_12.txt");
    while (std::getline(f, line)) {
      bmat8_enum.push_back(BMat8(std::stoul(line)));
    }
    f.close();
    REQUIRE(bmat8_enum.size() == 472207);

    std::cout << "finished reading!" << std::endl;

    using Perm = typename Perm<8>::type;
    const std::vector<Perm> S8_gens
        = {Perm({1, 2, 3, 4, 5, 6, 7, 0}), Perm({1, 0, 2, 3, 4, 5, 6, 7})};

    FroidurePin<Perm> S8(S8_gens);
    REQUIRE(S8.size() == 40320);
    std::cout << "finished computing S8!" << std::endl;

    std::vector<BMat8> S8_bmats;
    for (Perm p : S8) {
      S8_bmats.push_back(BMat8_from_perm<8>(p));
    }

    size_t                                     i = 0;
    std::vector<std::vector<std::bitset<256>>> col_spaces(
        257, std::vector<std::bitset<256>>(0));
    std::vector<std::vector<size_t>> containment_count(257,
                                                       std::vector<size_t>(0));

    for (BMat8 gen : bmat7_gens) {
      for (BMat8 x : simple_prime_extensions(gen, 7)) {
        std::bitset<256> bitset = col_space_bitset(x);
        size_t           count  = bitset.count();
        if (!(count == 192 || count == 136 || count == 128)) {
          continue;
        }
        col_spaces[count].push_back(bitset);
        containment_count[count].push_back(0);
        for (BMat8 const& A : S8_bmats) {
          col_spaces[count].push_back(col_space_bitset(A * x));
          containment_count[count].push_back(0);
        }
        if (++i % 1000 == 0) {
          std::cout << i << std::endl;
        }
      }
    }
    
    std::cout << "finished producing bitsets!" << std::endl;

    std::vector<BMat8> filtered(0);
    size_t             count = 0;
    for (BMat8 bm : bmat8_enum) {
      std::bitset<256> bitset = col_space_bitset(bm);
      bool             found  = false;
      size_t           size   = bitset.count();
      for (size_t i = size + 1; i < 257 && !found; ++i) {
        for (size_t j = 0; j < col_spaces[i].size(); ++j) {
          std::bitset<256> bs = col_spaces[i][j];
          if ((bitset | bs) == bs) {
            found = true;
            std::cout << i <<", " << j << std::endl;
            containment_count[i][j]++;
            break;
          }
        }
      }
      if (!found) {
        filtered.push_back(bm);
      }
      if (++count % 10000 == 0) {
        std::cout << "on " << count << ", keeping " << filtered.size()
                  << std::endl;
      }
    }

    for (size_t i = 0; i < 256; ++i) {
      for (size_t j = 0; j < containment_count[i].size(); ++j) {
        if (containment_count[i][j]) {
          std::cout << "col space " << i << "," << j
                    << " killed: " << containment_count[i][j] << std::endl;
        }
      }
    }

    std::ofstream o;
    o.open("bmat_filtered_8_13.txt", std::ios::out | std::ios::trunc);
    for (BMat8 i : filtered) {
      o << i.to_int() << "\n";
    }
    o.close();
    REQUIRE(filtered.size() == 472207);
  }
  
  LIBSEMIGROUPS_TEST_CASE("BMat8 enum",
                          "173",
                          "filter 8 by special 7 gens simple extensions",
                          "[extreme]") {
    std::vector<BMat8> bmat7_gens(0);
    std::vector<BMat8> bmat8_enum(0);
    std::ifstream      f("bmat_gens_7.txt");
    std::string        line;
    while (std::getline(f, line)) {
      bmat7_gens.push_back(BMat8(std::stoul(line)));
    }
    f.close();
    f = std::ifstream("bmat_filtered_8_12.txt");
    while (std::getline(f, line)) {
      bmat8_enum.push_back(BMat8(std::stoul(line)));
    }
    f.close();
    REQUIRE(bmat8_enum.size() == 472207);

    std::cout << "finished reading!" << std::endl;

    using Perm = typename Perm<8>::type;
    const std::vector<Perm> S8_gens
        = {Perm({1, 2, 3, 4, 5, 6, 7, 0}), Perm({1, 0, 2, 3, 4, 5, 6, 7})};

    FroidurePin<Perm> S8(S8_gens);
    REQUIRE(S8.size() == 40320);
    std::cout << "finished computing S8!" << std::endl;

    std::vector<BMat8> S8_bmats;
    for (Perm p : S8) {
      S8_bmats.push_back(BMat8_from_perm<8>(p));
    }

    size_t                                     i = 0;
    std::vector<std::vector<std::bitset<256>>> row_spaces(
        257, std::vector<std::bitset<256>>(0));
    std::vector<std::vector<size_t>> containment_count(257,
                                                       std::vector<size_t>(0));

    for (BMat8 gen : bmat7_gens) {
      for (BMat8 x : simple_prime_extensions(gen, 7)) {
        std::bitset<256> bitset = row_space_bitset(x);
        size_t           count  = bitset.count();
        if (count == 256) {
          continue;
        }
        row_spaces[count].push_back(bitset);
        containment_count[count].push_back(0);
        if (++i % 10000 == 0) {
          std::cout << i << std::endl;
        }
      }
    }

    std::cout << "finished producing bitsets!" << std::endl;

    std::vector<BMat8> filtered(0);
    size_t             count = 0;
    for (BMat8 bm : bmat8_enum) {
      std::bitset<256> bitset = row_space_bitset(bm);
      bool             found  = false;
      size_t           size   = bitset.count();
      for (size_t i = size + 1; i < 257 && !found; ++i) {
        for (size_t j = 0; j < row_spaces[i].size(); ++j) {
          std::bitset<256> bs = row_spaces[i][j];
          if ((bitset | bs) == bs) {
            found = true;
            std::cout << i <<", " << j << std::endl;
            containment_count[i][j]++;
            break;
          }
        }
      }
      if (!found) {
        filtered.push_back(bm);
      }
      if (++count % 10000 == 0) {
        std::cout << "on " << count << ", keeping " << filtered.size()
                  << std::endl;
      }
    }

    for (size_t i = 0; i < 256; ++i) {
      for (size_t j = 0; j < containment_count[i].size(); ++j) {
        if (containment_count[i][j]) {
          std::cout << "row space " << i << "," << j
                    << " killed: " << containment_count[i][j] << std::endl;
        }
      }
    }
    REQUIRE(filtered.size() == 472207);
  }
  
  LIBSEMIGROUPS_TEST_CASE("BMat8 enum",
                          "174",
                          "filter 8 by special 7 gens simple col extensions",
                          "[extreme]") {
    std::vector<BMat8> bmat7_gens(0);
    std::vector<BMat8> bmat8_enum(0);
    std::ifstream      f("bmat_gens_7.txt");
    std::string        line;
    while (std::getline(f, line)) {
      bmat7_gens.push_back(BMat8(std::stoul(line)));
    }
    f.close();
    f = std::ifstream("bmat_filtered_8_12.txt");
    while (std::getline(f, line)) {
      bmat8_enum.push_back(BMat8(std::stoul(line)));
    }
    f.close();
    REQUIRE(bmat8_enum.size() == 472207);

    std::cout << "finished reading!" << std::endl;

    using Perm = typename Perm<8>::type;
    const std::vector<Perm> S8_gens
        = {Perm({1, 2, 3, 4, 5, 6, 7, 0}), Perm({1, 0, 2, 3, 4, 5, 6, 7})};

    FroidurePin<Perm> S8(S8_gens);
    REQUIRE(S8.size() == 40320);
    std::cout << "finished computing S8!" << std::endl;

    std::vector<BMat8> S8_bmats;
    for (Perm p : S8) {
      S8_bmats.push_back(BMat8_from_perm<8>(p));
    }

    size_t                                     i = 0;
    std::vector<std::vector<std::bitset<256>>> row_spaces(
        257, std::vector<std::bitset<256>>(0));
    std::vector<std::vector<size_t>> containment_count(257,
                                                       std::vector<size_t>(0));

    for (BMat8 gen : bmat7_gens) {
      for (BMat8 x : simple_prime_extensions(gen.transpose(), 7)) {
        std::bitset<256> bitset = row_space_bitset(x.transpose());
        size_t           count  = bitset.count();
        if (count == 256) {
          continue;
        }
        row_spaces[count].push_back(bitset);
        containment_count[count].push_back(0);
        if (++i % 10000 == 0) {
          std::cout << i << std::endl;
        }
      }
    }

    std::cout << "finished producing bitsets!" << std::endl;

    std::vector<BMat8> filtered(0);
    size_t             count = 0;
    for (BMat8 bm : bmat8_enum) {
      std::bitset<256> bitset = row_space_bitset(bm);
      bool             found  = false;
      size_t           size   = bitset.count();
      for (size_t i = size + 1; i < 257 && !found; ++i) {
        for (size_t j = 0; j < row_spaces[i].size(); ++j) {
          std::bitset<256> bs = row_spaces[i][j];
          if ((bitset | bs) == bs) {
            found = true;
            std::cout << i <<", " << j << std::endl;
            containment_count[i][j]++;
            break;
          }
        }
      }
      if (!found) {
        filtered.push_back(bm);
      }
      if (++count % 10000 == 0) {
        std::cout << "on " << count << ", keeping " << filtered.size()
                  << std::endl;
      }
    }

    for (size_t i = 0; i < 256; ++i) {
      for (size_t j = 0; j < containment_count[i].size(); ++j) {
        if (containment_count[i][j]) {
          std::cout << "row space " << i << "," << j
                    << " killed: " << containment_count[i][j] << std::endl;
        }
      }
    }
    REQUIRE(filtered.size() == 472207);
  }
  
  LIBSEMIGROUPS_TEST_CASE("BMat8 enum",
                          "175",
                          "filter 8 by special 7 gens simple col extensions",
                          "[extreme]") {
    std::vector<BMat8> bmat7_gens(0);
    std::vector<BMat8> bmat8_enum(0);
    std::ifstream      f("bmat_gens_7.txt");
    std::string        line;
    while (std::getline(f, line)) {
      bmat7_gens.push_back(BMat8(std::stoul(line)));
    }
    f.close();
    f = std::ifstream("bmat_filtered_8_12.txt");
    while (std::getline(f, line)) {
      bmat8_enum.push_back(BMat8(std::stoul(line)));
    }
    f.close();
    REQUIRE(bmat8_enum.size() == 472207);

    std::cout << "finished reading!" << std::endl;

    using Perm = typename Perm<8>::type;
    const std::vector<Perm> S8_gens
        = {Perm({1, 2, 3, 4, 5, 6, 7, 0}), Perm({1, 0, 2, 3, 4, 5, 6, 7})};

    FroidurePin<Perm> S8(S8_gens);
    REQUIRE(S8.size() == 40320);
    std::cout << "finished computing S8!" << std::endl;

    std::vector<BMat8> S8_bmats;
    for (Perm p : S8) {
      S8_bmats.push_back(BMat8_from_perm<8>(p));
    }

    size_t                                     i = 0;
    std::vector<std::vector<std::bitset<256>>> row_spaces(
        257, std::vector<std::bitset<256>>(0));
    std::vector<std::vector<size_t>> containment_count(257,
                                                       std::vector<size_t>(0));

    for (BMat8 gen : bmat7_gens) {
      for (BMat8 x : simple_prime_extensions(gen, 7)) {
        std::bitset<256> bitset = row_space_bitset(x);
        size_t           count  = bitset.count();
        if (!(count == 86)) {
          continue;
        }
        row_spaces[count].push_back(bitset);
        containment_count[count].push_back(0);
        for (BMat8 const& A : S8_bmats) {
          row_spaces[count].push_back(row_space_bitset(x * A));
          containment_count[count].push_back(0);
        }
        if (++i % 1000 == 0) {
          std::cout << i << std::endl;
        }
      }
    }
    
    std::cout << "finished producing bitsets!" << std::endl;

    std::vector<BMat8> filtered(0);
    size_t             count = 0;
    for (BMat8 bm : bmat8_enum) {
      std::bitset<256> bitset = row_space_bitset(bm);
      bool             found  = false;
      size_t           size   = bitset.count();
      for (size_t i = size + 1; i < 257 && !found; ++i) {
        for (size_t j = 0; j < row_spaces[i].size(); ++j) {
          std::bitset<256> bs = row_spaces[i][j];
          if ((bitset | bs) == bs) {
            found = true;
            std::cout << i <<", " << j << std::endl;
            containment_count[i][j]++;
            break;
          }
        }
      }
      if (!found) {
        filtered.push_back(bm);
      }
      if (++count % 1000 == 0) {
        std::cout << "on " << count << ", keeping " << filtered.size()
                  << std::endl;
      }
    }

    for (size_t i = 0; i < 256; ++i) {
      for (size_t j = 0; j < containment_count[i].size(); ++j) {
        if (containment_count[i][j]) {
          std::cout << "row space " << i << "," << j
                    << " killed: " << containment_count[i][j] << std::endl;
        }
      }
    }

    std::ofstream o;
    o.open("bmat_filtered_8_14.txt", std::ios::out | std::ios::trunc);
    for (BMat8 i : filtered) {
      o << i.to_int() << "\n";
    }
    o.close();
    REQUIRE(filtered.size() == 472207);
  }
  
  LIBSEMIGROUPS_TEST_CASE("BMat8 enum",
                          "176",
                          "filter 8",
                          "[extreme]") {
    std::vector<BMat8> bmat7_gens(0);
    std::vector<BMat8> bmat8_enum(0);
    std::ifstream      f("bmat_gens_7.txt");
    std::string        line;
    while (std::getline(f, line)) {
      bmat7_gens.push_back(BMat8(std::stoul(line)));
    }
    f.close();
    f = std::ifstream("bmat_filtered_8_14.txt");
    while (std::getline(f, line)) {
      bmat8_enum.push_back(BMat8(std::stoul(line)));
    }
    f.close();
    REQUIRE(bmat8_enum.size() == 472207);

    std::cout << "finished reading!" << std::endl;

    using Perm = typename Perm<8>::type;
    const std::vector<Perm> S8_gens
        = {Perm({1, 2, 3, 4, 5, 6, 7, 0}), Perm({1, 0, 2, 3, 4, 5, 6, 7})};

    FroidurePin<Perm> S8(S8_gens);
    REQUIRE(S8.size() == 40320);
    std::cout << "finished computing S8!" << std::endl;

    std::vector<BMat8> S8_bmats;
    for (Perm p : S8) {
      S8_bmats.push_back(BMat8_from_perm<8>(p));
    }

    size_t                                     i = 0;
    std::vector<std::vector<std::bitset<256>>> row_spaces(
        257, std::vector<std::bitset<256>>(0));
    std::vector<std::vector<size_t>> containment_count(257,
                                                       std::vector<size_t>(0));

    for (BMat8 x : bmat8_enum) {
        std::bitset<256> bitset = row_space_bitset(x);
        size_t           count  = bitset.count();
        if (!(count == 86)) {
          continue;
        }
        row_spaces[count].push_back(bitset);
        containment_count[count].push_back(0);
        for (BMat8 const& A : S8_bmats) {
          row_spaces[count].push_back(row_space_bitset(x * A));
          containment_count[count].push_back(0);
        }
        if (++i % 1000 == 0) {
          std::cout << i << std::endl;
        }
      }
    
    std::cout << "finished producing bitsets!" << std::endl;

    std::vector<BMat8> filtered(0);
    size_t             count = 0;
    for (BMat8 bm : bmat8_enum) {
      std::bitset<256> bitset = row_space_bitset(bm);
      bool             found  = false;
      size_t           size   = bitset.count();
      for (size_t i = size + 1; i < 257 && !found; ++i) {
        for (size_t j = 0; j < row_spaces[i].size(); ++j) {
          std::bitset<256> bs = row_spaces[i][j];
          if ((bitset | bs) == bs) {
            found = true;
            std::cout << i <<", " << j << std::endl;
            containment_count[i][j]++;
            break;
          }
        }
      }
      if (!found) {
        filtered.push_back(bm);
      }
      if (++count % 1000 == 0) {
        std::cout << "on " << count << ", keeping " << filtered.size()
                  << std::endl;
      }
    }

    for (size_t i = 0; i < 256; ++i) {
      for (size_t j = 0; j < containment_count[i].size(); ++j) {
        if (containment_count[i][j]) {
          std::cout << "row space " << i << "," << j
                    << " killed: " << containment_count[i][j] << std::endl;
        }
      }
    }

    std::ofstream o;
    o.open("bmat_filtered_8_15.txt", std::ios::out | std::ios::trunc);
    for (BMat8 i : filtered) {
      o << i.to_int() << "\n";
    }
    o.close();
    REQUIRE(filtered.size() == 472142);
  }
  
  LIBSEMIGROUPS_TEST_CASE("BMat8 enum",
                          "177",
                          "filter 8",
                          "[extreme]") {
    std::vector<BMat8> bmat7_gens(0);
    std::vector<BMat8> bmat8_enum(0);
    std::ifstream      f("bmat_gens_7.txt");
    std::string        line;
    while (std::getline(f, line)) {
      bmat7_gens.push_back(BMat8(std::stoul(line)));
    }
    f.close();
    f = std::ifstream("bmat_filtered_8_15.txt");
    while (std::getline(f, line)) {
      bmat8_enum.push_back(BMat8(std::stoul(line)));
    }
    f.close();
    REQUIRE(bmat8_enum.size() == 472142);

    std::cout << "finished reading!" << std::endl;

    using Perm = typename Perm<8>::type;
    const std::vector<Perm> S8_gens
        = {Perm({1, 2, 3, 4, 5, 6, 7, 0}), Perm({1, 0, 2, 3, 4, 5, 6, 7})};

    FroidurePin<Perm> S8(S8_gens);
    REQUIRE(S8.size() == 40320);
    std::cout << "finished computing S8!" << std::endl;

    std::vector<BMat8> S8_bmats;
    for (Perm p : S8) {
      S8_bmats.push_back(BMat8_from_perm<8>(p));
    }

    size_t                                     i = 0;
    std::vector<std::vector<std::bitset<256>>> row_spaces(
        257, std::vector<std::bitset<256>>(0));
    std::vector<std::vector<size_t>> containment_count(257,
                                                       std::vector<size_t>(0));

    for (BMat8 x : bmat8_enum) {
        std::bitset<256> bitset = row_space_bitset(x);
        size_t           count  = bitset.count();
        if (!(count == 80 || count == 84)) {
          continue;
        }
        row_spaces[count].push_back(bitset);
        containment_count[count].push_back(0);
        for (BMat8 const& A : S8_bmats) {
          row_spaces[count].push_back(row_space_bitset(x * A));
          containment_count[count].push_back(0);
        }
        if (++i % 1000 == 0) {
          std::cout << i << std::endl;
        }
      }
    
    std::cout << "finished producing bitsets!" << std::endl;

    std::vector<BMat8> filtered(0);
    size_t             count = 0;
    for (BMat8 bm : bmat8_enum) {
      std::bitset<256> bitset = row_space_bitset(bm);
      bool             found  = false;
      size_t           size   = bitset.count();
      for (size_t i = size + 1; i < 257 && !found; ++i) {
        for (size_t j = 0; j < row_spaces[i].size(); ++j) {
          std::bitset<256> bs = row_spaces[i][j];
          if ((bitset | bs) == bs) {
            found = true;
            std::cout << i <<", " << j << std::endl;
            containment_count[i][j]++;
            break;
          }
        }
      }
      if (!found) {
        filtered.push_back(bm);
      }
      if (++count % 1000 == 0) {
        std::cout << "on " << count << ", keeping " << filtered.size()
                  << std::endl;
      }
    }

    for (size_t i = 0; i < 256; ++i) {
      for (size_t j = 0; j < containment_count[i].size(); ++j) {
        if (containment_count[i][j]) {
          std::cout << "row space " << i << "," << j
                    << " killed: " << containment_count[i][j] << std::endl;
        }
      }
    }

    std::ofstream o;
    o.open("bmat_filtered_8_16.txt", std::ios::out | std::ios::trunc);
    for (BMat8 i : filtered) {
      o << i.to_int() << "\n";
    }
    o.close();
    REQUIRE(filtered.size() == 472001);
  }
  
  LIBSEMIGROUPS_TEST_CASE("BMat8 enum",
                          "178",
                          "filter 8",
                          "[extreme]") {
    std::vector<BMat8> bmat7_gens(0);
    std::vector<BMat8> bmat8_enum(0);
    std::string        line;
    std::ifstream f("bmat_filtered_8_16.txt");
    while (std::getline(f, line)) {
      bmat8_enum.push_back(BMat8(std::stoul(line)));
    }
    f.close();
    REQUIRE(bmat8_enum.size() == 472001);

    std::cout << "finished reading!" << std::endl;

    using Perm = typename Perm<8>::type;
    const std::vector<Perm> S8_gens
        = {Perm({1, 2, 3, 4, 5, 6, 7, 0}), Perm({1, 0, 2, 3, 4, 5, 6, 7})};

    FroidurePin<Perm> S8(S8_gens);
    REQUIRE(S8.size() == 40320);
    std::cout << "finished computing S8!" << std::endl;

    std::vector<BMat8> S8_bmats;
    for (Perm p : S8) {
      S8_bmats.push_back(BMat8_from_perm<8>(p));
    }

    size_t                                     i = 0;
    std::vector<std::vector<std::bitset<256>>> row_spaces(
        257, std::vector<std::bitset<256>>(0));
    std::vector<std::vector<size_t>> containment_count(257,
                                                       std::vector<size_t>(0));

    for (BMat8 x : bmat8_enum) {
        std::bitset<256> bitset = row_space_bitset(x);
        size_t           count  = bitset.count();
        if (!(count == 82)) {
          continue;
        }
        row_spaces[count].push_back(bitset);
        containment_count[count].push_back(0);
        for (BMat8 const& A : S8_bmats) {
          row_spaces[count].push_back(row_space_bitset(x * A));
          containment_count[count].push_back(0);
        }
        if (++i % 1000 == 0) {
          std::cout << i << std::endl;
        }
      }
    
    std::cout << "finished producing bitsets!" << std::endl;

    std::vector<BMat8> filtered(0);
    size_t             count = 0;
    for (BMat8 bm : bmat8_enum) {
      std::bitset<256> bitset = row_space_bitset(bm);
      bool             found  = false;
      size_t           size   = bitset.count();
      for (size_t i = size + 1; i < 257 && !found; ++i) {
        for (size_t j = 0; j < row_spaces[i].size(); ++j) {
          std::bitset<256> bs = row_spaces[i][j];
          if ((bitset | bs) == bs) {
            found = true;
            std::cout << i <<", " << j << std::endl;
            containment_count[i][j]++;
            break;
          }
        }
      }
      if (!found) {
        filtered.push_back(bm);
      }
      if (++count % 1000 == 0) {
        std::cout << "on " << count << ", keeping " << filtered.size()
                  << std::endl;
      }
    }

    for (size_t i = 0; i < 256; ++i) {
      for (size_t j = 0; j < containment_count[i].size(); ++j) {
        if (containment_count[i][j]) {
          std::cout << "row space " << i << "," << j
                    << " killed: " << containment_count[i][j] << std::endl;
        }
      }
    }

    std::ofstream o;
    o.open("bmat_filtered_8_17.txt", std::ios::out | std::ios::trunc);
    for (BMat8 i : filtered) {
      o << i.to_int() << "\n";
    }
    o.close();
    REQUIRE(filtered.size() == 472001);
  }

  LIBSEMIGROUPS_TEST_CASE("BMat8 enum",
                          "180",
                          "produce extended B7 primes",
                          "[quick]"){
    std::vector<BMat8> bmat7_gens(0);
    std::vector<BMat8> bmat8_enum(0);
    std::ifstream      f("bmat_gens_7.txt");
    std::string        line;
    std::ofstream o;
    o.open("bmat8_filterers.txt", std::ios::out | std::ios::trunc);
    while (std::getline(f, line)) {
      o << BMat8(std::stoul(line) | 1).to_int() << "\n";
    }
    f.close();
    o.close();
  }

  LIBSEMIGROUPS_TEST_CASE("BMat8 enum",
                          "190",
                          "bmat7 gens extended row space sizes",
                          "[extreme]") {
    auto               rg = ReportGuard();
    std::vector<BMat8> bmat_enum;
    std::ifstream      f("bmat_gens_7.txt");
    std::string        line;
    while (std::getline(f, line)) {
      bmat_enum.push_back(BMat8(std::stoul(line)));
    }
    f.close();
    std::vector<size_t> row_space_count(257, 0);
    for (BMat8 x : bmat_enum) {
        row_space_count[x.row_space_size()]++;
    }
    for (size_t i = 0; i < 257; ++i) {
      std::cout << "row spaces of size " << i << ": "
                << row_space_count[i] << std::endl;
    }
  }

  LIBSEMIGROUPS_TEST_CASE("BMat8 enum", 
                          "200",
                          "print elem",
                          "[quick]") {
    for (size_t i = 0; i < 9; ++i) {
    std::cout << std::hex << BMAT8_ELEMENTARY[i].to_int() << std::endl << std::dec;
    }
  }
  
  LIBSEMIGROUPS_TEST_CASE("BMat8 enum", 
                          "300",
                          "count fully indecomposable",
                          "[quick]") {
  
    std::vector<BMat8> bmat8_enum;
    std::ifstream      f("bmat_trim_enum_8.txt");
    std::string        line;
    while (std::getline(f, line)) {
      bmat8_enum.push_back(BMat8(std::stoul(line)));
    }
    f.close();
    REQUIRE(bmat8_enum.size() == 17120845);

    std::cout << "finished reading!" << std::endl;
    
    std::ofstream o;
    o.open("fully_indecomposable_8.txt", std::ios::out | std::ios::trunc);
    size_t count = 0;
    for (BMat8 bm : bmat8_enum) {
      if (is_fully_indecomposable(bm, 8)) {
        count++;
        o << bm.to_int() << "\n";
      }
    }
    o.close();
  }
  
  LIBSEMIGROUPS_TEST_CASE("BMat8 enum", 
                          "301",
                          "count fully indecomposable",
                          "[quick]") {
  
    std::vector<BMat8> bmat8_enum;
    std::ifstream      f("bmat_filtered_8_17.txt");
    std::string        line;
    while (std::getline(f, line)) {
      bmat8_enum.push_back(BMat8(std::stoul(line)));
    }
    f.close();
    REQUIRE(bmat8_enum.size() == 471912);

    std::cout << "finished reading!" << std::endl;
    size_t count = 0;
    for (BMat8 bm : bmat8_enum) {
      count = is_fully_indecomposable(bm, 8) ? count + 1 : count;
    }
    std::cout << +count << " filtered fully indecomposable trim reps of dim 8"
              << std::endl;
  }
}  // namespace libsemigroups
