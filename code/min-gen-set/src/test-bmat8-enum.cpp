#include <bitset>
#include <fstream>
#include <iostream>
#include <unordered_set>  // for unordered_set
#include <vector>         // for vector

#include "catch.hpp"
#include "test-main.hpp"

#include "bmat8.hpp"
#include "element-helper.hpp"
#include "element.hpp"
#include "froidure-pin.hpp"
#include "hpcombi.hpp"
#include "runner.hpp"
#include "string.hpp"
#include "timer.hpp"
#include "transf.hpp"

#include "../extern/bliss-0.73/graph.hh"

namespace libsemigroups {
  typedef bliss::Digraph bliss_digraph;

  bliss::Stats stats;

  bliss_digraph bliss_digraph_from_BMat8(HPCombi::BMat8 bm, size_t dim = 8) {
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

  HPCombi::BMat8 permuted_BMat8(HPCombi::BMat8      bm,
                                size_t              dim,
                                const unsigned int* perm) {
    HPCombi::epu8 row_perm;
    HPCombi::epu8 col_perm;
    for (size_t i = 0; i < dim; ++i) {
      row_perm[perm[i]]             = i;
      col_perm[perm[i + dim] - dim] = i;
    }
    for (size_t i = dim; i < 8; ++i) {
      row_perm[i] = i;
      col_perm[i] = i;
    }

    return bm.row_permuted(row_perm).col_permuted(col_perm);
  }

  HPCombi::BMat8 canonical_BMat8(HPCombi::BMat8 bm, size_t dim) {
    bliss_digraph dg = bliss_digraph_from_BMat8(bm, dim);
    return permuted_BMat8(
        bm, dim, dg.canonical_form(stats, &bliss_hook_function, nullptr));
  }

  bool is_row_reduced(HPCombi::BMat8 bm) {
    return bm.nr_rows() == bm.row_space_basis().nr_rows();
  }

  bool is_col_reduced(HPCombi::BMat8 bm) {
    return is_row_reduced(bm.transpose());
  }

  std::bitset<256> row_space_bitset(HPCombi::BMat8 bm) {
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

  std::bitset<256> col_space_bitset(HPCombi::BMat8 bm) {
    return row_space_bitset(bm.transpose());
  }

  HPCombi::BMat8 right_mult(HPCombi::BMat8 pt, HPCombi::BMat8 x) {
    return pt * x;
  }

  HPCombi::BMat8 left_mult(HPCombi::BMat8 pt, HPCombi::BMat8 x) {
    return x * pt;
  }

  std::vector<uint8_t> row_space_vector(HPCombi::BMat8 bm) {
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

  bool is_row_trim(HPCombi::BMat8 bm, size_t dim = 8) {
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

  bool is_col_trim(HPCombi::BMat8 bm, size_t dim = 8) {
    return is_row_trim(bm.transpose(), dim);
  }

  bool is_trim(HPCombi::BMat8 bm, size_t dim = 8) {
    return is_row_trim(bm, dim) && is_col_trim(bm, dim);
  }

  int nr_ones(uint8_t x) {
    return __builtin_popcount(x);
  }

  // given a HPCombi::BMat8 bm of dimension dim, find all the BMat8s of
  // dimension one higher which are extensions of bm with 0s in the new column
  // and maximal row spaces
  // TODO: obvious optimisations
  std::vector<HPCombi::BMat8> simple_prime_extensions(HPCombi::BMat8 bm,
                                                      size_t         dim) {
    std::vector<uint8_t>        rows = row_space_vector(bm);
    std::vector<HPCombi::BMat8> extensions(0);

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
            HPCombi::BMat8(bm.to_int() | (((new_row << 1) | 1) << (7 - dim))));
      }
    }
    return extensions;
  }

  bool is_fully_indecomposable(HPCombi::BMat8 bm, size_t dim) {
    std::vector<uint8_t> rows = bm.rows();
    for (size_t x = 1; x < (1 << (dim - 1)); ++x) {
      uint8_t res   = 0;
      size_t  count = nr_ones(x);
      for (size_t i = 0; i < dim; ++i) {
        res = ((x >> i) & 1) ? res | rows[i] : res;
      }
      if (static_cast<size_t>(nr_ones(res)) < count + 1) {
        return false;
      }
    }
    return true;
  }

  class BMatEnumerator : public ::libsemigroups::Runner {
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
          _trim(trim),
          _finished(false) {}

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
            _rows[k]          = row;
            auto bm = bmat8_helpers::make<HPCombi::BMat8>(_rows.cbegin(),
                                                          _rows.cend());
            // move the matrix to the right place
            bm = HPCombi::BMat8(bm.to_int() << (8 - _n));
            if ((_trim && is_col_trim(bm, _n))
                || (!_trim && is_col_reduced(bm))) {
              HPCombi::BMat8 canon = canonical_BMat8(bm, _n);
              if (_set.find(canon) == _set.end()) {
                _set.insert(canon);
                _out.push_back(canon);
              }
            }
          }
        next_loop_last_row:
          // uint8_t overflow!!!
          if (row == _max) {
            break;
          }
        }
        if (report()) {
          REPORT_DEFAULT("found %d reps so far, currently on %s \n",
                         _out.size(),
                         detail::to_string(bmat8_helpers::make<HPCombi::BMat8>(
                             _rows.cbegin(), _rows.cend())));
        }
        _rows[k] = 0;
      }
    }

   public:
    void run_impl() {
      for (bool& x : _row_seen) {
        x = false;
      }
      _row_seen[0] = true;

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
      _out.push_back(
          canonical_BMat8(HPCombi::BMat8(static_cast<size_t>(1) << 63), _n));
      _out.push_back(HPCombi::BMat8(0));
      _finished = true;
      report_why_we_stopped();
    }

    bool finished_impl() const {
      return _finished;
    }

    const std::vector<HPCombi::BMat8>& reps() {
      if (!started()) {
        run();
      }
      return _out;
    }

   private:
    size_t                             _n;
    size_t                             _max;
    std::vector<uint8_t>               _rows;
    std::unordered_set<HPCombi::BMat8> _set;
    std::vector<HPCombi::BMat8>        _out;
    std::array<bool, 256>              _row_seen;
    std::vector<std::vector<uint8_t>>  _row_orb_by_row;
    size_t                             _first_row;
    int                                _min_ones;
    bool                               _trim;
    bool                               _finished;
  };

  template <size_t _dim>
  class Filterer : public ::libsemigroups::Runner {
   public:
    Filterer(std::string                 in,
             std::string                 out,
             std::vector<HPCombi::BMat8> filterers,
             bool                        permute,
             bool                        col_filt = false)
        : _in(in),
          _out(out),
          _filtered(0),
          _filterers(std::move(filterers)),
          _permute(permute),
          _finished(false) {
      if (col_filt) {
        _act      = left_mult;
        _bit_func = col_space_bitset;
      } else {
        _act      = right_mult;
        _bit_func = row_space_bitset;
      }
    }
    void run_impl() {
      std::vector<HPCombi::BMat8> bmat_enum;
      std::ifstream               f(_in);
      std::string                 line;
      while (std::getline(f, line)) {
        bmat_enum.push_back(HPCombi::BMat8(std::stoul(line)));
      }
      f.close();

      using Perm = typename PermHelper<_dim>::type;
      HPCombi::epu8 cycle;
      HPCombi::epu8 transposition;
      // TODO: fix for tiny perms
      transposition[0] = 1;
      transposition[1] = 0;
      cycle[0]         = 1;
      for (uint8_t i = 2; i < _dim; ++i) {
        cycle[i - 1]     = i;
        transposition[i] = i;
      }
      cycle[_dim - 1] = 0;

      std::vector<Perm> S_gens = {Perm(cycle), Perm(transposition)};

      FroidurePin<Perm> S(S_gens);
      S.run();

      std::vector<HPCombi::BMat8> S_bmats;
      for (Perm p : S) {
        S_bmats.push_back(
            bmat8_helpers::make<_dim, typename PermHelper<_dim>::type,
                                HPCombi::BMat8>(p));
      }

      // TODO evil copying
      if (_filterers.size() == 0) {
        _filterers = bmat_enum;
        _filterers.push_back(bmat8_helpers::elementary<HPCombi::BMat8>(_dim));
      }

      std::vector<std::vector<std::bitset<256>>> spaces(
          (1 << _dim) + 1, std::vector<std::bitset<256>>());
      for (HPCombi::BMat8 x : _filterers) {
        std::bitset<256> bitset = _bit_func(x);
        size_t           count  = bitset.count();
        spaces[count].push_back(bitset);
        if (_permute) {
          for (HPCombi::BMat8 y : S_bmats) {
            bitset = _bit_func(_act(x, y));
            spaces[count].push_back(bitset);
          }
        }
      }

      for (size_t i = 0; i < bmat_enum.size(); ++i) {
        HPCombi::BMat8   bm     = bmat_enum[i];
        std::bitset<256> bitset = _bit_func(bm);
        bool             found  = false;
        if (bitset.count() == (1 << _dim)) {
          continue;
        }
        for (size_t i = (1 << _dim) - 1; i > bitset.count(); --i) {
          for (std::bitset<256> bs : spaces[i]) {
            if ((bitset | bs) == bs) {
              found = true;
              break;
            }
          }
        }
        if (!found) {
          _filtered.push_back(bm);
        }
        if (report()) {
          REPORT_DEFAULT("On %d out of %d, keeping %d.\n",
                         i + 1,
                         bmat_enum.size(),
                         _filtered.size());
        }
      }

      std::ofstream o;
      o.open(_out, std::ios::out | std::ios::trunc);
      for (HPCombi::BMat8 i : _filtered) {
        o << i.to_int() << "\n";
      }
      o.close();

      _finished = true;
      report_why_we_stopped();
    }

    const std::vector<HPCombi::BMat8>& reps() {
      if (!started()) {
        run();
      }
      return _filtered;
    }

    size_t size() {
      if (!started()) {
        run();
      }
      return _filtered.size();
    }

    bool finished_impl() const {
      return _finished;
    }

   private:
    std::string                                                   _in;
    std::string                                                   _out;
    std::vector<HPCombi::BMat8>                                   _filtered;
    std::vector<HPCombi::BMat8>                                   _filterers;
    bool                                                          _permute;
    std::function<HPCombi::BMat8(HPCombi::BMat8, HPCombi::BMat8)> _act;
    std::function<std::bitset<256>(HPCombi::BMat8)>               _bit_func;
    bool                                                          _finished;
  };

  LIBSEMIGROUPS_TEST_CASE("BMat8 enum",
                          "001",
                          "test bliss canonicalisation",
                          "[quick][bliss]") {
    HPCombi::BMat8 x   = HPCombi::BMat8({{0, 0, 0}, {0, 0, 1}, {0, 1, 1}});
    HPCombi::BMat8 y   = HPCombi::BMat8({{0, 0, 0}, {0, 0, 1}, {1, 0, 1}});
    bliss_digraph  dgx = bliss_digraph_from_BMat8(x, 3);
    bliss_digraph  dgy = bliss_digraph_from_BMat8(y, 3);

    dgx.write_dot("dgx.dot");

    bliss_digraph dgx2 = bliss_digraph(6);
    dgx2.change_color(3, 1);
    dgx2.change_color(4, 1);
    dgx2.change_color(5, 1);
    dgx2.add_edge(1, 5);
    dgx2.add_edge(2, 4);
    dgx2.add_edge(2, 5);

    REQUIRE(dgx.cmp(dgx2) == 0);

    HPCombi::BMat8 canonx = permuted_BMat8(
        x, 3, dgx.canonical_form(stats, &bliss_hook_function, nullptr));
    HPCombi::BMat8 canony = permuted_BMat8(
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
    for (HPCombi::BMat8 i : enumerator.reps()) {
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
    REQUIRE(enumerator.reps().size() == 32);
    o.open("bmat_enum_trim_5.txt", std::ios::out | std::ios::trunc);
    for (HPCombi::BMat8 i : enumerator.reps()) {
      o << i.to_int() << "\n";
    }
    o.close();
  }

  LIBSEMIGROUPS_TEST_CASE("BMat8 enum",
                          "004",
                          "enumerate trim B6",
                          "[standard][enumerate]") {
    std::cout << "=> Finding minimal generating set for n = 6 ...\n";
    std::ofstream  o;
    auto           rg = ReportGuard(false);
    BMatEnumerator enumerator(6, true);
    REQUIRE(enumerator.reps().size() == 394);
    std::cout << "=> found " << enumerator.reps().size()
             << " trim matrices\n";

    o.open("build/output/bmat_enum_trim_6.txt", std::ios::out | std::ios::trunc);
    for (HPCombi::BMat8 i : enumerator.reps()) {
      o << i.to_int() << "\n";
    }
    o.close();
    std::cout << "=> writing file build/output/bmat_enum_trim_6.txt\n";

    std::cout << "=> filtering trim matrices ...\n";
    Filterer<6> f("build/output/bmat_enum_trim_6.txt",
                  "build/output/bmat_filtered_6.txt", {}, true);
    f.run();
    std::vector<HPCombi::BMat8> filtered = f.reps();
    REQUIRE(f.reps().size() == 64);
    std::ofstream ofs;
    ofs.open("build/output/bmat_filtered_6.txt",
             std::ofstream::out | std::ofstream::app);
    ofs << bmat8_helpers::elementary<HPCombi::BMat8>(6).to_int() << "\n";
    ofs << bmat8_helpers::one<HPCombi::BMat8>(5).to_int() << "\n";
    ofs << BMat8({{0, 1}, {1, 0}}).to_int() << "\n";
    ofs << BMat8({{0, 1, 0, 0, 0, 0},
                  {0, 0, 1, 0, 0, 0},
                  {0, 0, 0, 1, 0, 0},
                  {0, 0, 0, 0, 1, 0},
                  {0, 0, 0, 0, 0, 1},
                  {1, 0, 0, 0, 0, 0}})
               .to_int() << "\n";
    ofs.close();
    std::cout << "=> min. generating set has size " << f.reps().size() + 4
              << "\n";
    std::cout << "=> writing file build/output/bmat_filtered_6.txt\n";
  }

  LIBSEMIGROUPS_TEST_CASE("BMat8 enum",
                          "800",
                          "enumerate trim B8",
                          "[extreme][enumerate]") {
    std::ofstream  o;
    auto           rg = ReportGuard();
    BMatEnumerator enumerator_8_trim(8, true);
    o.open("bmat_trim_enum_8.txt", std::ios::out | std::ios::trunc);
    for (HPCombi::BMat8 i : enumerator_8_trim.reps()) {
      o << i.to_int() << "\n";
    }
    o.close();
    REQUIRE(enumerator_8_trim.reps().size() == 17120845);
  }

  LIBSEMIGROUPS_TEST_CASE("BMat8 enum", "006", "filter 5", "[quick]") {
    Filterer<5> f("bmat_enum_trim_5.txt", "bmat_filtered_5.txt", {}, true);
    f.run();
    std::vector<HPCombi::BMat8> filtered = f.reps();
    REQUIRE(filtered.size() == 9);
  }

  LIBSEMIGROUPS_TEST_CASE("BMat8 enum", "007", "filter 6", "[quick]") {
  }

  LIBSEMIGROUPS_TEST_CASE("BMat8 enum", "008", "filter 7", "[extreme]") {
    auto        rg = ReportGuard();
    Filterer<7> f("bmat_trim_enum_7.txt", "bmat_filtered_7.txt", {}, true);
    f.run();
    std::vector<HPCombi::BMat8> filtered = f.reps();
    REQUIRE(filtered.size() == 2139);
  }

  /////////////////////////////////////////////////////////////////////////////
  //
  // Filter the trim matrices in B8 by running tests 551 - 561
  // (after generating the trim B8 file!)
  // Or run the test executable with tag [filter8]
  // if things go wrong it might be because the tests aren't run in order
  // but you can force them to be with a flag to the executable (run with
  // --help)
  //
  /////////////////////////////////////////////////////////////////////////////

  LIBSEMIGROUPS_TEST_CASE("BMat8 enum",
                          "551",
                          "filter 8 by elementary containment - test Filterer",
                          "[extreme][filter8]") {
    HPCombi::BMat8 E({{1, 1, 0, 0, 0, 0, 0, 0},
                      {0, 1, 0, 0, 0, 0, 0, 0},
                      {0, 0, 1, 0, 0, 0, 0, 0},
                      {0, 0, 0, 1, 0, 0, 0, 0},
                      {0, 0, 0, 0, 1, 0, 0, 0},
                      {0, 0, 0, 0, 0, 1, 0, 0},
                      {0, 0, 0, 0, 0, 0, 1, 0},
                      {0, 0, 0, 0, 0, 0, 0, 1}});

    Filterer<8> filterer(
        "bmat_trim_enum_8.txt", "bmat_trim_8_filterered_1.txt", {E}, true);
    filterer.run();
    REQUIRE(filterer.size() == 16761162);
  }

  LIBSEMIGROUPS_TEST_CASE("BMat8 enum",
                          "552",
                          "filter 8 by extended 7 gens - by Filterer",
                          "[extreme][filter8]") {
    std::vector<HPCombi::BMat8> bmat7_gens(0);
    std::ifstream               f("bmat_filtered_7.txt");
    std::string                 line;
    while (std::getline(f, line)) {
      bmat7_gens.push_back(HPCombi::BMat8(std::stoul(line) | 1));
    }
    f.close();

    HPCombi::BMat8 F({{1, 0, 0, 0, 0, 0, 0, 0},
                      {0, 1, 0, 0, 0, 0, 0, 0},
                      {0, 0, 1, 0, 0, 0, 0, 0},
                      {0, 0, 0, 1, 0, 0, 0, 0},
                      {0, 0, 0, 0, 1, 0, 0, 0},
                      {0, 0, 0, 0, 0, 1, 0, 0},
                      {0, 0, 0, 0, 0, 0, 1, 0},
                      {0, 0, 0, 0, 0, 0, 0, 0}});

    bmat7_gens.push_back(F);

    Filterer<8> filterer("bmat_trim_8_filterered_1.txt",
                         "bmat_trim_8_filterered_2.txt",
                         bmat7_gens,
                         false);
    filterer.run();
    REQUIRE(filterer.size() == 15336159);
  }

  LIBSEMIGROUPS_TEST_CASE("BMat8 enum",
                          "553",
                          "filter 8 by special extended 7 gens - by Filterer",
                          "[extreme][filter8]") {
    std::vector<HPCombi::BMat8> bmat7_gens(0);
    std::ifstream               f("bmat_filtered_7.txt");
    std::string                 line;
    while (std::getline(f, line)) {
      HPCombi::BMat8 bm(HPCombi::BMat8(std::stoul(line) | 1));
      if (row_space_bitset(bm).count() == 160) {
        bmat7_gens.push_back(bm);
      }
    }
    f.close();

    Filterer<8> filterer("bmat_trim_8_filterered_2.txt",
                         "bmat_trim_8_filterered_3.txt",
                         bmat7_gens,
                         true);
    filterer.run();
    REQUIRE(filterer.size() == 2925520);
  }

  LIBSEMIGROUPS_TEST_CASE("BMat8 enum",
                          "554",
                          "filter 8 by special extended 7 gens - by Filterer",
                          "[extreme][filter8]") {
    std::vector<HPCombi::BMat8> bmat7_gens(0);
    std::ifstream               f("bmat_filtered_7.txt");
    std::string                 line;
    while (std::getline(f, line)) {
      HPCombi::BMat8 bm(HPCombi::BMat8(std::stoul(line) | 1));
      size_t         count = row_space_bitset(bm).count();
      if (count == 160 || count == 136 || count == 144) {
        bmat7_gens.push_back(bm);
      }
    }
    f.close();

    Filterer<8> filterer("bmat_trim_8_filterered_3.txt",
                         "bmat_trim_8_filterered_4.txt",
                         bmat7_gens,
                         true);
    filterer.run();
    REQUIRE(filterer.size() == 2066807);
  }

  LIBSEMIGROUPS_TEST_CASE("BMat8 enum",
                          "555",
                          "filter 8 by special extended 7 gens - by Filterer",
                          "[extreme][filter8]") {
    std::vector<HPCombi::BMat8> bmat7_gens(0);
    std::ifstream               f("bmat_filtered_7.txt");
    std::string                 line;
    while (std::getline(f, line)) {
      HPCombi::BMat8 bm(HPCombi::BMat8(std::stoul(line) | 1));
      size_t         count = row_space_bitset(bm).count();
      if (count == 112 || count == 116 || count == 120) {
        bmat7_gens.push_back(bm);
      }
    }
    f.close();

    Filterer<8> filterer("bmat_trim_8_filterered_4.txt",
                         "bmat_trim_8_filterered_5.txt",
                         bmat7_gens,
                         true);
    filterer.run();
    REQUIRE(filterer.size() == 1437821);
  }

  LIBSEMIGROUPS_TEST_CASE("BMat8 enum",
                          "556",
                          "filter 8 by special extended 7 gens - by Filterer",
                          "[extreme][filter8]") {
    std::vector<HPCombi::BMat8> bmat7_gens(0);
    std::ifstream               f("bmat_filtered_7.txt");
    std::string                 line;
    while (std::getline(f, line)) {
      HPCombi::BMat8 bm(HPCombi::BMat8(std::stoul(line) | 1));
      size_t         count = row_space_bitset(bm).count();
      if (count == 104 || count == 108) {
        bmat7_gens.push_back(bm);
      }
    }
    f.close();

    Filterer<8> filterer("bmat_trim_8_filterered_5.txt",
                         "bmat_trim_8_filterered_6.txt",
                         bmat7_gens,
                         true);
    filterer.run();
    REQUIRE(filterer.size() == 1121872);
  }

  LIBSEMIGROUPS_TEST_CASE("BMat8 enum",
                          "557",
                          "filter 8 by special extended 7 gens - by Filterer",
                          "[extreme][filter8]") {
    std::vector<HPCombi::BMat8> bmat7_gens(0);
    std::ifstream               f("bmat_filtered_7.txt");
    std::string                 line;
    while (std::getline(f, line)) {
      HPCombi::BMat8 bm(HPCombi::BMat8(std::stoul(line) | 1));
      size_t         count = row_space_bitset(bm).count();
      if (count == 100 || count == 102) {
        bmat7_gens.push_back(bm);
      }
    }
    f.close();

    Filterer<8> filterer("bmat_trim_8_filterered_6.txt",
                         "bmat_trim_8_filterered_7.txt",
                         bmat7_gens,
                         true);
    filterer.run();
    REQUIRE(filterer.size() == 968568);
  }

  LIBSEMIGROUPS_TEST_CASE("BMat8 enum",
                          "558",
                          "filter 8 by special extended 7 gens - by Filterer",
                          "[extreme][filter8]") {
    std::vector<HPCombi::BMat8> bmat7_gens(0);
    std::ifstream               f("bmat_filtered_7.txt");
    std::string                 line;
    while (std::getline(f, line)) {
      HPCombi::BMat8 bm(HPCombi::BMat8(std::stoul(line) | 1));
      size_t         count = row_space_bitset(bm).count();
      if (count == 96) {
        bmat7_gens.push_back(bm);
      }
    }
    f.close();

    Filterer<8> filterer("bmat_trim_8_filterered_7.txt",
                         "bmat_trim_8_filterered_8.txt",
                         bmat7_gens,
                         true);
    filterer.run();
    REQUIRE(filterer.size() == 900863);
  }

  LIBSEMIGROUPS_TEST_CASE("BMat8 enum",
                          "559",
                          "filter 8 by extended 7 gens columns - by Filterer",
                          "[extreme][filter8]") {
    std::vector<HPCombi::BMat8> bmat7_gens(0);
    std::ifstream               f("bmat_filtered_7.txt");
    std::string                 line;
    while (std::getline(f, line)) {
      HPCombi::BMat8 bm(HPCombi::BMat8(std::stoul(line) | 1));
      bmat7_gens.push_back(bm);
    }
    f.close();

    Filterer<8> filterer("bmat_trim_8_filterered_8.txt",
                         "bmat_trim_8_filterered_9.txt",
                         bmat7_gens,
                         false,
                         true);
    filterer.run();
    REQUIRE(filterer.size() == 892502);
  }

  LIBSEMIGROUPS_TEST_CASE("BMat8 enum",
                          "560",
                          "filter 8 by extended 7 gens columns - by Filterer",
                          "[extreme][filter8]") {
    std::vector<HPCombi::BMat8> bmat7_gens(0);
    std::ifstream               f("bmat_filtered_7.txt");
    std::string                 line;
    while (std::getline(f, line)) {
      HPCombi::BMat8 bm(HPCombi::BMat8(std::stoul(line) | 1));
      size_t         count = col_space_bitset(bm).count();
      if (count == 144 || count == 160) {
        bmat7_gens.push_back(bm);
      }
    }
    f.close();

    Filterer<8> filterer("bmat_trim_8_filterered_9.txt",
                         "bmat_trim_8_filterered_10.txt",
                         bmat7_gens,
                         true,
                         true);
    filterer.run();
    REQUIRE(filterer.size() == 570132);
  }

  LIBSEMIGROUPS_TEST_CASE("BMat8 enum",
                          "561",
                          "filter 8 by extended 7 gens columns - by Filterer",
                          "[extreme][filter8]") {
    std::vector<HPCombi::BMat8> bmat7_gens(0);
    std::ifstream               f("bmat_filtered_7.txt");
    std::string                 line;
    while (std::getline(f, line)) {
      HPCombi::BMat8 bm(HPCombi::BMat8(std::stoul(line) | 1));
      size_t         count = col_space_bitset(bm).count();
      if (count == 128 || count == 136) {
        bmat7_gens.push_back(bm);
      }
    }
    f.close();

    Filterer<8> filterer("bmat_trim_8_filterered_10.txt",
                         "bmat_trim_8_filterered_11.txt",
                         bmat7_gens,
                         true,
                         true);
    filterer.run();
    REQUIRE(filterer.size() == 523556);
  }

  LIBSEMIGROUPS_TEST_CASE("BMat8 enum",
                          "562",
                          "filter 8 by simple 7 gens columns - by Filterer",
                          "[extreme][filter8]") {
    std::vector<HPCombi::BMat8> bmat7_gens(0);
    std::ifstream               f("bmat_filtered_7.txt");
    std::string                 line;

    while (std::getline(f, line)) {
      HPCombi::BMat8 bm(HPCombi::BMat8(std::stoul(line)));
      size_t         count = col_space_bitset(bm).count();
      if (count >= 60) {
        bmat7_gens.push_back(bm);
      }
    }
    f.close();

    std::cout << "processing with " << bmat7_gens.size();

    Filterer<8> filterer("bmat_trim_8_filterered_11.txt",
                         "bmat_trim_8_filterered_12.txt",
                         bmat7_gens,
                         true,
                         true);
    filterer.run();
    REQUIRE(filterer.size() == 472207);
  }

  ///////////////////////////////////////////////////////////////////////
  // END FILTERING B8
  ///////////////////////////////////////////////////////////////////////

  ///////////////////////////////////////////////////////////////////////
  // OUTPUT ROW SPACES
  ///////////////////////////////////////////////////////////////////////

  LIBSEMIGROUPS_TEST_CASE("BMat8 enum",
                          "571",
                          "bmat8 row space digraphs 8",
                          "[standard]") {
    std::vector<HPCombi::BMat8> bmat_enum;
    std::ifstream               f("bmat_full_filtered_8.txt");
    std::string                 line;
    while (std::getline(f, line)) {
      bmat_enum.push_back(HPCombi::BMat8(std::stoul(line)));
    }
    f.close();

    std::ofstream o;
    o.open("row_space_numbers_8.txt", std::ios::out | std::ios::trunc);
    for (auto it = bmat_enum.cbegin(); it < bmat_enum.cend(); it++) {
      std::vector<uint8_t> row_space = row_space_vector(*it);
      o << "[";
      for (size_t i = 0; i < row_space.size() - 1; ++i) {
        o << row_space[i] << ", ";
      }
      o << row_space.back() << "]\n";
    }
    o.close();
  }

  ////////////////////////////////////////////////////////////////////

  LIBSEMIGROUPS_TEST_CASE("BMat8 enum",
                          "129",
                          "bmat7 gens row space sizes",
                          "[extreme]") {
    auto                        rg = ReportGuard();
    std::vector<HPCombi::BMat8> bmat_enum;
    std::ifstream               f("bmat_gens_test_6.txt");
    std::string                 line;
    while (std::getline(f, line)) {
      bmat_enum.push_back(HPCombi::BMat8(std::stoul(line)));
    }
    f.close();
    for (HPCombi::BMat8 x : bmat_enum) {
      std::cout << x.row_space_size() << std::endl;
    }
  }

  LIBSEMIGROUPS_TEST_CASE("BMat8 enum",
                          "130",
                          "bmat7 gens row space sizes",
                          "[extreme]") {
    auto                        rg = ReportGuard();
    std::vector<HPCombi::BMat8> bmat_enum;
    std::ifstream               f("bmat_gens_7.txt");
    std::string                 line;
    while (std::getline(f, line)) {
      bmat_enum.push_back(HPCombi::BMat8(std::stoul(line)));
    }
    f.close();
    std::vector<size_t> row_space_count(129, 0);
    for (HPCombi::BMat8 x : bmat_enum) {
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
    std::vector<HPCombi::BMat8> bmat_enum;
    std::ifstream               f("bmat_trim_enum_6.txt");
    std::string                 line;
    while (std::getline(f, line)) {
      bmat_enum.push_back(HPCombi::BMat8(std::stoul(line)));
    }
    f.close();

    std::ofstream o;
    o.open("row_spaces_6.txt", std::ios::out | std::ios::trunc);
    for (HPCombi::BMat8 bm : bmat_enum) {
      std::vector<uint8_t> row_space = row_space_vector(bm);
      o << "{";
      for (uint8_t x : row_space) {
        o << +x << ", ";
      }
      o << "0}" << std::endl;
      ;
    }
    o.close();
  }

  LIBSEMIGROUPS_TEST_CASE("BMat8 enum",
                          "132",
                          "bmat8 row space digraphs 5",
                          "[quick]") {
    std::vector<HPCombi::BMat8> bmat_enum;
    std::ifstream               f("bmat_trim_enum_5.txt");
    std::string                 line;
    while (std::getline(f, line)) {
      bmat_enum.push_back(HPCombi::BMat8(std::stoul(line)));
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
    std::vector<HPCombi::BMat8> bmat_enum;
    std::ifstream               f("bmat_trim_enum_6.txt");
    std::string                 line;
    while (std::getline(f, line)) {
      bmat_enum.push_back(HPCombi::BMat8(std::stoul(line)));
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
    std::vector<HPCombi::BMat8> bmat_enum;
    std::ifstream               f("bmat_trim_enum_7.txt");
    std::string                 line;
    while (std::getline(f, line)) {
      bmat_enum.push_back(HPCombi::BMat8(std::stoul(line)));
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
    std::vector<HPCombi::BMat8> bmat_enum;
    std::ifstream               f("bmat_gens_test_5.txt");
    std::string                 line;
    while (std::getline(f, line)) {
      bmat_enum.push_back(HPCombi::BMat8(std::stoul(line)));
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
    std::vector<HPCombi::BMat8> bmat_enum;
    std::ifstream               f("bmat_gens_test_6.txt");
    std::string                 line;
    while (std::getline(f, line)) {
      bmat_enum.push_back(HPCombi::BMat8(std::stoul(line)));
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
    std::vector<HPCombi::BMat8> bmat_enum;
    std::ifstream               f("bmat_gens_test_7.txt");
    std::string                 line;
    while (std::getline(f, line)) {
      bmat_enum.push_back(HPCombi::BMat8(std::stoul(line)));
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
    std::cout << HPCombi::BMat8(13853107707017691136ull) << std::endl;
    std::cout << HPCombi::BMat8(4620710844303409152) << std::endl;
    std::cout << HPCombi::BMat8(4647750068672397312) << std::endl;
    std::cout << HPCombi::BMat8(9241421688590303232ull) << std::endl;

    std::cout << bmat8_helpers::elementary(6).to_int() << std::endl;
    std::cout << bmat8_helpers::one(5).to_int() << std::endl;
  }

  LIBSEMIGROUPS_TEST_CASE("BMat8 enum", "170", "prime extensions", "[quick]") {
    std::vector<HPCombi::BMat8> bmat_enum;
    std::ifstream               f("bmat_gens_7.txt");
    std::string                 line;
    while (std::getline(f, line)) {
      bmat_enum.push_back(HPCombi::BMat8(std::stoul(line)));
    }
    f.close();
    HPCombi::BMat8 x = bmat_enum[1021];
    for (HPCombi::BMat8& bm : simple_prime_extensions(x, 7)) {
      std::cout << bm << std::endl;
    }
  }

  LIBSEMIGROUPS_TEST_CASE("BMat8 enum",
                          "180",
                          "produce extended B7 primes",
                          "[quick]") {
    std::vector<HPCombi::BMat8> bmat7_gens(0);
    std::vector<HPCombi::BMat8> bmat8_enum(0);
    std::ifstream               f("bmat_gens_7.txt");
    std::string                 line;
    std::ofstream               o;
    o.open("bmat8_filterers.txt", std::ios::out | std::ios::trunc);
    while (std::getline(f, line)) {
      o << HPCombi::BMat8(std::stoul(line) | 1).to_int() << "\n";
    }
    f.close();
    o.close();
  }

  LIBSEMIGROUPS_TEST_CASE("BMat8 enum",
                          "190",
                          "bmat7 gens extended row space sizes",
                          "[extreme]") {
    auto                        rg = ReportGuard();
    std::vector<HPCombi::BMat8> bmat_enum;
    std::ifstream               f("bmat_gens_7.txt");
    std::string                 line;
    while (std::getline(f, line)) {
      bmat_enum.push_back(HPCombi::BMat8(std::stoul(line)));
    }
    f.close();
    std::vector<size_t> row_space_count(257, 0);
    for (HPCombi::BMat8 x : bmat_enum) {
      row_space_count[x.row_space_size()]++;
    }
    for (size_t i = 0; i < 257; ++i) {
      std::cout << "row spaces of size " << i << ": " << row_space_count[i]
                << std::endl;
    }
  }

  LIBSEMIGROUPS_TEST_CASE("BMat8 enum",
                          "300",
                          "count fully indecomposable",
                          "[quick]") {
    std::vector<HPCombi::BMat8> bmat8_enum;
    std::ifstream               f("bmat_trim_enum_8.txt");
    std::string                 line;
    while (std::getline(f, line)) {
      bmat8_enum.push_back(HPCombi::BMat8(std::stoul(line)));
    }
    f.close();
    REQUIRE(bmat8_enum.size() == 17120845);

    std::cout << "finished reading!" << std::endl;

    std::ofstream o;
    o.open("fully_indecomposable_8.txt", std::ios::out | std::ios::trunc);
    size_t count = 0;
    for (HPCombi::BMat8 bm : bmat8_enum) {
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
    std::vector<HPCombi::BMat8> bmat8_enum;
    std::ifstream               f("bmat_filtered_8_17.txt");
    std::string                 line;
    while (std::getline(f, line)) {
      bmat8_enum.push_back(HPCombi::BMat8(std::stoul(line)));
    }
    f.close();
    REQUIRE(bmat8_enum.size() == 471912);

    std::cout << "finished reading!" << std::endl;
    size_t count = 0;
    for (HPCombi::BMat8 bm : bmat8_enum) {
      count = is_fully_indecomposable(bm, 8) ? count + 1 : count;
    }
    std::cout << +count << " filtered fully indecomposable trim reps of dim 8"
              << std::endl;
  }
  } // namespace libsemigroups
