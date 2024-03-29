#include <bitset>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <set>            // for set
#include <unordered_set>  // for unordered_set
#include <vector>         // for vector

#include <algorithm>
#include <random>

#include "catch.hpp"
#include "test-main.hpp"

#include "bmat8.hpp"
#include "element-helper.hpp"
#include "element.hpp"
#include "froidure-pin.hpp"
#include "hpcombi.hpp"
#include "libsemigroups-exception.hpp"  // for LIBSEMIGROUPS_EXCEPTION
#include "runner.hpp"
#include "string.hpp"
#include "timer.hpp"

#include "../extern/bliss-0.73/graph.hh"

#include <omp.h>
namespace libsemigroups {
  typedef bliss_digraphs::Digraph bliss_digraph;

  bliss_digraphs::Stats stats;

  bliss_digraph &bliss_digraph_from_BMat8(HPCombi::BMat8 bm, size_t dim = 8) {
    static bliss_digraph out = bliss_digraph(2 * dim);
    out.clear();
    size_t x = bm.to_int();
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

  bliss_digraph &bliss_locked_digraph_from_BMat8(HPCombi::BMat8 bm,
                                                 size_t         dim = 8) {
    static bliss_digraph out = bliss_digraph(3 * dim);
    out.clear();
    size_t x = bm.to_int();
    for (size_t i = 0; i < dim; ++i) {
      out.change_color(i, 0);
      out.change_color(dim + i, 1);
      out.change_color(2 * dim + i, 2);
      out.add_edge(i, 2 * dim + i);
      out.add_edge(2 * dim + i, i);
      out.add_edge(dim + i, 2 * dim + i);
      out.add_edge(2 * dim + i, dim + i);
      for (size_t j = 0; j < dim; ++j) {
        if ((x >> (63 - 8 * i - j)) & 1) {
          out.add_edge(i, dim + j);
        }
      }
    }
    out.write_dot("../output/aaaaa.dot");
    return out;
  }

  void bliss_hook_function(void *, const unsigned int, const unsigned int *) {}

  HPCombi::BMat8
  permuted_BMat8(HPCombi::BMat8                                   bm,
                 size_t                                           dim,
                 bliss_digraphs::uint_pointer_to_const_substitute perm) {
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

  HPCombi::BMat8
  locked_permuted_BMat8(HPCombi::BMat8                                   bm,
                        size_t                                           dim,
                        bliss_digraphs::uint_pointer_to_const_substitute perm) {
    HPCombi::epu8 row_perm;
    for (size_t i = 0; i < dim; ++i) {
      row_perm[perm[i]] = i;
    }
    for (size_t i = dim; i < 8; ++i) {
      row_perm[i] = i;
    }

    return bm.row_permuted(row_perm).col_permuted(row_perm);
  }

  HPCombi::BMat8 canonical_BMat8(HPCombi::BMat8 bm, size_t dim) {
    bliss_digraph &dg = bliss_digraph_from_BMat8(bm, dim);
    return permuted_BMat8(
        bm, dim, dg.canonical_form(stats, &bliss_hook_function, nullptr));
  }

  HPCombi::BMat8 canonical_locked_BMat8(HPCombi::BMat8 bm, size_t dim) {
    bliss_digraph &dg = bliss_locked_digraph_from_BMat8(bm, dim);
    return locked_permuted_BMat8(
        bm, dim, dg.canonical_form(stats, &bliss_hook_function, nullptr));
  }

  bool is_row_reduced(HPCombi::BMat8 bm) {
    return bm.nr_rows() == bm.row_space_basis().nr_rows();
  }

  bool is_col_reduced(HPCombi::BMat8 bm) {
    return is_row_reduced(bm.transpose());
  }

  std::pair<std::bitset<256>, std::vector<uint8_t>>
  row_space_impl(HPCombi::BMat8 bm) {
    std::bitset<256>     lookup;
    std::vector<uint8_t> row_vec = bm.row_space_basis().rows();
    auto                 last = std::remove(row_vec.begin(), row_vec.end(), 0);
    row_vec.erase(last, row_vec.end());
    row_vec.push_back(0);  // Maybe push front?
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
    return std::make_pair(lookup, row_space);
  }

  std::bitset<256> row_space_bitset(HPCombi::BMat8 bm) {
    return row_space_impl(bm).first;
  }

  std::vector<uint8_t> row_space_vector(HPCombi::BMat8 bm) {
    return row_space_impl(bm).second;
  }
  
  size_t row_space_size(HPCombi::BMat8 bm) {
    return row_space_vector(bm).size();
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

  int nr_ones(uint8_t x) {
    return __builtin_popcount(x);
  }

  // not compiler safe!
  int nr_ones(HPCombi::BMat8 x) {
    return __builtin_popcountl(x.to_int());
  }

  bool is_row_trim(HPCombi::BMat8 bm, size_t dim = 8) {
    static std::array<uint8_t, 8> r;
    r.fill(0);
    r = bm.row_array();
    for (size_t i = 0; i < dim; ++i) {
      for (size_t j = 0; j < dim; ++j) {
        if (r[i] && (i != j) && ((r[i] | r[j]) == r[j])) {
          return false;
        }
      }
    }
    return true;
  }

  bool is_col_trim(HPCombi::BMat8 bm, size_t dim = 8) {
    return is_row_trim(bm.transpose(), dim);
  }

  // given a HPCombi::BMat8 bm of dimension dim, find all the BMat8s of
  // dimension one higher which are extensions of bm with 0s in the new column
  // and maximal row spaces
  // TODO: obvious optimisations
  std::vector<HPCombi::BMat8> simple_prime_extensions(HPCombi::BMat8 bm,
                                                      size_t const   dim) {
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

  HPCombi::BMat8 reflect_vertically(HPCombi::BMat8 bm, size_t const dim) {
    static std::array<uint8_t, 256> const uint8_reflects
        = {0,   128, 64,  192, 32,  160, 96,  224, 16,  144, 80,  208, 48,  176,
           112, 240, 8,   136, 72,  200, 40,  168, 104, 232, 24,  152, 88,  216,
           56,  184, 120, 248, 4,   132, 68,  196, 36,  164, 100, 228, 20,  148,
           84,  212, 52,  180, 116, 244, 12,  140, 76,  204, 44,  172, 108, 236,
           28,  156, 92,  220, 60,  188, 124, 252, 2,   130, 66,  194, 34,  162,
           98,  226, 18,  146, 82,  210, 50,  178, 114, 242, 10,  138, 74,  202,
           42,  170, 106, 234, 26,  154, 90,  218, 58,  186, 122, 250, 6,   134,
           70,  198, 38,  166, 102, 230, 22,  150, 86,  214, 54,  182, 118, 246,
           14,  142, 78,  206, 46,  174, 110, 238, 30,  158, 94,  222, 62,  190,
           126, 254, 1,   129, 65,  193, 33,  161, 97,  225, 17,  145, 81,  209,
           49,  177, 113, 241, 9,   137, 73,  201, 41,  169, 105, 233, 25,  153,
           89,  217, 57,  185, 121, 249, 5,   133, 69,  197, 37,  165, 101, 229,
           21,  149, 85,  213, 53,  181, 117, 245, 13,  141, 77,  205, 45,  173,
           109, 237, 29,  157, 93,  221, 61,  189, 125, 253, 3,   131, 67,  195,
           35,  163, 99,  227, 19,  147, 83,  211, 51,  179, 115, 243, 11,  139,
           75,  203, 43,  171, 107, 235, 27,  155, 91,  219, 59,  187, 123, 251,
           7,   135, 71,  199, 39,  167, 103, 231, 23,  151, 87,  215, 55,  183,
           119, 247, 15,  143, 79,  207, 47,  175, 111, 239, 31,  159, 95,  223,
           63,  191, 127, 255};
    std::vector<uint8_t> rows = bm.rows();
    for (auto it = rows.begin(); it < rows.end(); ++it) {
      *it = uint8_reflects[*it] << (8 - dim);
    }
    return bmat8_helpers::make<HPCombi::BMat8>(rows.cbegin(), rows.cend());
  }

  bool cmp_by_fewer_ones(HPCombi::BMat8 x, HPCombi::BMat8 y) {
    return nr_ones(x) > nr_ones(y);
  }

  ////////////////////////////////////////////////////////////////////////
  // JDM helper functions
  ////////////////////////////////////////////////////////////////////////

  std::vector<HPCombi::BMat8> read_bmat_file(std::string const &fn) {
    std::vector<HPCombi::BMat8> out;
    std::ifstream               f(fn);
    std::string                 line;
    while (std::getline(f, line)) {
      out.push_back(HPCombi::BMat8(std::stoul(line)));
    }
    f.close();
    return out;
  }

  void write_bmat_file(std::string const &                fn,
                       std::vector<HPCombi::BMat8> const &vec) {
    std::ofstream f(fn, std::ios::out | std::ios::trunc);
    for (auto x : vec) {
      f << x.to_int() << "\n";
    }
    f.close();
  }

  void append_bmat_file(std::string const &fn, HPCombi::BMat8 x) {
    std::ofstream f(fn, std::ios::out | std::ios::app);
    f << x.to_int() << "\n";
    f.close();
  }

  ////////////////////////////////////////////////////////////////////////
  // Heavy lifters
  ////////////////////////////////////////////////////////////////////////

  class BMatEnumerator : public libsemigroups::Runner {
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
            _rows[k] = row;
            auto bm  = bmat8_helpers::make<HPCombi::BMat8>(_rows.cbegin(),
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
      for (bool &x : _row_seen) {
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

    std::vector<HPCombi::BMat8> const &reps() {
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

  class BMatReflexiveTrimStraightEnumerator : public libsemigroups::Runner {
   public:
    BMatReflexiveTrimStraightEnumerator(size_t dim)
        : _n(dim),
          _max((1 << _n) - 1),
          _rows(8, 0),
          _out(),
          _set(),
          _finished(false) {}

   private:
    void dive(size_t k) {
      for (uint8_t row = 1 << k; row < _max; ++row) {
        _rows[k] = row;
        if (!(row & (1 << k))) {
          goto next_loop;
        }
        for (size_t i = 0; i < k; ++i) {
          if ((_rows[i] | row) == row || ((_rows[i] | row) == _rows[i])) {
            goto next_loop;
          }
        }
        if (k == _n - 1) {
          auto bm = bmat8_helpers::make<HPCombi::BMat8>(_rows.cbegin(),
                                                        _rows.cend());
          // move the matrix to the right place
          bm = HPCombi::BMat8(bm.to_int() << (8 - _n));

          // flip the matrix so it's actually reflexive
          bm = reflect_vertically(bm, _n);
          if (is_col_trim(bm, _n)) {
            HPCombi::BMat8 canon = canonical_locked_BMat8(bm, _n);
            if (_set.find(canon) == _set.end()) {
              _set.insert(canon);
              _out.push_back(canon);
            }
          }
          if (report()) {
            REPORT_DEFAULT("found %d reps so far, currently on %s \n",
                           _out.size(),
                           detail::to_string(bm));
          }
        } else {
          // DIVE
          dive(k + 1);
        }
      next_loop:;
      }
    }

   public:
    void run_impl() {
      dive(0);
      HPCombi::BMat8 one = bmat8_helpers::one<HPCombi::BMat8>(_n);
      if (!_set.count(one)) {
        _out.push_back(one);
      }
      _finished = true;
      report_why_we_stopped();
    }

    bool finished_impl() const {
      return _finished;
    }

    std::vector<HPCombi::BMat8> const &reps() {
      if (!started()) {
        run();
      }
      return _out;
    }

   private:
    size_t                             _n;
    size_t                             _max;
    std::vector<uint8_t>               _rows;
    std::vector<HPCombi::BMat8>        _out;
    std::unordered_set<HPCombi::BMat8> _set;
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
          _finished(false),
          _mtx() {
      if (col_filt) {
        _act      = left_mult;
        _bit_func = col_space_bitset;
      } else {
        _act      = right_mult;
        _bit_func = row_space_bitset;
      }
    }
    void run_impl() {
      auto                        rg   = ReportGuard();
      std::vector<HPCombi::BMat8> bmat_enum;
      std::ifstream               f(_in);
      std::string                 line;
      while (std::getline(f, line)) {
        bmat_enum.push_back(HPCombi::BMat8(std::stoul(line)));
      }
      f.close();
      REPORT_DEFAULT("read " + std::to_string(bmat_enum.size())
                     + " matrices to filter from " + _in + "\n");
      REPORT_DEFAULT("filtering by " + std::to_string(_filterers.size())
                     + " matrices" + (_permute ? " with permutation" : "")
                     + "\n");
      using Perm = typename PermHelper<_dim>::type;
      HPCombi::epu8 cycle;
      HPCombi::epu8 transposition;
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
        S_bmats.push_back(bmat8_helpers::make<_dim,
                                              typename PermHelper<_dim>::type,
                                              HPCombi::BMat8>(p));
      }

      if (_filterers.size() == 0) {
        _filterers = bmat_enum;
        _filterers.push_back(bmat8_helpers::elementary<HPCombi::BMat8>(_dim));
      }

      std::vector<std::vector<std::bitset<256>>> spaces(
          (1 << _dim) + 1, std::vector<std::bitset<256>>());
      
      std::atomic<size_t> count{0};
#pragma omp parallel for
      for (size_t i = 0; i < _filterers.size(); ++i) {
        HPCombi::BMat8 x = _filterers[i];
        std::vector<std::bitset<256>> thread_bitsets;
        std::bitset<256> bitset = _bit_func(x);
        size_t           bit_count  = bitset.count();
        thread_bitsets.push_back(bitset);
        if (_permute) {
          for (HPCombi::BMat8 y : S_bmats) {
            bitset = _bit_func(_act(x, y));
            thread_bitsets.push_back(bitset);
          }
        }
        {
          std::lock_guard<std::mutex> lg(_mtx);
          spaces[bit_count].insert(spaces[bit_count].end(),
                                   thread_bitsets.begin(),
                                   thread_bitsets.end());
        }
        count++;
        if (report()) {
          REPORT_DEFAULT("Generating bitsets: on %d out of %d.\n",
                         count.load(),
                         _filterers.size());
        }
      }

      count = 0;
#pragma omp parallel for
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
          std::lock_guard<std::mutex> lg(_mtx);
          _filtered.push_back(bm);
        }
        count++;
        if (report()) {
          REPORT_DEFAULT("On %d (overall %d out of %d), keeping %d.\n",
                         i + 1,
                         count.load(),
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

    const std::vector<HPCombi::BMat8> &reps() {
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
    std::mutex                  _mtx;
  };

  std::vector<std::vector<uint8_t>> intersections(HPCombi::BMat8 bm,
                                                  size_t         dim) {
    std::vector<std::vector<uint8_t>> out(dim, std::vector<uint8_t>(0));
    std::vector<uint8_t>              rows = bm.rows();
    std::vector<bool>                 lookup(256, false);
    for (size_t i = 1; i < (1 << dim); ++i) {
      size_t  choice = i << (8 - dim);
      uint8_t res    = -1;
      for (size_t j = 0; (j < dim) && res; ++j) {
        if ((128 >> j) & choice) {
          res &= rows[j];
        }
      }
      if (!lookup[res]) {
        lookup[res] = true;
        for (size_t j = 0; j < dim; ++j) {
          if ((128 >> j) & res) {
            out[j].push_back(res);
          }
        }
      }
    }
    return out;
  }

  // Returns a candidate left multiplier bm in the reflexive monoid such that
  // bm * x = y
  // NOTE: Only a candidate! (but if bm * x != y, no such bm actually
  // exists)
  HPCombi::BMat8 left_reflexive_multiplier(HPCombi::BMat8 x,
                                           HPCombi::BMat8 y,
                                           size_t         dim) {
    HPCombi::BMat8       out(0);
    std::vector<uint8_t> x_rows = x.rows();
    std::vector<uint8_t> y_rows = y.rows();
    for (size_t i = 0; i < dim; ++i) {
      for (size_t j = 0; j < dim; ++j) {
        if ((x_rows[j] | y_rows[i]) == y_rows[i]) {
          out.set(i, j, true);
        }
      }
      out.set(i, i, true);
    }
    return out;
  }

  bool increment_tuple(std::vector<uint8_t> &    tup,
                       const std::vector<size_t> maxes) {
    for (auto it = tup.rbegin(); it < tup.rend(); ++it) {
      if (*it < maxes[maxes.size() - 1 - std::distance(tup.rbegin(), it)]) {
        (*it)++;
        return true;
      } else {
        *it = 0;
      }
    }
    return false;
  }

  bool is_maximal_reflexive_bmat(HPCombi::BMat8 x, size_t dim) {
    std::vector<uint8_t>              tup(dim, 0);
    std::vector<std::vector<uint8_t>> intersects = intersections(x, dim);
    std::vector<size_t>               maxes(0);
    for (auto x : intersects) {
      maxes.push_back(x.size() - 1);
    }

    do {
      size_t o = 0;
      for (size_t i = 0; i < dim; ++i) {
        o |= (static_cast<size_t>(intersects[i][tup[i]]) << ((7 - i) * 8));
      }
      HPCombi::BMat8 bm(o);
      HPCombi::BMat8 mult = left_reflexive_multiplier(bm, x, dim);
      if (mult != x && bm != x && mult * bm == x) {
        return false;
      }
    } while (increment_tuple(tup, maxes));
    return true;
  }

  bool in_left_reflexive_ideal(HPCombi::BMat8 x, HPCombi::BMat8 y, size_t dim) {
    return left_reflexive_multiplier(x, y, dim) * x == y;
  }

  template <size_t _dim>
  class ReflexiveFilterer : public ::libsemigroups::Runner {
   public:
    ReflexiveFilterer(std::string in, std::string out, std::string disc)
        : _in(in),
          _out(out),
          _discarded(disc),
          _filtered(0),
          _finished(false),
          _mtx(),
          _remaining() {
      std::ifstream f(_in);
      std::string   line;
      while (std::getline(f, line)) {
        _remaining.insert(HPCombi::BMat8(std::stoul(line)));
      }
      f.close();
    }
    void run_impl() {
      auto rg = ReportGuard();

      std::vector<HPCombi::BMat8> bmat_enum(_remaining.begin(),
                                            _remaining.end());
      REPORT_DEFAULT("read " + std::to_string(bmat_enum.size())
                     + " matrices to filter from " + _in + "\n");
      HPCombi::BMat8 one = bmat8_helpers::one<HPCombi::BMat8>(_dim);

      // clear the files
      write_bmat_file(_out, _filtered);
      write_bmat_file(_discarded, _filtered);


      // shuffle so the hardest matrices aren't clumped together
      std::random_device rd;
      std::mt19937       g(rd());
      std::shuffle(bmat_enum.begin(), bmat_enum.end(), g);

      std::atomic<size_t> count{0};
      bool break_flag = false;
#pragma omp parallel for schedule(dynamic, 1)
      for (size_t i = 0; i < bmat_enum.size(); ++i) {
        // can't break in omp loops
        if (break_flag) {
          continue;
        }
        HPCombi::BMat8 bm = bmat_enum[i];

        std::vector<uint8_t>              tup(_dim, 0);
        std::vector<std::vector<uint8_t>> intersects = intersections(bm, _dim);
        std::vector<size_t>               maxes(0);
        bool                              is_maximal = true;
        for (auto x : intersects) {
          maxes.push_back(x.size() - 1);
        }
        do {
          size_t o = 0;
          for (size_t i = 0; i < _dim; ++i) {
            o |= (static_cast<size_t>(intersects[i][tup[i]]) << ((7 - i) * 8));
          }
          HPCombi::BMat8 bm2(o);
          HPCombi::BMat8 mult = left_reflexive_multiplier(bm2, bm, _dim);
          if (mult != bm && bm2 != bm && mult * bm2 == bm) {
            is_maximal = false;
            break;
          }
        } while (increment_tuple(tup, maxes) && !stopped());


        if (stopped()) {
          break_flag = true;
          continue;
        }

        if (is_maximal || bm == one) {
          std::lock_guard<std::mutex> lg(_mtx);
          _filtered.push_back(bm);
          append_bmat_file(_out, bm);
          _remaining.erase(bm);
        } else {
          std::lock_guard<std::mutex> lg(_mtx);
          append_bmat_file(_discarded, bm);
          _remaining.erase(bm);
        }
        count++;
        if (report()) {
          REPORT_DEFAULT("On %d (overall %d out of %d), keeping %d.\n",
                         i + 1,
                         count.load(),
                         bmat_enum.size(),
                         _filtered.size());
        }
      }

      _finished = true;
      report_why_we_stopped();
    }

    const std::vector<HPCombi::BMat8> &reps() {
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

    size_t nr_remaining() {
      std::lock_guard<std::mutex> lg(_mtx);
      return _remaining.size();
    }

    std::vector<HPCombi::BMat8> unprocessed_reps() {
      return std::vector<HPCombi::BMat8>(_remaining.begin(), _remaining.end());
    }

   private:
    std::string                        _in;
    std::string                        _out;
    std::string                        _discarded;
    std::vector<HPCombi::BMat8>        _filtered;
    bool                               _finished;
    std::mutex                         _mtx;
    std::unordered_set<HPCombi::BMat8> _remaining;
  };

  template <size_t _dim>
  class ReflexiveBruteFilterer : public ::libsemigroups::Runner {
   public:
    ReflexiveBruteFilterer(std::vector<HPCombi::BMat8> filtees,
                           std::vector<HPCombi::BMat8> filters)
        : _filtees(filtees),
          _filters(filters),
          _filtered(0),
          _finished(false),
          _mtx() {}

    void run_impl() {
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
        S_bmats.push_back(bmat8_helpers::make<_dim,
                                              typename PermHelper<_dim>::type,
                                              HPCombi::BMat8>(p));
      }

      HPCombi::BMat8 one = bmat8_helpers::one<HPCombi::BMat8>(_dim);

      for (size_t i = 0; i < _filtees.size(); ++i) {
        HPCombi::BMat8              x = _filtees[i];
        std::vector<HPCombi::BMat8> permuted;
        for (HPCombi::BMat8 p : S_bmats) {
          permuted.push_back(p.transpose() * x * p);
        }
        bool found = false;
        for (size_t j = 0; j < _filters.size() && !found; ++j) {
          HPCombi::BMat8 y = _filters[j];
          if (x == y | y == one) {
            continue;
          }
          for (HPCombi::BMat8 z : permuted) {
            if ((z.to_int() | y.to_int() == z.to_int())
                && left_reflexive_multiplier(y, z, _dim) * y == z) {
              found = true;
              break;
            }
          }
        }
        if (!found) {
          _filtered.push_back(x);
        }
        if (report()) {
          REPORT_DEFAULT("On %d out of %d, keeping %d.\n",
                         i + 1,
                         _filtees.size(),
                         _filtered.size());
        }
      }
      _finished = true;
      report_why_we_stopped();
    }

    const std::vector<HPCombi::BMat8> &reps() {
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
    std::vector<HPCombi::BMat8> _filtees;
    std::vector<HPCombi::BMat8> _filters;
    std::vector<HPCombi::BMat8> _filtered;
    bool                        _finished;
    std::mutex                  _mtx;
  };

  template <size_t _dim>
  class ReflexiveOrbiter : public ::libsemigroups::Runner {
   public:
    ReflexiveOrbiter(std::string in, std::string out_pref)
        : _in(in), _out_pref(out_pref), _out(0), _finished() {}
    void run_impl() {
      auto rg = ReportGuard();
      std::vector<HPCombi::BMat8> bmat_enum;
      std::ifstream               f(_in);
      std::string                 line;
      while (std::getline(f, line)) {
        bmat_enum.push_back(HPCombi::BMat8(std::stoul(line)));
      }
      f.close();

      REPORT_DEFAULT("read " + std::to_string(bmat_enum.size())
                     + " matrices to orbit from " + _in + "\n");

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
      std::vector<HPCombi::BMat8> S_bmats_transpose;
      for (Perm p : S) {
        HPCombi::BMat8 x = bmat8_helpers::
            make<_dim, typename PermHelper<_dim>::type, HPCombi::BMat8>(p);
        S_bmats.push_back(x);
        S_bmats_transpose.push_back(x.transpose());
      }

      std::atomic<size_t> count{0};

#pragma omp parallel for schedule(dynamic, 100)
      for (size_t i = 0; i < bmat_enum.size(); ++i) {
        std::unordered_set<HPCombi::BMat8> set(0);
        std::vector<HPCombi::BMat8>        tmp(0);
        HPCombi::BMat8                     bm = bmat_enum[i];
        for (size_t j = 0; j < S_bmats.size(); ++j) {
          HPCombi::BMat8 x = S_bmats[j] * bm * S_bmats_transpose[j];
          if (set.find(x) == set.end()) {
            set.insert(x);
            tmp.push_back(x);
          }
        }
        count++;
        write_bmat_file(
            _out_pref + "_" + std::to_string(i) + ".txt", tmp);

        if (report()) {
          REPORT_DEFAULT("Orbited %d (overall %d out of %d).\n",
                         i + 1,
                         count.load(),
                         bmat_enum.size());
        }
      }

      for (size_t i = 0; i < bmat_enum.size(); ++i) {
        std::vector<HPCombi::BMat8> tmp = read_bmat_file(
            _out_pref + "_" + std::to_string(i) + ".txt");
        _out.insert(_out.end(), tmp.begin(), tmp.end());
        std::remove((_out_pref + "_" + std::to_string(i) + ".txt").c_str());
        if (report()) {
          REPORT_DEFAULT("Gathered %d out of %d; %d generators found.\n",
                         i + 1,
                         bmat_enum.size(),
                         _out.size());
        }
      }
      write_bmat_file(_out_pref + ".txt", _out);
      _finished = true;
      report_why_we_stopped();
    }

    const std::vector<HPCombi::BMat8> &reps() {
      if (!started()) {
        run();
      }
      return _out;
    }

    size_t size() {
      if (!started()) {
        run();
      }
      return _out.size();
    }

    bool finished_impl() const {
      return _finished;
    }

   private:
    std::string                 _in;
    std::string                 _out_pref;
    std::vector<HPCombi::BMat8> _out;
    bool                        _finished;
  };

  ////////////////////////////////////////////////////////////////////////
  // The function actually performing the calculation
  ////////////////////////////////////////////////////////////////////////

  std::vector<HPCombi::BMat8> trim_bmats(size_t n) {
    auto           rg = ReportGuard();
    BMatEnumerator enumerator(n, true);
    std::string    tf
        = "build/output/bmat_enum_trim_" + detail::to_string(n) + ".txt";
    write_bmat_file(tf, enumerator.reps());
    append_bmat_file(tf, bmat8_helpers::elementary<HPCombi::BMat8>(n));
    std::cout << enumerator.reps().size() + 1 << " trim matrices written to:";
    std::cout << "    " << tf << "\n";
    std::vector<HPCombi::BMat8> reps = enumerator.reps();
    reps.push_back(bmat8_helpers::elementary<HPCombi::BMat8>(n));
    return reps;
  }

  void write_row_spaces(std::string const &                out,
                        std::vector<HPCombi::BMat8> const &reps) {
    std::ofstream o(out, std::ios::out | std::ios::trunc);
    for (auto const &x : reps) {
      std::vector<uint8_t> row_space = row_space_vector(x);
      o << "[";
      for (size_t i = 0; i < row_space.size() - 1; ++i) {
        o << size_t(row_space[i]) << ", ";
      }
      o << size_t(row_space.back()) << "]\n";
    }
    o.close();
    std::cout << reps.size() << " row spaces written to:";
    std::cout << "       " << out << "\n";
  }

  void full_candidates_and_row_spaces(size_t n) {
    auto                        rg   = ReportGuard();
    std::vector<HPCombi::BMat8> reps = trim_bmats(n);
    std::string                 rsn
        = "build/output/row_space_numbers_" + detail::to_string(n) + ".txt";
    write_row_spaces(rsn, reps);
  }

  template <size_t n>
  std::vector<HPCombi::BMat8> full_prefilter() {
    auto rg = ReportGuard();
    std::vector<HPCombi::BMat8> reps = read_bmat_file(
        "build/output/bmat_enum_trim_" + detail::to_string(n) + ".txt");

    // the number of largest row/column space sizes to filter on
    size_t NR_SPACE_SIZES = 13;

    std::vector<HPCombi::BMat8> ext_gens;
    std::ifstream f("build/output/bmat-int-gens-" + detail::to_string(n - 1)
                    + ".txt");
    std::string   line;
    if (!f.good()) {
      LIBSEMIGROUPS_EXCEPTION(
          "tried to pre-filter generators for B_" + detail::to_string(n)
          + " without having generators for B_" + detail::to_string(n - 1)
          + "; please compute those generators first \n");
    }
    while (std::getline(f, line)) {
      ext_gens.push_back(HPCombi::BMat8(std::stoul(line) | 1));
    }
    f.close();

    std::set<size_t> space_sizes;
    for (auto it = ext_gens.begin(); it != ext_gens.end(); ++it) {
      space_sizes.insert(row_space_size(*it));
    }

    size_t i = 0;
    std::set<size_t> filtering_space_sizes;
    for (auto it = space_sizes.rbegin();
         it != space_sizes.rend() && i < NR_SPACE_SIZES;
         ++it, ++i) {
      filtering_space_sizes.insert(*it);
    }
   
    std::vector<HPCombi::BMat8> filts;
    for (auto bm : ext_gens) {
      if (filtering_space_sizes.count(row_space_size(bm))) {
        filts.push_back(bm);
      }
    }

    std::string prefilt0
        = "build/output/bmat_enum_trim_" + detail::to_string(n) + ".txt";
    std::string prefilt1
        = "build/output/bmat_full_" + detail::to_string(n) + "_prefilter_1.txt";
    std::string prefilt2
        = "build/output/bmat_full_" + detail::to_string(n) + "_prefilter_2.txt";
    std::string prefilt3
        = "build/output/bmat_full_" + detail::to_string(n) + "_prefilter_3.txt";
    std::string prefilt4
        = "build/output/bmat_full_" + detail::to_string(n) + "_prefiltered.txt";

    {
      Filterer<n> filterer1(prefilt0, prefilt1, filts, false, false);
      filterer1.run();
    }
    {
      Filterer<n> filterer2(prefilt1, prefilt2, filts, false, true);
      filterer2.run();
    }
    {
      Filterer<n> filterer3(prefilt2, prefilt3, filts, true, false);
      filterer3.run();
    }
    Filterer<n> filterer4(prefilt3, prefilt4, filts, true, true);
    filterer4.run();

    return filterer4.reps();
  }
  
  template <size_t n>
  void full_run_no_enum() {
    std::string strn = detail::to_string(n);
    std::string logf
        = "build/output/log-full-" + strn + ".txt";
    std::ofstream o(logf, std::ios::out | std::ios::app);

    auto timer = detail::Timer();
    std::vector<HPCombi::BMat8> reps = full_prefilter<n>();
    std::string timestr = timer.string();
    std::string                 logstr
        = "prefiltered and wrote "
          + detail::to_string(read_bmat_file("build/output/bmat_full_" + strn
                                             + "_prefiltered.txt")
                                  .size())
          + " candidates to build/output/bmat_full_" + strn
          + "_prefiltered.txt in " + timestr + "\n";
    o << logstr << std::endl;
    ;
    std::cout << logstr;

    std::string rsn = "build/output/row_space_numbers_" + detail::to_string(n);
    write_row_spaces(rsn + "_prefiltered.txt", reps);
    write_row_spaces(
        rsn + ".txt",
        read_bmat_file("build/output/bmat_enum_trim_" + strn + ".txt"));
    o << "row spaces written" << std::endl
      << "CPP done" << std::endl
      << std::endl;
  }

  template <size_t n>
  void full_run() {
    auto rg = ReportGuard();
    std::string strn = detail::to_string(n);
    std::string logf
        = "build/output/log-full-" + strn + ".txt";
    std::ofstream o(logf, std::ios::out | std::ios::trunc);
    o << "computing candidates for generators of B_" + strn << std::endl;
    auto timer = detail::Timer();
    trim_bmats(n);
    std::string timestr = timer.string();
    std::string logstr
        = "wrote "
          + detail::to_string(
                read_bmat_file("build/output/bmat_enum_trim_" + strn + ".txt")
                    .size())
          + " trim reps to build/output/bmat_enum_trim_" + strn + ".txt in "
          + timestr + "\n";
    o << logstr;
    std::cout << logstr;
    o.close();
    full_run_no_enum<n>();
  }

  template <size_t n>
  void reflexive_run_no_enum() {
    std::string strn = detail::to_string(n);
    std::string logf
        = "build/output/log-reflexive-" + strn + ".txt";
    std::ofstream o(logf, std::ios::out | std::ios::app);
    std::string   candf = "build/output/bmat_reflexive_candidates_"
                        + detail::to_string(n) + ".txt";
    std::string filtf
        = "build/output/bmat_reflexive_filtered_" + detail::to_string(n) + ".txt";
    std::string discf
        = "build/output/bmat_reflexive_discarded_" + detail::to_string(n) + ".txt";
    std::string gensf
        = "build/output/bmat_reflexive_gens_" + detail::to_string(n) + ".txt";
    std::string asocf = "build/output/bmat_reflexive_associates_" + strn;

    std::string timestr;
    ReflexiveFilterer<n> filterer(candf, filtf, discf);
    auto timer = detail::Timer();
    if (n < 6) {
      filterer.run();
    } else {
      filterer.run_until(
          [&filterer]() { return filterer.nr_remaining() < 12; });
      std::vector<HPCombi::BMat8> remainders = filterer.unprocessed_reps();
      std::cout << "processing " << detail::to_string(remainders.size())
                << " leftovers" << std::endl;
      std::vector<HPCombi::BMat8> filters = read_bmat_file(filtf);
      filters.insert(filters.end(), remainders.begin(), remainders.end());
      ReflexiveBruteFilterer<n> brute_filterer(remainders, filters);
      std::vector<HPCombi::BMat8> extra_reps = brute_filterer.reps();
      for (auto it = extra_reps.begin(); it != extra_reps.end(); ++it) {
        append_bmat_file(filtf, *it);
      }
    }
    timestr = timer.string();
    std::string logstr
        = "filtered " + detail::to_string(read_bmat_file(filtf).size())
          + " reflexive reps in " + timestr + ", written to " + filtf + "\n";
    o << logstr;
    std::cout << logstr;
    timer.reset();
    ReflexiveOrbiter<n> orbiter(filtf, asocf);
    orbiter.run();
    timestr = timer.string();
    logstr  = "orbited "
             + detail::to_string(read_bmat_file(asocf + ".txt").size())
             + " reflexive associates in " + timestr + ", written to " + asocf
             + ".txt\n";
    o << logstr;
    std::cout << logstr;
    write_bmat_file(gensf, orbiter.reps());
    for (size_t i = 0; i < n; ++i) {
      for (size_t j = 0; j < n; ++j) {
        if (i != j) {
          HPCombi::BMat8 bm = bmat8_helpers::one<HPCombi::BMat8>(n);
          bm.set(i, j, true);
          append_bmat_file(gensf, bm);
        }
      }
    }
    logstr = "wrote " + detail::to_string(read_bmat_file(gensf).size())
             + " reflexive generators to " + gensf + "\n";
    o << logstr;
    std::cout << logstr;
  }

  template <size_t n>
  void reflexive_minimal_generating_set() {
    auto        rg    = ReportGuard();
    std::string strn = detail::to_string(n);
    std::string logf
        = "build/output/log-reflexive-" + strn + ".txt";
    std::ofstream o(logf, std::ios::out | std::ios::trunc);
    o << "computing orbit representatives for reflexive maximal elements of B_"
             + strn
      << std::endl;
    std::string candf = "build/output/bmat_reflexive_candidates_"
                        + detail::to_string(n) + ".txt";
    auto                                timer = detail::Timer();
    BMatReflexiveTrimStraightEnumerator enumerator(n);
    enumerator.run();
    std::string                         timestr = timer.string();
    write_bmat_file(candf, enumerator.reps());
    std::string logstr
        = "computed " + detail::to_string(read_bmat_file(candf).size())
          + " reflexive reps in " + timestr + ", written to " + candf + "\n";
    o << logstr;
    std::cout << logstr;
    o.close();
    reflexive_run_no_enum<n>();
  }
  

  /*
  template <size_t n>
  void reflexive_filter_and_orbit() {
    auto        rg    = ReportGuard();
    std::string candf = "build/output/bmat_reflexive_candidates_"
                        + detail::to_string(n) + ".txt";
    std::string filtf = "build/output/bmat_reflexive_filtered_"
                        + detail::to_string(n) + "_filt1.txt";
    std::string discf = "build/output/bmat_reflexive_discarded_"
                        + detail::to_string(n) + "_filt1.txt";
    std::string gensf
        = "build/output/bmat_reflexive_gens_" + detail::to_string(n) + ".txt";
    std::string          asocf = "build/output/bmat_reflexive";
    ReflexiveFilterer<n> filterer(candf, filtf, discf);
    filterer.run();
    ReflexiveOrbiter<n> orbiter(filtf, asocf);
    write_bmat_file(gensf, orbiter.reps());
    for (size_t i = 0; i < n; ++i) {
      for (size_t j = 0; j < n; ++j) {
        if (i != j) {
          HPCombi::BMat8 bm = bmat8_helpers::one<HPCombi::BMat8>(n);
          bm.set(i, j, true);
          append_bmat_file(gensf, bm);
        }
      }
    }
    std::cout << orbiter.size() + (n * n - n) << " generators written to:";
    std::cout << "    " << gensf << "\n";
  }
  */

  ////////////////////////////////////////////////////////////////////////

  LIBSEMIGROUPS_TEST_CASE("BMat8 enum",
                          "3",
                          "enumerate minimal generating set for n = 3",
                          "[standard][enumerate][full-3]") {
    full_run<3>();
  }

  LIBSEMIGROUPS_TEST_CASE("BMat8 enum",
                          "4",
                          "enumerate minimal generating set for n = 4",
                          "[standard][enumerate][full-4]") {
    full_run<4>();
  }

  LIBSEMIGROUPS_TEST_CASE("BMat8 enum",
                          "5",
                          "enumerate minimal generating set for n = 5",
                          "[standard][enumerate][full-5]") {
    full_run<5>();
  }

  LIBSEMIGROUPS_TEST_CASE("BMat8 enum",
                          "6",
                          "enumerate minimal generating set for n = 6",
                          "[standard][enumerate][full-6]") {
    full_run<6>();
  }

  LIBSEMIGROUPS_TEST_CASE("BMat8 enum",
                          "7",
                          "enumerate minimal generating set for n = 7",
                          "[standard][enumerate][full-7]") {
    full_run<7>();
  }

  LIBSEMIGROUPS_TEST_CASE("BMat8 enum",
                          "8",
                          "enumerate minimal generating set for n = 8",
                          "[standard][enumerate][full-8]") {
    full_run<8>();
  }

  /////////////////////////////////////////////////////////////////////////////
  //
  // Reflexive BMat8 Generators
  //
  /////////////////////////////////////////////////////////////////////////////

  LIBSEMIGROUPS_TEST_CASE("BMat8 enum",
                          "13",
                          "enumerate minimal generating set for reflexive "
                          "boolean mat monoid with n = 3",
                          "[standard][enumerate][ref-3]") {
    reflexive_minimal_generating_set<3>();
  }

  LIBSEMIGROUPS_TEST_CASE("BMat8 enum",
                          "14",
                          "enumerate minimal generating set for reflexive "
                          "boolean mat monoid with n = 4",
                          "[standard][enumerate][ref-4]") {
    reflexive_minimal_generating_set<4>();
  }

  LIBSEMIGROUPS_TEST_CASE("BMat8 enum",
                          "15",
                          "enumerate minimal generating set for reflexive "
                          "boolean mat monoid with n = 5",
                          "[standard][enumerate][ref-5]") {
    reflexive_minimal_generating_set<5>();
  }

  LIBSEMIGROUPS_TEST_CASE("BMat8 enum",
                          "16",
                          "enumerate minimal generating set for reflexive "
                          "boolean mat monoid with n = 6",
                          "[standard][enumerate][ref-6]") {
    reflexive_minimal_generating_set<6>();
  }

  LIBSEMIGROUPS_TEST_CASE("BMat8 enum",
                          "16",
                          "enumerate minimal generating set for reflexive "
                          "boolean mat monoid with n = 7",
                          "[standard][enumerate][ref-7]") {
    reflexive_minimal_generating_set<7>();
  }

  /*
  LIBSEMIGROUPS_TEST_CASE("BMat8 enum",
                          "17",
                          "enumerate minimal generating set for reflexive "
                          "boolean mat monoid with n = 7",
                          "[standard][enumerate]") {
    reflexive_minimal_generating_set<7>();
  }

  LIBSEMIGROUPS_TEST_CASE("BMat8 enum",
                          "18",
                          "filter reflexive 7 reps - delete this when done"
                          "boolean mat monoid with n = 7",
                          "[standard][enumerate]") {
    reflexive_filter_and_orbit<7>();
  }

  LIBSEMIGROUPS_TEST_CASE(
      "BMat8 enum",
      "20",
      "filter the leftover matrices from the reflexive monoid"
      "boolean mat monoid with n = 7",
      "[standard][enumerate]") {
    std::vector<HPCombi::BMat8> filters
        = read_bmat_file("../output/saved/bmat_reflexive_part_filtered_7.txt");
    write_bmat_file("../output/saved/bmat_reflexive_fully_filtered_7.txt",
                    filters); 
    std::vector<HPCombi::BMat8> filtees =
        read_bmat_file("../output/saved/bmat_reflexive_leftovers_7.txt");
    filters.insert(filters.end(), filtees.begin(), filtees.end());
    ReflexiveBruteFilterer<7> filt(filtees, filters);
    filt.run();
    write_bmat_file("../output/bmat_reflexive_leftovers_filtered_7.txt",
                    filt.reps());
    for (HPCombi::BMat8 bm : filt.reps()) {
      append_bmat_file("../output/saved/bmat_reflexive_fully_filtered_7.txt",
                       bm);
    }
  }

  LIBSEMIGROUPS_TEST_CASE("BMat8 enum",
                          "21",
                          "orbit the reflexive 7 reps"
                          "boolean mat monoid with n = 7",
                          "[standard][enumerate]") {
    std::string gensf = "../output/saved/bmat_reflexive_gens_7.txt";
    ReflexiveOrbiter<7> orbiter(
        "../output/saved/bmat_reflexive_fully_filtered_7.txt",
        "../output/saved/bmat_reflexive");
    write_bmat_file(gensf, orbiter.reps());
    for (size_t i = 0; i < 7; ++i) {
      for (size_t j = 0; j < 7; ++j) {
        if (i != j) {
          HPCombi::BMat8 bm = bmat8_helpers::one<HPCombi::BMat8>(7);
          bm.set(i, j, true);
          append_bmat_file(gensf, bm);
        }
      }
    }
    std::cout << orbiter.size() << std::endl;
    std::cout << orbiter.size() + (7 * 6)
              << " generators written to:";
    std::cout << "    " << gensf << "\n";
  }
*/

  /////////////////////////////////////////////////////////////////////////////
  //
  // Filter the trim matrices in Bn by brute force
  //
  /////////////////////////////////////////////////////////////////////////////

  LIBSEMIGROUPS_TEST_CASE("BMat8 enum",
                          "33",
                          "brute filter B3 - with prefiltering",
                          "[extreme][brute-filter-3]") {
    std::string bmatf = "build/output/bmat_full_3_prefiltered.txt";
    std::string outf = "build/output/bmat_full_3_brute_prefiltered.txt";
    std::vector<HPCombi::BMat8> filts = read_bmat_file(bmatf);

    Filterer<3> filterer(bmatf,
                         outf,
                         filts,
                         true);
    filterer.run();
    std::cout << filterer.reps().size() << " trim matrices written to:";
    std::cout << "    " << outf << "\n";
  }
  
  LIBSEMIGROUPS_TEST_CASE("BMat8 enum",
                          "34",
                          "brute filter B4 - with prefiltering",
                          "[extreme][brute-filter-4]") {
    std::string bmatf = "build/output/bmat_full_4_prefiltered.txt";
    std::string outf = "build/output/bmat_full_4_brute_prefiltered.txt";
    std::vector<HPCombi::BMat8> filts = read_bmat_file(bmatf);

    Filterer<4> filterer(bmatf,
                         outf,
                         filts,
                         true);
    filterer.run();
    std::cout << filterer.reps().size() << " trim matrices written to:";
    std::cout << "    " << outf << "\n";
  }
  
  LIBSEMIGROUPS_TEST_CASE("BMat8 enum",
                          "35",
                          "brute filter B5 - with prefiltering",
                          "[extreme][brute-filter-5]") {
    std::string bmatf = "build/output/bmat_full_5_prefiltered.txt";
    std::string outf = "build/output/bmat_full_5_brute_prefiltered.txt";
    std::vector<HPCombi::BMat8> filts = read_bmat_file(bmatf);

    Filterer<5> filterer(bmatf,
                         outf,
                         filts,
                         true);
    filterer.run();
    std::cout << filterer.reps().size() << " trim matrices written to:";
    std::cout << "    " << outf << "\n";
  }
  
  LIBSEMIGROUPS_TEST_CASE("BMat8 enum",
                          "36",
                          "brute filter B6 - with prefiltering",
                          "[extreme][brute-filter-6]") {
    std::string bmatf = "build/output/bmat_full_6_prefiltered.txt";
    std::string outf = "build/output/bmat_full_6_brute_prefiltered.txt";
    std::vector<HPCombi::BMat8> filts = read_bmat_file(bmatf);

    Filterer<6> filterer(bmatf,
                         outf,
                         filts,
                         true);
    filterer.run();
    std::cout << filterer.reps().size() << " trim matrices written to:";
    std::cout << "    " << outf << "\n";
  }

  LIBSEMIGROUPS_TEST_CASE("BMat8 enum",
                          "37",
                          "brute filter B7 - with prefiltering",
                          "[extreme][brute-filter-7]") {
    auto                        rg   = ReportGuard();
    std::string bmatf = "build/output/bmat_full_7_prefiltered.txt";
    std::string outf = "build/output/bmat_full_7_brute_prefiltered.txt";
    std::vector<HPCombi::BMat8> filts = read_bmat_file(bmatf);

    Filterer<7> filterer(bmatf,
                         outf,
                         filts,
                         true);
    filterer.run();
    std::cout << filterer.reps().size() << " trim matrices written to:";
    std::cout << "    " << outf << "\n";
  }

  LIBSEMIGROUPS_TEST_CASE("BMat8 enum",
                          "43",
                          "brute filter B3 - with prefiltering",
                          "[extreme][brute-filter-3]") {
    std::string bmatf = "build/output/bmat_enum_trim_3.txt";
    std::string outf = "build/output/bmat_full_3_brute.txt";
    std::vector<HPCombi::BMat8> filts = read_bmat_file(bmatf);

    Filterer<3> filterer(bmatf,
                         outf,
                         filts,
                         true);
    filterer.run();
    std::cout << filterer.reps().size() << " trim matrices written to:";
    std::cout << "    " << outf << "\n";
  }
  
  LIBSEMIGROUPS_TEST_CASE("BMat8 enum",
                          "44",
                          "brute filter B4 - with prefiltering",
                          "[extreme][brute-filter-4]") {
    std::string bmatf = "build/output/bmat_enum_trim_4.txt";
    std::string outf = "build/output/bmat_full_4_brute.txt";
    std::vector<HPCombi::BMat8> filts = read_bmat_file(bmatf);

    Filterer<4> filterer(bmatf,
                         outf,
                         filts,
                         true);
    filterer.run();
    std::cout << filterer.reps().size() << " trim matrices written to:";
    std::cout << "    " << outf << "\n";
  }
  
  LIBSEMIGROUPS_TEST_CASE("BMat8 enum",
                          "45",
                          "brute filter B5 - with prefiltering",
                          "[extreme][brute-filter-5]") {
    std::string bmatf = "build/output/bmat_enum_trim_5.txt";
    std::string outf = "build/output/bmat_full_5_brute.txt";
    std::vector<HPCombi::BMat8> filts = read_bmat_file(bmatf);

    Filterer<5> filterer(bmatf,
                         outf,
                         filts,
                         true);
    filterer.run();
    std::cout << filterer.reps().size() << " trim matrices written to:";
    std::cout << "    " << outf << "\n";
  }
  
  LIBSEMIGROUPS_TEST_CASE("BMat8 enum",
                          "46",
                          "brute filter B6 - with prefiltering",
                          "[extreme][brute-filter-6]") {
    std::string bmatf = "build/output/bmat_enum_trim_6.txt";
    std::string outf = "build/output/bmat_full_6_brute.txt";
    std::vector<HPCombi::BMat8> filts = read_bmat_file(bmatf);

    Filterer<6> filterer(bmatf,
                         outf,
                         filts,
                         true);
    filterer.run();
    std::cout << filterer.reps().size() << " trim matrices written to:";
    std::cout << "    " << outf << "\n";
  }

  LIBSEMIGROUPS_TEST_CASE("BMat8 enum",
                          "47",
                          "brute filter B7 - with prefiltering",
                          "[extreme][brute-filter-7]") {
    auto                        rg   = ReportGuard();
    std::string bmatf = "build/output/bmat_enum_trim_7.txt";
    std::string outf = "build/output/bmat_full_7_brute.txt";
    std::vector<HPCombi::BMat8> filts = read_bmat_file(bmatf);

    Filterer<7> filterer(bmatf,
                         outf,
                         filts,
                         true);
    filterer.run();
    std::cout << filterer.reps().size() << " trim matrices written to:";
    std::cout << "    " << outf << "\n";
  }
}  // namespace libsemigroups
