

  bool is_trim(HPCombi::BMat8 bm, size_t dim = 8) {
    return is_row_trim(bm, dim) && is_col_trim(bm, dim);
  }

  ////////////////////////////////////////////////////////////////////////
  // The next test cases are just tests that things work properly
  ////////////////////////////////////////////////////////////////////////

  LIBSEMIGROUPS_TEST_CASE("BMat8 enum",
                          "001",
                          "test bliss canonicalisation",
                          "[quick][bliss]") {
    auto           rg  = ReportGuard(false);
    HPCombi::BMat8 x   = HPCombi::BMat8({{0, 0, 0}, {0, 0, 1}, {0, 1, 1}});
    HPCombi::BMat8 y   = HPCombi::BMat8({{0, 0, 0}, {0, 0, 1}, {1, 0, 1}});
    bliss_digraph  dgx = bliss_digraph_from_BMat8(x, 3);
    bliss_digraph  dgy = bliss_digraph_from_BMat8(y, 3);

    // dgx.write_dot("dgx.dot");

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
    auto           rg = ReportGuard(false);
    BMatEnumerator enumerator(4, false);
    enumerator.report_every(std::chrono::nanoseconds(1));
    REQUIRE(enumerator.reps().size() == 60);

    std::ofstream o;
    o.open("build/output/bmat_enum_4.txt", std::ios::out | std::ios::trunc);
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
    o.open("build/output/bmat_enum_trim_5.txt",
           std::ios::out | std::ios::trunc);
    for (HPCombi::BMat8 i : enumerator.reps()) {
      o << i.to_int() << "\n";
    }
    o.close();
  }

  LIBSEMIGROUPS_TEST_CASE("BMat8 enum",
                          "800",
                          "enumerate trim B8",
                          "[extreme][enumerate]") {
    std::ofstream  o;
    auto           rg = ReportGuard();
    BMatEnumerator enumerator_8_trim(8, true);
    o.open("build/output/bmat_trim_enum_8.txt",
           std::ios::out | std::ios::trunc);
    for (HPCombi::BMat8 i : enumerator_8_trim.reps()) {
      o << i.to_int() << "\n";
    }
    o.close();
    REQUIRE(enumerator_8_trim.reps().size() == 17120845);
  }

  LIBSEMIGROUPS_TEST_CASE("BMat8 enum", "006", "filter 5", "[quick]") {
    Filterer<5> f("build/output/bmat_enum_trim_5.txt",
                  "build/output/bmat_filtered_5.txt",
                  {},
                  true);
    f.run();
    std::vector<HPCombi::BMat8> filtered = f.reps();
    REQUIRE(filtered.size() == 9);
  }

  LIBSEMIGROUPS_TEST_CASE("BMat8 enum", "007", "filter 6", "[quick]") {}

  LIBSEMIGROUPS_TEST_CASE("BMat8 enum", "008", "filter 7", "[extreme]") {
    auto        rg = ReportGuard();
    Filterer<7> f("build/output/bmat_trim_enum_7.txt",
                  "build/output/bmat_filtered_7.txt",
                  {},
                  true);
    f.run();
    std::vector<HPCombi::BMat8> filtered = f.reps();
    REQUIRE(filtered.size() == 2139);
  }

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
    o.open("build/output/row_spaces_digraphs_5.txt",
           std::ios::out | std::ios::trunc);
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
    o.open("build/output/row_spaces_digraphs_6.txt",
           std::ios::out | std::ios::trunc);
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
    o.open("build/output/row_spaces_digraphs_7.txt",
           std::ios::out | std::ios::trunc);
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
    o.open("build/output/row_spaces_digraphs_primes_5.txt",
           std::ios::out | std::ios::trunc);
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
    o.open("build/output/row_spaces_digraphs_primes_6.txt",
           std::ios::out | std::ios::trunc);
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
    o.open("build/output/row_spaces_digraphs_primes_7.txt",
           std::ios::out | std::ios::trunc);
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

  // LIBSEMIGROUPS_TEST_CASE("BMat8 enum", "160", "print matrices", "[quick]") {
  //   std::cout << HPCombi::BMat8(13853107707017691136ull) << std::endl;
  //   std::cout << HPCombi::BMat8(4620710844303409152) << std::endl;
  //   std::cout << HPCombi::BMat8(4647750068672397312) << std::endl;
  //   std::cout << HPCombi::BMat8(9241421688590303232ull) << std::endl;

  //   std::cout << bmat8_helpers::elementary(6).to_int() << std::endl;
  //   std::cout << bmat8_helpers::one(5).to_int() << std::endl;
  // }

  LIBSEMIGROUPS_TEST_CASE("BMat8 enum", "170", "prime extensions", "[fails]") {
    std::vector<HPCombi::BMat8> bmat_enum;
    std::ifstream               f("bmat_gens_7.txt");
    std::string                 line;
    while (std::getline(f, line)) {
      bmat_enum.push_back(HPCombi::BMat8(std::stoul(line)));
    }
    f.close();
    HPCombi::BMat8 x = bmat_enum[1021];
    for (HPCombi::BMat8 &bm : simple_prime_extensions(x, 7)) {
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
    o.open("build/output/bmat8_filterers.txt", std::ios::out | std::ios::trunc);
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
                          "[fails]") {
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
                          "[fails]") {
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
