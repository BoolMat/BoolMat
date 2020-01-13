

  bool is_trim(HPCombi::BMat8 bm, size_t dim = 8) {
    return is_row_trim(bm, dim) && is_col_trim(bm, dim);
  }
