#include "A_matrix.hpp"
#include <algorithm>
#include <iostream>
#include <memory>
#include <set>
#include <vector>

#if defined (SIMPLEX_DENSE)

struct A_matrix::impl {
  impl(unsigned n) : n(n), r(0), data(n, std::vector<int>(n+1, 0)) {}

  /* Data */

  unsigned n;
  unsigned r;
  std::vector<std::vector<int>> data;

  /* Methods */

  int entry(unsigned j, unsigned h) const { return data[j][h]; }

  void add_col(unsigned h, unsigned k) {
    for (unsigned j = 0; j < n; j++) {
      data[j][h] ^= data[j][k];
    }
  }

  void add_row(unsigned j, unsigned k) {
    for (unsigned h = 0; h < r; h++) {
      data[j][h] ^= data[k][h];
    }
  }

  unsigned row_weight(unsigned j) const {
    unsigned c = 0;
    for (unsigned h = 0; h < r; h++) {
      if (data[j][h]) c++;
    }
    return c;
  }

  unsigned col_weight(unsigned h) const {
    unsigned c = 0;
    for (unsigned j = 0; j < n; j++) {
      if (data[j][h]) c++;
    }
    return c;
  }

  void swap_col(unsigned h) {
    const unsigned r1 = r - 1;
    for (unsigned j = 0; j < n; j++) {
      std::vector<int>& A_j = data[j];
      std::iter_swap(A_j.begin() + h, A_j.begin() + r1);
    }
  }

  void zero_append_basis_col(unsigned j) {
    std::vector<int>& A_j = data[j];
    for (unsigned h = 0; h < r; h++) {
      A_j[h] = 0;
    }
    for (unsigned k = 0; k < n; k++) {
      data[k][r] = 0;
    }
    A_j[r] = 1;
    r++;
  }

  const std::set<unsigned> cols_where_one(unsigned j) const {
    const std::vector<int>& A_j = data[j];
    std::set<unsigned> H;
    for (unsigned h = 0; h < r; h++) {
      if (A_j[h]) {
        H.insert(h);
      }
    }
    return H;
  }

  const std::set<unsigned> cols_where_one(unsigned j, unsigned k) const {
    const std::vector<int>& A_j = data[j];
    const std::vector<int>& A_k = data[k];
    std::set<unsigned> H;
    for (unsigned h = 0; h < r; h++) {
      if (A_j[h] & A_k[h]) {
        H.insert(h);
      }
    }
    return H;
  }

  void drop_final_col() { r--; }
};

#else

struct A_matrix::impl {
  impl(unsigned n) : n(n), r(0), rows(n), cols(n+1) {}

  /* Data */

  unsigned n;
  unsigned r;
  std::vector<std::set<unsigned>> rows;
  std::vector<std::set<unsigned>> cols;

  /* Methods */

  int entry(unsigned j, unsigned h) const { return rows[j].contains(h); }

  void add_col(unsigned h, unsigned k) {
    std::set<unsigned> newcol;
    for (unsigned j : cols[h]) {
      if (!cols[k].contains(j)) {
        newcol.insert(j);
      }
    }
    for (unsigned j : cols[k]) {
      if (!cols[h].contains(j)) {
        newcol.insert(j);
      }
    }
    for (unsigned j : cols[h]) {
      if (!newcol.contains(j)) {
        rows[j].erase(h);
      }
    }
    for (unsigned j : newcol) {
      rows[j].insert(h);
    }
    cols[h] = newcol;
  }

  void add_row(unsigned j, unsigned k) {
    std::set<unsigned> newrow;
    for (unsigned h : rows[j]) {
      if (!rows[k].contains(h)) {
        newrow.insert(h);
      }
    }
    for (unsigned h : rows[k]) {
      if (!rows[j].contains(h)) {
        newrow.insert(h);
      }
    }
    for (unsigned h : rows[j]) {
      if (!newrow.contains(h)) {
        cols[h].erase(j);
      }
    }
    for (unsigned h : newrow) {
      cols[h].insert(j);
    }
    rows[j] = newrow;
  }

  unsigned row_weight(unsigned j) const {
    return rows[j].size();
  }

  unsigned col_weight(unsigned h) const {
    return cols[h].size();
  }

  void swap_col(unsigned h) {
    for (unsigned j : cols[h]) {
      if (!cols[r - 1].contains(j)) {
        rows[j].erase(h);
        rows[j].insert(r - 1);
      }
    }
    for (unsigned j : cols[r - 1]) {
      if (!cols[h].contains(j)) {
        rows[j].erase(r - 1);
        rows[j].insert(h);
      }
    }
    std::iter_swap(cols.begin() + h, cols.begin() + (r - 1));
  }

  void zero_append_basis_col(unsigned j) {
    for (unsigned h : rows[j]) {
      cols[h].erase(j);
    }
    rows[j].clear();
    rows[j].insert(r);
    cols[r].insert(j);
    r++;
  }

  const std::set<unsigned> cols_where_one(unsigned j) const {
    return rows[j];
  }

  const std::set<unsigned> cols_where_one(unsigned j, unsigned k) const {
    std::set<unsigned> l;
    for (unsigned h : rows[j]) {
      if (rows[k].contains(h)) {
        l.insert(h);
      }
    }
    return l;
  }

  void drop_final_col() {
    for (unsigned j : cols[r - 1]) {
      rows[j].erase(r - 1);
    }
    cols[r - 1].clear();
    r--;
  }
};

#endif

A_matrix::A_matrix(unsigned n) : pImpl(std::make_unique<impl>(n)) {}
A_matrix::~A_matrix() = default;
A_matrix::A_matrix(const A_matrix& other)
  : pImpl(std::make_unique<impl>(*other.pImpl)) {}
A_matrix::A_matrix(A_matrix&& other) = default;
A_matrix& A_matrix::operator=(const A_matrix& other) {
  return *this = A_matrix(other);
}
A_matrix& A_matrix::operator=(A_matrix&& other) = default;

int A_matrix::entry(unsigned j, unsigned h) const { return pImpl->entry(j, h); }
void A_matrix::add_col(unsigned h, unsigned k) { pImpl->add_col(h, k); }
void A_matrix::add_row(unsigned j, unsigned k) { pImpl->add_row(j, k); }
unsigned A_matrix::row_weight(unsigned j) const { return pImpl->row_weight(j); }
unsigned A_matrix::col_weight(unsigned j) const { return pImpl->col_weight(j); }
void A_matrix::swap_col(unsigned h) { pImpl->swap_col(h); }
void A_matrix::zero_append_basis_col(unsigned j) {
  pImpl->zero_append_basis_col(j);
}
const std::set<unsigned> A_matrix::cols_where_one(unsigned j) const {
  return pImpl->cols_where_one(j);
}
const std::set<unsigned> A_matrix::cols_where_one(unsigned j, unsigned k) const {
  return pImpl->cols_where_one(j, k);
}
void A_matrix::drop_final_col() { pImpl->drop_final_col(); }

std::ostream& operator<<(std::ostream& os, const A_matrix& A) {
  for (unsigned j = 0; j < A.pImpl->n; j++) {
    os << "[ ";
    for (unsigned h = 0; h < A.pImpl->r; h++) {
      os << A.pImpl->entry(j, h) << " ";
    }
    os << "]" << std::endl;
  }
  return os;
}
