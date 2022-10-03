#pragma once

#include <iostream>
#include <memory>
#include <set>

class A_matrix {
public:
  A_matrix(unsigned n);

  ~A_matrix();
  A_matrix(const A_matrix& other);
  A_matrix(A_matrix&& other);
  A_matrix& operator=(const A_matrix& other);
  A_matrix& operator=(A_matrix&& other);

  int entry(unsigned j, unsigned h) const;

  // XOR column k into column h
  void add_col(unsigned h, unsigned k);

  // XOR row k into row j
  void add_row(unsigned j, unsigned k);

  // Number of elements in row j containing 1
  unsigned row_weight(unsigned j) const;

  // Number of elements in column h containing 1
  unsigned col_weight(unsigned h) const;

  // Swap column h with column r-1
  void swap_col(unsigned h);

  // Zero row j and append new column e_j
  void zero_append_basis_col(unsigned j);

  // Set of column indices h s.t. A[j,h] = 1
  const std::set<unsigned> cols_where_one(unsigned j) const;

  // Set of column indices h s.t. A[j,h] = A[j,k] = 1
  const std::set<unsigned> cols_where_one(unsigned j, unsigned k) const;

  void drop_final_col();

  friend std::ostream& operator<<(std::ostream& os, const A_matrix& A);

private:
  struct impl;
  std::unique_ptr<impl> pImpl;
};
