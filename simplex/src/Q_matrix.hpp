#pragma once

#include <iostream>
#include <memory>
#include <set>

// A symmetric 0,1 off-diagonal matrix
class Q_matrix {
public:
  Q_matrix(unsigned n);

  ~Q_matrix();
  Q_matrix(const Q_matrix& other);
  Q_matrix(Q_matrix&& other);
  Q_matrix& operator=(const Q_matrix& other);
  Q_matrix& operator=(Q_matrix&& other);

  int entry(unsigned h1, unsigned h2) const;

  // XOR row/column k into row/column h
  void add_rowcol(unsigned h, unsigned k);

  // Swap row/column with row/column r-1
  void swap_rowcol(unsigned h);

  // Set of h s.t. Q[h][r-1] = 1
  const std::set<unsigned> rows_with_terminal_1() const;

  // Flip the (h1, h2) entry for all h1, h2 in H
  void flip_submatrix(const std::set<unsigned>& H);

  // Flip the (h1, h2) entry for all h1 in H1, h2 in H2 (preserving symmetry)
  void flip_submatrix(
    const std::set<unsigned>& H1, const std::set<unsigned>& H2);

  // Append new row/column (given as set of indices where 1)
  void append_rowcol(const std::set<unsigned>& H);

  // Whether a given row/column is all-zero
  bool rowcol_is_zero(unsigned h) const;

  void drop_final_rowcol();

  unsigned r() const;

  friend std::ostream& operator<<(std::ostream& os, const Q_matrix& Q);

private:
  struct impl;
  std::unique_ptr<impl> pImpl;
};
