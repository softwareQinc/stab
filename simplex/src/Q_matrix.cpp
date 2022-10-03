#include "Q_matrix.hpp"
#include <algorithm>
#include <iostream>
#include <memory>
#include <set>
#include <vector>

#if defined (SIMPLEX_DENSE)

struct Q_matrix::impl {
  impl(unsigned n) : n(n), r(0), data(n+1, std::vector<int>(n+1, 0)) {}

  /* Data */

  unsigned n;
  unsigned r;
  std::vector<std::vector<int>> data;

  /* Methods */

  int entry(unsigned h1, unsigned h2) const { return data[h1][h2]; }

  void add_rowcol(unsigned h, unsigned k) {
    for (unsigned j = 0; j < r; j++) {
      data[j][h] ^= data[j][k];
    }
    for (unsigned j = 0; j < r; j++) {
      data[h][j] ^= data[k][j];
    }
  }

  void swap_rowcol(unsigned h) {
    const unsigned r1 = r - 1;
    for (unsigned j = 0; j < r; j++) {
      int t = data[j][h];
      data[j][h] = data[j][r1];
      data[j][r1] = t;
    }
    for (unsigned j = 0; j < r; j++) {
      int t = data[h][j];
      data[h][j] = data[r1][j];
      data[r1][j] = t;
    }
  }

  const std::set<unsigned> rows_with_terminal_1() const {
    std::set<unsigned> H;
    unsigned r1 = r - 1;
    for (unsigned h = 0; h < r1; h++) {
      if (data[h][r1]) H.insert(h);
    }
    return H;
  }

  void flip_submatrix(const std::set<unsigned>& H) {
    for (unsigned h1 : H) {
      for (unsigned h2 : H) {
        if (h1 != h2) {
          data[h1][h2] ^= 1;
        }
      }
    }
  }

  void flip_submatrix(
    const std::set<unsigned>& H1, const std::set<unsigned>& H2)
  {
    for (unsigned h1 : H1) {
      for (unsigned h2 : H2) {
        data[h1][h2] ^= 1;
        data[h2][h1] ^= 1;
      }
    }
  }

  void append_rowcol(const std::set<unsigned>& H) {
    for (unsigned h = 0; h < r; h++) {
      data[h][r] = data[r][h] = 0;
    }
    for (unsigned h : H) {
      data[h][r] = data[r][h] = 1;
    }
    r++;
  }

  bool rowcol_is_zero(unsigned h) const {
    for (unsigned j = 0; j < r; j++) {
      if (j != h && data[h][j]) {
        return false;
      }
    }
    return true;
  }

  void drop_final_rowcol() { r--; }
};

#else

struct Q_matrix::impl {
  impl(unsigned n) : n(n), r(0), rows(n+1) {}

  /* Data */

  unsigned n;
  unsigned r;
  std::vector<std::set<unsigned>> rows;
  std::vector<std::set<unsigned>> cols;

  /* Methods */

  int entry(unsigned h1, unsigned h2) const { return rows[h1].contains(h2); }

  void add_rowcol(unsigned h, unsigned k) {
    // 1. Compute what the new row will be:
    std::set<unsigned> newrow;
    for (unsigned j : rows[h]) {
      if (!rows[k].contains(j)) {
        newrow.insert(j);
      }
    }
    for (unsigned j : rows[k]) {
      if (j != h && !rows[h].contains(j)) {
        newrow.insert(j);
      }
    }
    // 2. Fix up the other rows for symmetry:
    for (unsigned j : rows[h]) {
      if (!newrow.contains(j)) {
        rows[j].erase(h);
      }
    }
    for (unsigned j : newrow) {
      rows[j].insert(h);
    }
    // 3. Replace the row:
    rows[h] = newrow;
  }

  void swap_rowcol(unsigned h) {
    std::set<unsigned> D1(rows[h]);
    for (unsigned k : rows[r - 1]) {
      D1.erase(k);
    }
    D1.erase(h); D1.erase(r - 1);
    std::set<unsigned> D2(rows[r - 1]);
    for (unsigned k : rows[h]) {
      D2.erase(k);
    }
    D2.erase(h); D2.erase(r - 1);
    for (unsigned k : D1) {
      rows[h].erase(k); rows[k].erase(h);
      rows[r - 1].insert(k); rows[k].insert(r - 1);
    }
    for (unsigned k : D2) {
      rows[r - 1].erase(k); rows[k].erase(r - 1);
      rows[h].insert(k); rows[k].insert(h);
    }
  }

  const std::set<unsigned> rows_with_terminal_1() const {
    return rows[r - 1];
  }

  void flip_submatrix(const std::set<unsigned>& H) {
    for (unsigned h1 : H) {
      for (unsigned h2 : H) {
        if (h1 < h2) {
          if (rows[h1].contains(h2)) {
            rows[h1].erase(h2);
            rows[h2].erase(h1);
          } else {
            rows[h1].insert(h2);
            rows[h2].insert(h1);
          }
        }
      }
    }
  }

  void flip_submatrix(
    const std::set<unsigned>& H1, const std::set<unsigned>& H2)
  {
    for (unsigned h1 : H1) {
      for (unsigned h2 : H2) {
        if (h1 != h2) {
          if (rows[h1].contains(h2)) {
            rows[h1].erase(h2);
            rows[h2].erase(h1);
          } else {
            rows[h1].insert(h2);
            rows[h2].insert(h1);
          }
        }
      }
    }
  }

  void append_rowcol(const std::set<unsigned>& H) {
    for (unsigned h : H) {
      rows[r].insert(h);
      rows[h].insert(r);
    }
    r++;
  }

  bool rowcol_is_zero(unsigned h) const {
    unsigned sz = rows[h].size();
    if (sz == 0) {
      return true;
    } else if (sz == 1) {
      return rows[h].contains(h);
    } else {
      return false;
    }
  }

  void drop_final_rowcol() {
    for (unsigned h : rows[r - 1]) {
      if (h != r - 1) {
        rows[h].erase(r - 1);
      }
    }
    rows[r - 1].clear();
    r--;
  }
};

#endif

Q_matrix::Q_matrix(unsigned n) : pImpl(std::make_unique<impl>(n)) {}
Q_matrix::~Q_matrix() = default;
Q_matrix::Q_matrix(const Q_matrix& other)
  : pImpl(std::make_unique<impl>(*other.pImpl)) {}
Q_matrix::Q_matrix(Q_matrix&& other) = default;
Q_matrix& Q_matrix::operator=(const Q_matrix& other) {
  return *this = Q_matrix(other);
}
Q_matrix& Q_matrix::operator=(Q_matrix&& other) = default;

int Q_matrix::entry(unsigned h1, unsigned h2) const {
  return pImpl->entry(h1, h2);
}
void Q_matrix::add_rowcol(unsigned h, unsigned k) { pImpl->add_rowcol(h, k); }
void Q_matrix::swap_rowcol(unsigned h) { pImpl->swap_rowcol(h); }
const std::set<unsigned> Q_matrix::rows_with_terminal_1() const {
  return pImpl->rows_with_terminal_1();
}
void Q_matrix::flip_submatrix(const std::set<unsigned>& H) {
  pImpl->flip_submatrix(H);
}
void Q_matrix::flip_submatrix(
  const std::set<unsigned>& H1, const std::set<unsigned>& H2) {
  pImpl->flip_submatrix(H1, H2);
}
void Q_matrix::append_rowcol(const std::set<unsigned>& H) {
  pImpl->append_rowcol(H);
}
bool Q_matrix::rowcol_is_zero(unsigned h) const {
  return pImpl->rowcol_is_zero(h);
}
void Q_matrix::drop_final_rowcol() {pImpl->drop_final_rowcol(); }
unsigned Q_matrix::r() const { return pImpl->r; }

std::ostream& operator<<(std::ostream& os, const Q_matrix& Q) {
  for (unsigned j = 0; j < Q.pImpl->r; j++) {
    os << "[ ";
    for (unsigned h = 0; h < Q.pImpl->r; h++) {
      os << Q.pImpl->entry(j, h) << " ";
    }
    os << "]" << std::endl;
  }
  return os;
}
