#include "simplex.hpp"
#include "A_matrix.hpp"
#include "Q_matrix.hpp"
#include "bimap.hpp"
#include "parse-stim.hpp"

#include <algorithm>
#include <iostream>
#include <memory>
#include <optional>
#include <random>
#include <set>
#include <vector>

/* Implementation */

class RBG {
public:
  RBG(int seed = 0) : gen(seed), distrib(0, 1) {}
  int get() { return distrib(gen); }
private:
  std::mt19937 gen;
  std::uniform_int_distribution<> distrib;
};

struct Simplex::impl {
  impl(unsigned n, int seed = 0)
    : n(n), r(0), A(n), b(n, 0), Q(n), R0(n+1, 0), R1(n+1, 0), p(), g(0),
    deterministic(true), rbg(seed) {}

  impl(struct instrs is, int seed = 0) : impl(is.n, seed) {
    for (const auto &op : is.ops) {
      switch (op.type) {
        case optype::X: SimulateX(op.qubits[0]); break;
        case optype::Y: SimulateY(op.qubits[0]); break;
        case optype::Z: SimulateZ(op.qubits[0]); break;
        case optype::H: SimulateH(op.qubits[0]); break;
        case optype::S: SimulateS(op.qubits[0]); break;
        case optype::Sdg: SimulateSdg(op.qubits[0]); break;
        case optype::CX: SimulateCX(op.qubits[0], op.qubits[1]); break;
        case optype::CZ: SimulateCZ(op.qubits[0], op.qubits[1]); break;
        case optype::MeasX: SimulateMeasX(op.qubits[0]); break;
        case optype::MeasY: SimulateMeasY(op.qubits[0]); break;
        case optype::MeasZ: SimulateMeasZ(op.qubits[0]); break;
        case optype::ResetX: SimulateResetX(op.qubits[0]); break;
        case optype::ResetY: SimulateResetY(op.qubits[0]); break;
        case optype::ResetZ: SimulateResetZ(op.qubits[0]); break;
        default:
          std::cerr << "Unrecognized operation" << std::endl;
          throw;
      }
    }
  }

  impl(const char *p, int seed = 0)
    : impl(parse_file(p), seed) {}

  /* Data */

  unsigned n;
  unsigned r;
  A_matrix A;
  std::vector<int> b;
  Q_matrix Q;
  std::vector<int> R0;
  std::vector<int> R1;
  Bimap p;
  int g;
  bool deterministic;
  RBG rbg;

  /* Methods */

  void ReindexSubtColumn(unsigned k, unsigned c) {
    if (k == c) return;
    A.add_col(k, c);
    int R0k = R0[k];
    int R0c = R0[c];
    int Qkc = Q.entry(k, c);
    Q.add_rowcol(k, c);
    if (R0c) {
      Q.flip_submatrix({k, c});
    }
    R0[k] ^= R0c;
    R1[k] ^= R1[c] ^ Qkc ^ (R0k & R0c);
  }

  void MakePrincipal(unsigned c, unsigned j) {
    if (A.entry(j, c)) {
      const std::set<unsigned> H = A.cols_where_one(j);
      for (unsigned k : H) {
        if (k != c) {
          // This modifies A[j][k] but no other entries in A_j:
          ReindexSubtColumn(k, c);
        }
      }
      p.make_match(c, j);
    }
  }

  void ReselectPrincipalRow(
    unsigned c, std::optional<unsigned> j = std::nullopt)
  {
    std::optional<unsigned> n0;
    unsigned j0;
    for (unsigned j1 = 0; j1 < n; j1++) {
      if (!j || j1 != *j) {
        if (A.entry(j1, c)) {
          unsigned n1 = A.row_weight(j1);
          if (!n0 || n1 < *n0) {
            j0 = j1;
            n0 = n1;
          }
        }
      }
    }
    if (n0) {
      MakePrincipal(c, j0);
    }
  }

  std::optional<unsigned> principate(unsigned j) {
    std::optional<unsigned> c = p.inv_at(j);
    if (c) {
      ReselectPrincipalRow(*c, j);
      if (j != *p.fwd_at(*c)) {
        c = std::nullopt;
      }
    }
    return c;
  }

  // Swap column k with column r-1
  void ReindexSwapColumn(unsigned k) {
    const unsigned r1 = r - 1;
    if (k == r1) return;
    A.swap_col(k);
    std::iter_swap(R0.begin() + k, R0.begin() + r1);
    std::iter_swap(R1.begin() + k, R1.begin() + r1);
    Q.swap_rowcol(k);
    p.swap_fwd(k, r1);
  }

  void expand(unsigned j, const std::set<unsigned>& H) {
    A.zero_append_basis_col(j);
    Q.append_rowcol(H);
    r++;
  }

  void contract() {
    A.drop_final_col();
    Q.drop_final_rowcol();
    p.fwd_erase(r - 1);
    r--;
  }

  void FixFinalBit(int z) {
    if (z) {
      const unsigned r1 = r - 1;
      for (unsigned j = 0; j < n; j++) {
        b[j] ^= A.entry(j, r1);
      }
      for (unsigned h = 0; h < r1; h++) {
        R1[h] ^= Q.entry(h, r1);
      }
      if (z) {
        if (R0[r1]) g += 2;
        if (R1[r1]) g += 4;
        g %= 8;
      }
    }
    contract();
  }

  void ZeroColumnElim(unsigned c) {
    ReindexSwapColumn(c);
    const std::set<unsigned> H = Q.rows_with_terminal_1();
    int u0 = R0[r - 1];
    int u1 = R1[r - 1];
    contract();
    if (u0) {
      Q.flip_submatrix(H);
      for (unsigned h : H) {
        R0[h] ^= 1;
        R1[h] ^= R0[h] ^ u1;
      }
      g += 2;
      if (u0) g += 7;
      if (u1) g += 6;
      g %= 8;
    } else if (!H.empty()) {
      unsigned l = *H.begin();
      for (unsigned h : H) {
        ReindexSubtColumn(h, l);
      }
      ReindexSwapColumn(l);
      FixFinalBit(u1);
    }
  }

  void new_principal_column(
    unsigned j,
    int r0, int r1,
    std::optional<unsigned> c = std::nullopt,
    const std::set<unsigned>& H = std::set<unsigned>())
  {
    expand(j, H);
    b[j] = 0;
    R0[r - 1] = r0;
    R1[r - 1] = r1;
    p.make_match(r - 1, j);
    if (c) {
      ZeroColumnElim(*c);
    }
  }

  void SimulateX(unsigned j) { b[j] ^= 1; }

  void SimulateY(unsigned j) {
    g += 2; g %= 8;
    SimulateZ(j); SimulateX(j);
  }

  void SimulateZ(unsigned j) {
    if (b[j]) {
      g += 4; g %= 8;
    }
    const std::set<unsigned> H = A.cols_where_one(j);
    for (unsigned h : H) {
      R1[h] ^= 1;
    }
  }

  void SimulateH(unsigned j) {
    std::optional<unsigned> c = principate(j);
    const std::set<unsigned> H = A.cols_where_one(j);
    new_principal_column(j, 0, b[j], c, H);
  }

  void SimulateS(unsigned j) {
    const std::set<unsigned> H = A.cols_where_one(j);
    Q.flip_submatrix(H);
    int z = b[j];
    for (unsigned h : H) {
      R1[h] ^= R0[h] ^ z;
      R0[h] ^= 1;
    }
    if (z) {
      g += 2; g %= 8;
    }
  }

  void SimulateSdg(unsigned j) {
    const std::set<unsigned> H = A.cols_where_one(j);
    Q.flip_submatrix(H);
    const int z = b[j];
    for (unsigned h : H) {
      R0[h] ^= 1;
      R1[h] ^= R0[h] ^ z;
    }
    if (z) {
      g += 6; g %= 8;
    }
  }

  void SimulateCX(unsigned j, unsigned k) {
    A.add_row(k, j);
    b[k] ^= b[j];
    std::optional<unsigned> c = p.inv_at(k);
    if (c) {
      ReselectPrincipalRow(*c);
    }
  }

  void SimulateCZ(unsigned j, unsigned k) {
    const std::set<unsigned> H_j = A.cols_where_one(j);
    const std::set<unsigned> H_k = A.cols_where_one(k);
    Q.flip_submatrix(H_j, H_k);
    const std::set<unsigned> H_jk = A.cols_where_one(j, k);
    for (unsigned h : H_jk) {
      R1[h] ^= 1;
    }
    const int z_j = b[j];
    const int z_k = b[k];
    for (unsigned h : H_j) {
      R1[h] ^= z_k;
    }
    for (unsigned h : H_k) {
      R1[h] ^= z_j;
    }
    if (z_j && z_k) {
      g += 4; g %= 8;
    }
  }

  int toss_coin(std::optional<int> coin) {
    deterministic = false;
    if (coin) {
      return *coin;
    } else {
      return rbg.get();
    }
  }

  int SimulateMeasX(unsigned j, std::optional<int> coin = std::nullopt) {
    int beta;
    std::optional<unsigned> c = principate(j);
    if (c && Q.rowcol_is_zero(*c)) {
      if (R0[*c] == 0) {
        return R1[*c];
      } else {
        beta = toss_coin(coin);
        R0[*c] = 0;
        R1[*c] = beta;
        return beta;
      }
    } else {
      beta = toss_coin(coin);
    }
    const std::set<unsigned> H = A.cols_where_one(j);
    for (unsigned h : H) {
      R1[h] ^= beta;
    }
    new_principal_column(j, 0, beta, c);
    return beta;
  }

  int SimulateMeasY(unsigned j, std::optional<int> coin = std::nullopt) {
    int beta;
    std::optional<unsigned> c = principate(j);
    if (c && Q.rowcol_is_zero(*c)) {
      if (R0[*c] == 1) {
        return R1[*c] ^ b[j];
      } else {
        beta = toss_coin(coin);
        R0[*c] = 1;
        R1[*c] = beta;
        return beta;
      }
    } else {
      beta = toss_coin(coin);
    }
    const std::set<unsigned> H = A.cols_where_one(j);
    Q.flip_submatrix(H);
    const int z = b[j] ^ beta;
    for (unsigned h : H) {
      R0[h] ^= 1;
      R1[h] ^= R0[h] ^ z;
    }
    new_principal_column(j, 1, beta, c);
    return beta;
  }

  int SimulateMeasZ(unsigned j, std::optional<int> coin = std::nullopt) {
    if (A.row_weight(j) == 0) {
      return b[j];
    } else {
      int beta = toss_coin(coin);
      const std::set<unsigned> H = A.cols_where_one(j);
      unsigned k;
      unsigned m = n + 1;
      for (unsigned h : H) {
        unsigned c = A.col_weight(h);
        if (c < m) {
          k = h;
          m = c;
        }
      }
      ReindexSwapColumn(k);
      MakePrincipal(r - 1, j);
      FixFinalBit(beta ^ b[j]);
      return beta;
    }
  }

  void SimulateResetX(unsigned j) {
    if (SimulateMeasX(j, 0)) {
      SimulateZ(j);
    }
  }

  void SimulateResetY(unsigned j) {
    if (SimulateMeasY(j, 0)) {
      SimulateX(j);
    }
  }

  void SimulateResetZ(unsigned j) {
    if (SimulateMeasZ(j, 0)) {
      SimulateX(j);
    }
  }

  int phase() const { return g; }

  bool is_deterministic() const { return deterministic; }
};

/* Public interface */

Simplex::Simplex(unsigned n, int seed)
  : pImpl(std::make_unique<impl>(n, seed)) {}

Simplex::Simplex(const char *p, int seed)
  : pImpl(std::make_unique<impl>(p, seed)) {}

Simplex::~Simplex() = default;
Simplex::Simplex(const Simplex& other)
  : pImpl(std::make_unique<impl>(*other.pImpl)) {}
Simplex::Simplex(Simplex&& other) = default;
Simplex& Simplex::operator=(const Simplex& other) {
  return *this = Simplex(other);
}
Simplex& Simplex::operator=(Simplex&& other) = default;

unsigned Simplex::n() const { return pImpl->n; }
void Simplex::X(unsigned j) { pImpl->SimulateX(j); }
void Simplex::Y(unsigned j) { pImpl->SimulateY(j); }
void Simplex::Z(unsigned j) { pImpl->SimulateZ(j); }
void Simplex::H(unsigned j) { pImpl->SimulateH(j); }
void Simplex::S(unsigned j) { pImpl->SimulateS(j); }
void Simplex::Sdg(unsigned j) { pImpl->SimulateSdg(j); }
void Simplex::CX(unsigned j, unsigned k) { pImpl->SimulateCX(j, k); }
void Simplex::CZ(unsigned j, unsigned k) { pImpl->SimulateCZ(j, k); }
int Simplex::MeasX(unsigned j, std::optional<int> coin) {
  return pImpl->SimulateMeasX(j, coin);
}
int Simplex::MeasY(unsigned j, std::optional<int> coin) {
  return pImpl->SimulateMeasY(j, coin);
}
int Simplex::MeasZ(unsigned j, std::optional<int> coin) {
  return pImpl->SimulateMeasZ(j, coin);
}
void Simplex::ResetX(unsigned j) { pImpl->SimulateResetX(j); }
void Simplex::ResetY(unsigned j) { pImpl->SimulateResetY(j); }
void Simplex::ResetZ(unsigned j) { pImpl->SimulateResetZ(j); }
int Simplex::phase() const { return pImpl->phase(); }
bool Simplex::is_deterministic() const { return pImpl->is_deterministic(); }

std::ostream& operator<<(std::ostream& os, const Simplex& S) {
  os << "n: " << S.n() << std::endl;
  os << "A:" << std::endl << S.pImpl->A;
  os << "b: [ ";
  for (unsigned j = 0; j < S.n(); j++) {
    os << S.pImpl->b[j] << " ";
  }
  os << "]" << std::endl;
  os << "Q:" << std::endl;
  unsigned r = S.pImpl->Q.r();
  for (unsigned j = 0; j < r; j++) {
    os << "[ ";
    int R = S.pImpl->R0[j] + 2 * S.pImpl->R1[j];
    for (unsigned k = 0; k < r; k++) {
      os << ((j == k) ? R : S.pImpl->Q.entry(j, k)) << " ";
    }
    os << "]" << std::endl;
  }
  os << "g: " << S.pImpl->g << std::endl;
  os << "p: " << S.pImpl->p << std::endl;
  return os;
}
