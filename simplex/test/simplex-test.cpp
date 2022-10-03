#include <simplex.hpp>
#include <iostream>
#include <vector>

#define CHECK(a) \
  do { \
    if (!(a)) { \
      std::cout << "Assertion '" << #a << "' (line " << __LINE__ << ") failed!" << std::endl; \
      return 1; \
    } \
  } while (0)

#define CHECK_OK(a) CHECK(!(a))

static int test_X() {
  Simplex S(2);
  S.X(0);
  int b0 = S.MeasZ(0);
  int b1 = S.MeasZ(1);
  CHECK(b0 == 1);
  CHECK(b1 == 0);
  CHECK(S.is_deterministic());
  return 0;
}

static int test_Y() {
  Simplex S(2);
  S.Y(0);
  int b0 = S.MeasZ(0);
  int b1 = S.MeasZ(1);
  CHECK(b0 == 1);
  CHECK(b1 == 0);
  CHECK(S.is_deterministic());
  return 0;
}

static int test_Z() {
  Simplex S(2);
  S.Z(0);
  int b0 = S.MeasZ(0);
  int b1 = S.MeasZ(1);
  CHECK(b0 == 0);
  CHECK(b1 == 0);
  CHECK(S.is_deterministic());
  return 0;
}

static int test_H() {
  Simplex S(2);
  S.H(1);
  int b0 = S.MeasZ(0);
  CHECK(b0 == 0);
  CHECK(S.is_deterministic());
  int b1 = S.MeasZ(1);
  CHECK(!S.is_deterministic());
  return 0;
}

static int test_S() {
  Simplex S(1);
  S.S(0);
  int b0 = S.MeasZ(0);
  CHECK(b0 == 0);
  CHECK(S.is_deterministic());
  return 0;
}

static int test_Sdg() {
  Simplex S(1);
  S.Sdg(0);
  int b0 = S.MeasZ(0);
  CHECK(b0 == 0);
  CHECK(S.is_deterministic());
  return 0;
}

static int test_CX() {
  Simplex S(2);
  S.X(0);
  S.CX(0, 1);
  int b0 = S.MeasZ(0);
  int b1 = S.MeasZ(1);
  CHECK(b0 == 1);
  CHECK(b1 == 1);
  CHECK(S.is_deterministic());
  return 0;
}

static int test_CZ() {
  Simplex S(2);
  S.X(0);
  S.CZ(0, 1);
  int b0 = S.MeasZ(0);
  int b1 = S.MeasZ(1);
  CHECK(b0 == 1);
  CHECK(b1 == 0);
  CHECK(S.is_deterministic());
  return 0;
}

static int test_MeasX() {
  Simplex S(1);
  S.H(0);
  int b0 = S.MeasX(0);
  CHECK(b0 == 0);
  CHECK(S.is_deterministic());
  return 0;
}

static int test_MeasY() {
  Simplex S(1);
  S.H(0);
  S.S(0);
  int b0 = S.MeasY(0);
  CHECK(b0 == 0);
  CHECK(S.is_deterministic());
  return 0;
}

static int test_MeasZ() {
  Simplex S(2);
  S.H(0);
  S.CX(0, 1);
  int b0 = S.MeasZ(0);
  int b1 = S.MeasZ(1);
  CHECK(b0 == b1);
  CHECK(!S.is_deterministic());
  return 0;
}

static int test_many_qubits() {
  Simplex S(30);
  S.H(0);
  for (unsigned i = 1; i < 30; i++) {
    S.CX(0, i);
  }
  int b0 = S.MeasZ(0);
  for (unsigned i = 1; i < 30; i++) {
    int b = S.MeasZ(i);
    CHECK(b == b0);
  }
  CHECK(!S.is_deterministic());
  return 0;
}

static int test_complicated_circuit_1() {
  Simplex S(3);
  S.H(0);
  S.H(1);
  S.H(2);
  S.CX(0,1);
  S.CZ(1,2);
  S.S(1);
  S.Y(2);
  S.CX(2,0);
  S.Z(1);
  S.MeasX(1);
  S.MeasZ(0);
  S.MeasY(2);
  S.S(0);
  S.S(1);
  S.S(2);
  S.CX(1,2);
  S.CX(0,1);
  S.H(0);
  S.X(1);
  S.Z(1);
  S.S(0);
  S.CX(0,2);
  S.CZ(1,0);
  S.MeasZ(0);
  S.MeasZ(1);
  S.MeasY(1);
  S.MeasZ(2);
  S.H(0);
  S.H(1);
  S.H(2);
  S.CX(0,1);
  S.CZ(1,2);
  S.S(1);
  S.Y(2);
  S.CX(2,0);
  S.Z(1);
  S.MeasX(1);
  return 0;
}

static int test_complicated_circuit_2() {
  Simplex S(3);
  S.H(0);
  S.H(1);
  S.H(2);
  S.CX(0,1);
  S.CZ(1,2);
  S.S(1);
  S.Y(2);
  S.CX(2,0);
  S.Z(1);
  S.S(0);
  S.S(1);
  S.S(2);
  S.CX(1,2);
  S.CX(0,1);
  S.H(0);
  return 0;
}

static int test_complicated_x_cx_circuit() {
  Simplex S(4);
  S.X(0);
  S.X(2);
  S.CX(0,3);
  S.CX(3,1);
  S.X(1);
  S.CX(0,1);
  S.X(3);
  S.CX(3,2);
  S.X(2);
  S.CX(1,2);
  S.X(1);
  S.CX(1,2);
  S.CX(1,3);
  int b0 = S.MeasZ(0);
  int b1 = S.MeasZ(1);
  int b2 = S.MeasZ(2);
  int b3 = S.MeasZ(3);
  CHECK(b0 == 1);
  CHECK(b1 == 0);
  CHECK(b2 == 1);
  CHECK(b3 == 0);
  CHECK(S.is_deterministic());
  return 0;
}

static int test_invert() {
  Simplex S(3);
  S.H(0);
  S.CX(0,1);
  S.S(1);
  S.Y(2);
  S.CZ(1,2);
  S.H(2);
  S.CX(2,1);
  S.Z(0);
  S.X(1);
  S.CZ(0,1);
  S.CZ(0,1);
  S.X(1);
  S.Z(0);
  S.CX(2,1);
  S.H(2);
  S.CZ(1,2);
  S.Y(2);
  S.Sdg(1);
  S.CX(0,1);
  S.H(0);
  int b0 = S.MeasZ(0);
  int b1 = S.MeasZ(1);
  int b2 = S.MeasZ(2);
  CHECK(b0 == 0);
  CHECK(b1 == 0);
  CHECK(b2 == 0);
  return 0;
}

static int test_measurements() {
  // Succeeding measurement in the same basis yields the same result.
  const  std::vector<std::pair<int, int>> twocoins = {
    {0,0}, {0,1}, {1,0}, {1,1}};
  // X,X
  for (auto coins : twocoins) {
    Simplex S(1);
    int b0 = S.MeasX(0, coins.first);
    int b1 = S.MeasX(0, coins.second);
    CHECK(b0 == b1);
  }
  // Y,Y
  for (auto coins : twocoins) {
    Simplex S(1);
    int b0 = S.MeasY(0, coins.first);
    int b1 = S.MeasY(0, coins.second);
    CHECK(b0 == b1);
  }
  // Z,Z
  for (auto coins : twocoins) {
    Simplex S(1);
    S.H(0);
    int b0 = S.MeasX(0, coins.first);
    int b1 = S.MeasX(0, coins.second);
    CHECK(b0 == b1);
  }
  // Succeeding orthogonal measurement yields 0 or 1 with equal probability.
  // X,Y
  for (int coin = 0; coin < 2; coin++) {
    Simplex S(1);
    S.MeasX(0, coin);
    Simplex T(S);
    int b0 = S.MeasY(0, 0);
    int b1 = T.MeasY(0, 1);
    CHECK(b0 != b1);
  }
  // X,Z
  for (int coin = 0; coin < 2; coin++) {
    Simplex S(1);
    S.MeasX(0, coin);
    Simplex T(S);
    int b0 = S.MeasZ(0, 0);
    int b1 = T.MeasZ(0, 1);
    CHECK(b0 != b1);
  }
  // Y,Z
  for (int coin = 0; coin < 2; coin++) {
    Simplex S(1);
    S.MeasY(0, coin);
    Simplex T(S);
    int b0 = S.MeasZ(0, 0);
    int b1 = T.MeasZ(0, 1);
    CHECK(b0 != b1);
  }
  // Y,X
  for (int coin = 0; coin < 2; coin++) {
    Simplex S(1);
    S.MeasY(0, coin);
    Simplex T(S);
    int b0 = S.MeasX(0, 0);
    int b1 = T.MeasX(0, 1);
    CHECK(b0 != b1);
  }
  // Z,X
  for (int coin = 0; coin < 2; coin++) {
    Simplex S(1);
    S.H(0);
    S.MeasZ(0, coin);
    Simplex T(S);
    int b0 = S.MeasX(0, 0);
    int b1 = T.MeasX(0, 1);
    CHECK(b0 != b1);
  }
  // Z,Y
  for (int coin = 0; coin < 2; coin++) {
    Simplex S(1);
    S.H(0);
    S.MeasZ(0, coin);
    Simplex T(S);
    int b0 = S.MeasY(0, 0);
    int b1 = T.MeasY(0, 1);
    CHECK(b0 != b1);
  }
  return 0;
}

static int test_mid_circ_meas() {
  { // [X; MZ; X; MZ] should yield 1 then 0
    Simplex S(1);
    S.X(0);
    CHECK(S.MeasZ(0) == 1);
    S.X(0);
    CHECK(S.MeasZ(0) == 0);
  }
  { // Replace first MZ with [H; MX; H]
    Simplex S(1);
    S.X(0);
    S.H(0);
    CHECK(S.MeasX(0) == 1);
    S.H(0);
    S.X(0);
    CHECK(S.MeasZ(0) == 0);
  }
  { // Replace first MZ with [H; S; MY; S*; H]
    Simplex S(1);
    S.X(0);
    S.H(0);
    S.S(0);
    CHECK(S.MeasY(0) == 1);
    S.Sdg(0);
    S.H(0);
    S.X(0);
    CHECK(S.MeasZ(0) == 0);
  }
  { // A more complicated circuit that ends in the (1,1) state
    Simplex S(2);
    S.X(0);
    S.CX(0, 1);
    S.S(1);
    S.CX(1, 0);
    S.CZ(0, 1);
    S.CX(1, 0);
    S.S(1);
    S.Z(0);
    {
      Simplex T(S);
      CHECK(T.MeasZ(0) == 1);
      CHECK(T.MeasZ(1) == 1);
      T.Z(0);
      T.Sdg(1);
      T.CX(1, 0);
      T.CZ(0, 1);
      T.CX(1, 0);
      T.Sdg(1);
      T.CX(0, 1);
      T.X(0);
      CHECK(T.MeasZ(0) == 0);
      CHECK(T.MeasZ(1) == 0);
    }
    { // Replace each first MZ with [H; MX; H]
      Simplex T(S);
      T.H(0);
      T.H(1);
      CHECK(T.MeasX(0) == 1);
      CHECK(T.MeasX(1) == 1);
      T.H(0);
      T.H(1);
      T.Z(0);
      T.Sdg(1);
      T.CX(1, 0);
      T.CZ(0, 1);
      T.CX(1, 0);
      T.Sdg(1);
      T.CX(0, 1);
      T.X(0);
      CHECK(T.MeasZ(0) == 0);
      CHECK(T.MeasZ(1) == 0);
    }
    { // Replace each first MZ with [H; S; MY; S*; H]
      Simplex T(S);
      T.H(0);
      T.S(0);
      T.H(1);
      T.S(1);
      CHECK(T.MeasY(0) == 1);
      CHECK(T.MeasY(1) == 1);
      T.Sdg(0);
      T.H(0);
      T.Sdg(1);
      T.H(1);
      T.Z(0);
      T.Sdg(1);
      T.CX(1, 0);
      T.CZ(0, 1);
      T.CX(1, 0);
      T.Sdg(1);
      T.CX(0, 1);
      T.X(0);
      CHECK(T.MeasZ(0) == 0);
      CHECK(T.MeasZ(1) == 0);
    }
  }
  return 0;
}

// [BEGIN] Autogenerated by `gen-test-cases` script

static int test_dist_0_0() {
  Simplex S(2);
  S.H(0);
  S.CX(0, 1);
  std::vector<unsigned> v(4, 0);
  for (unsigned m = 0; m < 4; m++) {
    Simplex T(S);
    unsigned c = 0;
    c += T.MeasX(0, (m >> 0) & 1) << 0;
    c += T.MeasX(1, (m >> 1) & 1) << 1;
    v[c]++;
  }
  CHECK(v[0] == 2);
  CHECK(v[1] == 0);
  CHECK(v[2] == 0);
  CHECK(v[3] == 2);
  return 0;
}
static int test_dist_0_1() {
  Simplex S(2);
  S.H(0);
  S.CX(0, 1);
  std::vector<unsigned> v(4, 0);
  for (unsigned m = 0; m < 4; m++) {
    Simplex T(S);
    unsigned c = 0;
    c += T.MeasY(0, (m >> 0) & 1) << 0;
    c += T.MeasY(1, (m >> 1) & 1) << 1;
    v[c]++;
  }
  CHECK(v[0] == 0);
  CHECK(v[1] == 2);
  CHECK(v[2] == 2);
  CHECK(v[3] == 0);
  return 0;
}
static int test_dist_0_2() {
  Simplex S(2);
  S.H(0);
  S.CX(0, 1);
  std::vector<unsigned> v(4, 0);
  for (unsigned m = 0; m < 4; m++) {
    Simplex T(S);
    unsigned c = 0;
    c += T.MeasZ(0, (m >> 0) & 1) << 0;
    c += T.MeasZ(1, (m >> 1) & 1) << 1;
    v[c]++;
  }
  CHECK(v[0] == 2);
  CHECK(v[1] == 0);
  CHECK(v[2] == 0);
  CHECK(v[3] == 2);
  return 0;
}
static int test_dist_0_3() {
  Simplex S(2);
  S.H(0);
  S.CX(0, 1);
  std::vector<unsigned> v(4, 0);
  for (unsigned m = 0; m < 4; m++) {
    Simplex T(S);
    unsigned c = 0;
    c += T.MeasX(0, (m >> 0) & 1) << 0;
    c += T.MeasY(1, (m >> 1) & 1) << 1;
    v[c]++;
  }
  CHECK(v[0] == 1);
  CHECK(v[1] == 1);
  CHECK(v[2] == 1);
  CHECK(v[3] == 1);
  return 0;
}
static int test_dist_0_4() {
  Simplex S(2);
  S.H(0);
  S.CX(0, 1);
  std::vector<unsigned> v(4, 0);
  for (unsigned m = 0; m < 4; m++) {
    Simplex T(S);
    unsigned c = 0;
    c += T.MeasY(0, (m >> 0) & 1) << 0;
    c += T.MeasZ(1, (m >> 1) & 1) << 1;
    v[c]++;
  }
  CHECK(v[0] == 1);
  CHECK(v[1] == 1);
  CHECK(v[2] == 1);
  CHECK(v[3] == 1);
  return 0;
}
static int test_dist_0_5() {
  Simplex S(2);
  S.H(0);
  S.CX(0, 1);
  std::vector<unsigned> v(4, 0);
  for (unsigned m = 0; m < 4; m++) {
    Simplex T(S);
    unsigned c = 0;
    c += T.MeasZ(0, (m >> 0) & 1) << 0;
    c += T.MeasX(1, (m >> 1) & 1) << 1;
    v[c]++;
  }
  CHECK(v[0] == 1);
  CHECK(v[1] == 1);
  CHECK(v[2] == 1);
  CHECK(v[3] == 1);
  return 0;
}
static int test_dist_1_0() {
  Simplex S(3);
  S.H(0);
  S.CX(0, 1);
  S.CX(0, 2);
  std::vector<unsigned> v(8, 0);
  for (unsigned m = 0; m < 8; m++) {
    Simplex T(S);
    unsigned c = 0;
    c += T.MeasX(0, (m >> 0) & 1) << 0;
    c += T.MeasX(1, (m >> 1) & 1) << 1;
    c += T.MeasX(2, (m >> 2) & 1) << 2;
    v[c]++;
  }
  CHECK(v[0] == 2);
  CHECK(v[1] == 0);
  CHECK(v[2] == 0);
  CHECK(v[3] == 2);
  CHECK(v[4] == 0);
  CHECK(v[5] == 2);
  CHECK(v[6] == 2);
  CHECK(v[7] == 0);
  return 0;
}
static int test_dist_1_1() {
  Simplex S(3);
  S.H(0);
  S.CX(0, 1);
  S.CX(0, 2);
  std::vector<unsigned> v(8, 0);
  for (unsigned m = 0; m < 8; m++) {
    Simplex T(S);
    unsigned c = 0;
    c += T.MeasY(0, (m >> 0) & 1) << 0;
    c += T.MeasY(1, (m >> 1) & 1) << 1;
    c += T.MeasY(2, (m >> 2) & 1) << 2;
    v[c]++;
  }
  CHECK(v[0] == 1);
  CHECK(v[1] == 1);
  CHECK(v[2] == 1);
  CHECK(v[3] == 1);
  CHECK(v[4] == 1);
  CHECK(v[5] == 1);
  CHECK(v[6] == 1);
  CHECK(v[7] == 1);
  return 0;
}
static int test_dist_1_2() {
  Simplex S(3);
  S.H(0);
  S.CX(0, 1);
  S.CX(0, 2);
  std::vector<unsigned> v(8, 0);
  for (unsigned m = 0; m < 8; m++) {
    Simplex T(S);
    unsigned c = 0;
    c += T.MeasZ(0, (m >> 0) & 1) << 0;
    c += T.MeasZ(1, (m >> 1) & 1) << 1;
    c += T.MeasZ(2, (m >> 2) & 1) << 2;
    v[c]++;
  }
  CHECK(v[0] == 4);
  CHECK(v[1] == 0);
  CHECK(v[2] == 0);
  CHECK(v[3] == 0);
  CHECK(v[4] == 0);
  CHECK(v[5] == 0);
  CHECK(v[6] == 0);
  CHECK(v[7] == 4);
  return 0;
}
static int test_dist_1_3() {
  Simplex S(3);
  S.H(0);
  S.CX(0, 1);
  S.CX(0, 2);
  std::vector<unsigned> v(8, 0);
  for (unsigned m = 0; m < 8; m++) {
    Simplex T(S);
    unsigned c = 0;
    c += T.MeasX(0, (m >> 0) & 1) << 0;
    c += T.MeasY(1, (m >> 1) & 1) << 1;
    c += T.MeasZ(2, (m >> 2) & 1) << 2;
    v[c]++;
  }
  CHECK(v[0] == 1);
  CHECK(v[1] == 1);
  CHECK(v[2] == 1);
  CHECK(v[3] == 1);
  CHECK(v[4] == 1);
  CHECK(v[5] == 1);
  CHECK(v[6] == 1);
  CHECK(v[7] == 1);
  return 0;
}
static int test_dist_1_4() {
  Simplex S(3);
  S.H(0);
  S.CX(0, 1);
  S.CX(0, 2);
  std::vector<unsigned> v(8, 0);
  for (unsigned m = 0; m < 8; m++) {
    Simplex T(S);
    unsigned c = 0;
    c += T.MeasY(0, (m >> 0) & 1) << 0;
    c += T.MeasZ(1, (m >> 1) & 1) << 1;
    c += T.MeasX(2, (m >> 2) & 1) << 2;
    v[c]++;
  }
  CHECK(v[0] == 1);
  CHECK(v[1] == 1);
  CHECK(v[2] == 1);
  CHECK(v[3] == 1);
  CHECK(v[4] == 1);
  CHECK(v[5] == 1);
  CHECK(v[6] == 1);
  CHECK(v[7] == 1);
  return 0;
}
static int test_dist_1_5() {
  Simplex S(3);
  S.H(0);
  S.CX(0, 1);
  S.CX(0, 2);
  std::vector<unsigned> v(8, 0);
  for (unsigned m = 0; m < 8; m++) {
    Simplex T(S);
    unsigned c = 0;
    c += T.MeasZ(0, (m >> 0) & 1) << 0;
    c += T.MeasX(1, (m >> 1) & 1) << 1;
    c += T.MeasY(2, (m >> 2) & 1) << 2;
    v[c]++;
  }
  CHECK(v[0] == 1);
  CHECK(v[1] == 1);
  CHECK(v[2] == 1);
  CHECK(v[3] == 1);
  CHECK(v[4] == 1);
  CHECK(v[5] == 1);
  CHECK(v[6] == 1);
  CHECK(v[7] == 1);
  return 0;
}
static int test_dist_2_0() {
  Simplex S(3);
  S.H(1);
  S.CX(1, 2);
  S.CZ(2, 0);
  S.S(0);
  S.CX(0, 2);
  S.H(1);
  S.Sdg(1);
  S.S(0);
  S.CX(1, 0);
  std::vector<unsigned> v(8, 0);
  for (unsigned m = 0; m < 8; m++) {
    Simplex T(S);
    unsigned c = 0;
    c += T.MeasX(0, (m >> 0) & 1) << 0;
    c += T.MeasX(1, (m >> 1) & 1) << 1;
    c += T.MeasX(2, (m >> 2) & 1) << 2;
    v[c]++;
  }
  CHECK(v[0] == 1);
  CHECK(v[1] == 1);
  CHECK(v[2] == 1);
  CHECK(v[3] == 1);
  CHECK(v[4] == 1);
  CHECK(v[5] == 1);
  CHECK(v[6] == 1);
  CHECK(v[7] == 1);
  return 0;
}
static int test_dist_2_1() {
  Simplex S(3);
  S.H(1);
  S.CX(1, 2);
  S.CZ(2, 0);
  S.S(0);
  S.CX(0, 2);
  S.H(1);
  S.Sdg(1);
  S.S(0);
  S.CX(1, 0);
  std::vector<unsigned> v(8, 0);
  for (unsigned m = 0; m < 8; m++) {
    Simplex T(S);
    unsigned c = 0;
    c += T.MeasY(0, (m >> 0) & 1) << 0;
    c += T.MeasY(1, (m >> 1) & 1) << 1;
    c += T.MeasY(2, (m >> 2) & 1) << 2;
    v[c]++;
  }
  CHECK(v[0] == 0);
  CHECK(v[1] == 2);
  CHECK(v[2] == 2);
  CHECK(v[3] == 0);
  CHECK(v[4] == 2);
  CHECK(v[5] == 0);
  CHECK(v[6] == 0);
  CHECK(v[7] == 2);
  return 0;
}
static int test_dist_2_2() {
  Simplex S(3);
  S.H(1);
  S.CX(1, 2);
  S.CZ(2, 0);
  S.S(0);
  S.CX(0, 2);
  S.H(1);
  S.Sdg(1);
  S.S(0);
  S.CX(1, 0);
  std::vector<unsigned> v(8, 0);
  for (unsigned m = 0; m < 8; m++) {
    Simplex T(S);
    unsigned c = 0;
    c += T.MeasZ(0, (m >> 0) & 1) << 0;
    c += T.MeasZ(1, (m >> 1) & 1) << 1;
    c += T.MeasZ(2, (m >> 2) & 1) << 2;
    v[c]++;
  }
  CHECK(v[0] == 2);
  CHECK(v[1] == 0);
  CHECK(v[2] == 0);
  CHECK(v[3] == 2);
  CHECK(v[4] == 2);
  CHECK(v[5] == 0);
  CHECK(v[6] == 0);
  CHECK(v[7] == 2);
  return 0;
}
static int test_dist_2_3() {
  Simplex S(3);
  S.H(1);
  S.CX(1, 2);
  S.CZ(2, 0);
  S.S(0);
  S.CX(0, 2);
  S.H(1);
  S.Sdg(1);
  S.S(0);
  S.CX(1, 0);
  std::vector<unsigned> v(8, 0);
  for (unsigned m = 0; m < 8; m++) {
    Simplex T(S);
    unsigned c = 0;
    c += T.MeasX(0, (m >> 0) & 1) << 0;
    c += T.MeasY(1, (m >> 1) & 1) << 1;
    c += T.MeasZ(2, (m >> 2) & 1) << 2;
    v[c]++;
  }
  CHECK(v[0] == 0);
  CHECK(v[1] == 2);
  CHECK(v[2] == 2);
  CHECK(v[3] == 0);
  CHECK(v[4] == 2);
  CHECK(v[5] == 0);
  CHECK(v[6] == 0);
  CHECK(v[7] == 2);
  return 0;
}
static int test_dist_2_4() {
  Simplex S(3);
  S.H(1);
  S.CX(1, 2);
  S.CZ(2, 0);
  S.S(0);
  S.CX(0, 2);
  S.H(1);
  S.Sdg(1);
  S.S(0);
  S.CX(1, 0);
  std::vector<unsigned> v(8, 0);
  for (unsigned m = 0; m < 8; m++) {
    Simplex T(S);
    unsigned c = 0;
    c += T.MeasY(0, (m >> 0) & 1) << 0;
    c += T.MeasZ(1, (m >> 1) & 1) << 1;
    c += T.MeasX(2, (m >> 2) & 1) << 2;
    v[c]++;
  }
  CHECK(v[0] == 2);
  CHECK(v[1] == 2);
  CHECK(v[2] == 0);
  CHECK(v[3] == 0);
  CHECK(v[4] == 0);
  CHECK(v[5] == 0);
  CHECK(v[6] == 2);
  CHECK(v[7] == 2);
  return 0;
}
static int test_dist_2_5() {
  Simplex S(3);
  S.H(1);
  S.CX(1, 2);
  S.CZ(2, 0);
  S.S(0);
  S.CX(0, 2);
  S.H(1);
  S.Sdg(1);
  S.S(0);
  S.CX(1, 0);
  std::vector<unsigned> v(8, 0);
  for (unsigned m = 0; m < 8; m++) {
    Simplex T(S);
    unsigned c = 0;
    c += T.MeasZ(0, (m >> 0) & 1) << 0;
    c += T.MeasX(1, (m >> 1) & 1) << 1;
    c += T.MeasY(2, (m >> 2) & 1) << 2;
    v[c]++;
  }
  CHECK(v[0] == 1);
  CHECK(v[1] == 1);
  CHECK(v[2] == 1);
  CHECK(v[3] == 1);
  CHECK(v[4] == 1);
  CHECK(v[5] == 1);
  CHECK(v[6] == 1);
  CHECK(v[7] == 1);
  return 0;
}
static int test_dist_3_0() {
  Simplex S(3);
  S.H(0);
  S.H(1);
  S.H(2);
  std::vector<unsigned> v(8, 0);
  for (unsigned m = 0; m < 8; m++) {
    Simplex T(S);
    unsigned c = 0;
    c += T.MeasX(0, (m >> 0) & 1) << 0;
    c += T.MeasX(1, (m >> 1) & 1) << 1;
    c += T.MeasX(2, (m >> 2) & 1) << 2;
    v[c]++;
  }
  CHECK(v[0] == 8);
  CHECK(v[1] == 0);
  CHECK(v[2] == 0);
  CHECK(v[3] == 0);
  CHECK(v[4] == 0);
  CHECK(v[5] == 0);
  CHECK(v[6] == 0);
  CHECK(v[7] == 0);
  return 0;
}
static int test_dist_3_1() {
  Simplex S(3);
  S.H(0);
  S.H(1);
  S.H(2);
  std::vector<unsigned> v(8, 0);
  for (unsigned m = 0; m < 8; m++) {
    Simplex T(S);
    unsigned c = 0;
    c += T.MeasY(0, (m >> 0) & 1) << 0;
    c += T.MeasY(1, (m >> 1) & 1) << 1;
    c += T.MeasY(2, (m >> 2) & 1) << 2;
    v[c]++;
  }
  CHECK(v[0] == 1);
  CHECK(v[1] == 1);
  CHECK(v[2] == 1);
  CHECK(v[3] == 1);
  CHECK(v[4] == 1);
  CHECK(v[5] == 1);
  CHECK(v[6] == 1);
  CHECK(v[7] == 1);
  return 0;
}
static int test_dist_3_2() {
  Simplex S(3);
  S.H(0);
  S.H(1);
  S.H(2);
  std::vector<unsigned> v(8, 0);
  for (unsigned m = 0; m < 8; m++) {
    Simplex T(S);
    unsigned c = 0;
    c += T.MeasZ(0, (m >> 0) & 1) << 0;
    c += T.MeasZ(1, (m >> 1) & 1) << 1;
    c += T.MeasZ(2, (m >> 2) & 1) << 2;
    v[c]++;
  }
  CHECK(v[0] == 1);
  CHECK(v[1] == 1);
  CHECK(v[2] == 1);
  CHECK(v[3] == 1);
  CHECK(v[4] == 1);
  CHECK(v[5] == 1);
  CHECK(v[6] == 1);
  CHECK(v[7] == 1);
  return 0;
}
static int test_dist_3_3() {
  Simplex S(3);
  S.H(0);
  S.H(1);
  S.H(2);
  std::vector<unsigned> v(8, 0);
  for (unsigned m = 0; m < 8; m++) {
    Simplex T(S);
    unsigned c = 0;
    c += T.MeasX(0, (m >> 0) & 1) << 0;
    c += T.MeasY(1, (m >> 1) & 1) << 1;
    c += T.MeasZ(2, (m >> 2) & 1) << 2;
    v[c]++;
  }
  CHECK(v[0] == 2);
  CHECK(v[1] == 0);
  CHECK(v[2] == 2);
  CHECK(v[3] == 0);
  CHECK(v[4] == 2);
  CHECK(v[5] == 0);
  CHECK(v[6] == 2);
  CHECK(v[7] == 0);
  return 0;
}
static int test_dist_3_4() {
  Simplex S(3);
  S.H(0);
  S.H(1);
  S.H(2);
  std::vector<unsigned> v(8, 0);
  for (unsigned m = 0; m < 8; m++) {
    Simplex T(S);
    unsigned c = 0;
    c += T.MeasY(0, (m >> 0) & 1) << 0;
    c += T.MeasZ(1, (m >> 1) & 1) << 1;
    c += T.MeasX(2, (m >> 2) & 1) << 2;
    v[c]++;
  }
  CHECK(v[0] == 2);
  CHECK(v[1] == 2);
  CHECK(v[2] == 2);
  CHECK(v[3] == 2);
  CHECK(v[4] == 0);
  CHECK(v[5] == 0);
  CHECK(v[6] == 0);
  CHECK(v[7] == 0);
  return 0;
}
static int test_dist_3_5() {
  Simplex S(3);
  S.H(0);
  S.H(1);
  S.H(2);
  std::vector<unsigned> v(8, 0);
  for (unsigned m = 0; m < 8; m++) {
    Simplex T(S);
    unsigned c = 0;
    c += T.MeasZ(0, (m >> 0) & 1) << 0;
    c += T.MeasX(1, (m >> 1) & 1) << 1;
    c += T.MeasY(2, (m >> 2) & 1) << 2;
    v[c]++;
  }
  CHECK(v[0] == 2);
  CHECK(v[1] == 2);
  CHECK(v[2] == 0);
  CHECK(v[3] == 0);
  CHECK(v[4] == 2);
  CHECK(v[5] == 2);
  CHECK(v[6] == 0);
  CHECK(v[7] == 0);
  return 0;
}
static int test_dist_4_0() {
  Simplex S(2);
  S.H(1);
  S.X(0);
  S.CX(1, 0);
  std::vector<unsigned> v(4, 0);
  for (unsigned m = 0; m < 4; m++) {
    Simplex T(S);
    unsigned c = 0;
    c += T.MeasX(0, (m >> 0) & 1) << 0;
    c += T.MeasX(1, (m >> 1) & 1) << 1;
    v[c]++;
  }
  CHECK(v[0] == 2);
  CHECK(v[1] == 0);
  CHECK(v[2] == 0);
  CHECK(v[3] == 2);
  return 0;
}
static int test_dist_4_1() {
  Simplex S(2);
  S.H(1);
  S.X(0);
  S.CX(1, 0);
  std::vector<unsigned> v(4, 0);
  for (unsigned m = 0; m < 4; m++) {
    Simplex T(S);
    unsigned c = 0;
    c += T.MeasY(0, (m >> 0) & 1) << 0;
    c += T.MeasY(1, (m >> 1) & 1) << 1;
    v[c]++;
  }
  CHECK(v[0] == 2);
  CHECK(v[1] == 0);
  CHECK(v[2] == 0);
  CHECK(v[3] == 2);
  return 0;
}
static int test_dist_4_2() {
  Simplex S(2);
  S.H(1);
  S.X(0);
  S.CX(1, 0);
  std::vector<unsigned> v(4, 0);
  for (unsigned m = 0; m < 4; m++) {
    Simplex T(S);
    unsigned c = 0;
    c += T.MeasZ(0, (m >> 0) & 1) << 0;
    c += T.MeasZ(1, (m >> 1) & 1) << 1;
    v[c]++;
  }
  CHECK(v[0] == 0);
  CHECK(v[1] == 2);
  CHECK(v[2] == 2);
  CHECK(v[3] == 0);
  return 0;
}
static int test_dist_4_3() {
  Simplex S(2);
  S.H(1);
  S.X(0);
  S.CX(1, 0);
  std::vector<unsigned> v(4, 0);
  for (unsigned m = 0; m < 4; m++) {
    Simplex T(S);
    unsigned c = 0;
    c += T.MeasX(0, (m >> 0) & 1) << 0;
    c += T.MeasY(1, (m >> 1) & 1) << 1;
    v[c]++;
  }
  CHECK(v[0] == 1);
  CHECK(v[1] == 1);
  CHECK(v[2] == 1);
  CHECK(v[3] == 1);
  return 0;
}
static int test_dist_4_4() {
  Simplex S(2);
  S.H(1);
  S.X(0);
  S.CX(1, 0);
  std::vector<unsigned> v(4, 0);
  for (unsigned m = 0; m < 4; m++) {
    Simplex T(S);
    unsigned c = 0;
    c += T.MeasY(0, (m >> 0) & 1) << 0;
    c += T.MeasZ(1, (m >> 1) & 1) << 1;
    v[c]++;
  }
  CHECK(v[0] == 1);
  CHECK(v[1] == 1);
  CHECK(v[2] == 1);
  CHECK(v[3] == 1);
  return 0;
}
static int test_dist_4_5() {
  Simplex S(2);
  S.H(1);
  S.X(0);
  S.CX(1, 0);
  std::vector<unsigned> v(4, 0);
  for (unsigned m = 0; m < 4; m++) {
    Simplex T(S);
    unsigned c = 0;
    c += T.MeasZ(0, (m >> 0) & 1) << 0;
    c += T.MeasX(1, (m >> 1) & 1) << 1;
    v[c]++;
  }
  CHECK(v[0] == 1);
  CHECK(v[1] == 1);
  CHECK(v[2] == 1);
  CHECK(v[3] == 1);
  return 0;
}
static int test_dist_5_0() {
  Simplex S(3);
  S.Y(0);
  S.H(1);
  S.CX(1, 2);
  S.Sdg(2);
  S.CZ(0, 2);
  S.Z(0);
  S.S(1);
  S.CX(1, 0);
  std::vector<unsigned> v(8, 0);
  for (unsigned m = 0; m < 8; m++) {
    Simplex T(S);
    unsigned c = 0;
    c += T.MeasX(0, (m >> 0) & 1) << 0;
    c += T.MeasX(1, (m >> 1) & 1) << 1;
    c += T.MeasX(2, (m >> 2) & 1) << 2;
    v[c]++;
  }
  CHECK(v[0] == 0);
  CHECK(v[1] == 2);
  CHECK(v[2] == 2);
  CHECK(v[3] == 0);
  CHECK(v[4] == 2);
  CHECK(v[5] == 0);
  CHECK(v[6] == 0);
  CHECK(v[7] == 2);
  return 0;
}
static int test_dist_5_1() {
  Simplex S(3);
  S.Y(0);
  S.H(1);
  S.CX(1, 2);
  S.Sdg(2);
  S.CZ(0, 2);
  S.Z(0);
  S.S(1);
  S.CX(1, 0);
  std::vector<unsigned> v(8, 0);
  for (unsigned m = 0; m < 8; m++) {
    Simplex T(S);
    unsigned c = 0;
    c += T.MeasY(0, (m >> 0) & 1) << 0;
    c += T.MeasY(1, (m >> 1) & 1) << 1;
    c += T.MeasY(2, (m >> 2) & 1) << 2;
    v[c]++;
  }
  CHECK(v[0] == 1);
  CHECK(v[1] == 1);
  CHECK(v[2] == 1);
  CHECK(v[3] == 1);
  CHECK(v[4] == 1);
  CHECK(v[5] == 1);
  CHECK(v[6] == 1);
  CHECK(v[7] == 1);
  return 0;
}
static int test_dist_5_2() {
  Simplex S(3);
  S.Y(0);
  S.H(1);
  S.CX(1, 2);
  S.Sdg(2);
  S.CZ(0, 2);
  S.Z(0);
  S.S(1);
  S.CX(1, 0);
  std::vector<unsigned> v(8, 0);
  for (unsigned m = 0; m < 8; m++) {
    Simplex T(S);
    unsigned c = 0;
    c += T.MeasZ(0, (m >> 0) & 1) << 0;
    c += T.MeasZ(1, (m >> 1) & 1) << 1;
    c += T.MeasZ(2, (m >> 2) & 1) << 2;
    v[c]++;
  }
  CHECK(v[0] == 0);
  CHECK(v[1] == 4);
  CHECK(v[2] == 0);
  CHECK(v[3] == 0);
  CHECK(v[4] == 0);
  CHECK(v[5] == 0);
  CHECK(v[6] == 4);
  CHECK(v[7] == 0);
  return 0;
}
static int test_dist_5_3() {
  Simplex S(3);
  S.Y(0);
  S.H(1);
  S.CX(1, 2);
  S.Sdg(2);
  S.CZ(0, 2);
  S.Z(0);
  S.S(1);
  S.CX(1, 0);
  std::vector<unsigned> v(8, 0);
  for (unsigned m = 0; m < 8; m++) {
    Simplex T(S);
    unsigned c = 0;
    c += T.MeasX(0, (m >> 0) & 1) << 0;
    c += T.MeasY(1, (m >> 1) & 1) << 1;
    c += T.MeasZ(2, (m >> 2) & 1) << 2;
    v[c]++;
  }
  CHECK(v[0] == 1);
  CHECK(v[1] == 1);
  CHECK(v[2] == 1);
  CHECK(v[3] == 1);
  CHECK(v[4] == 1);
  CHECK(v[5] == 1);
  CHECK(v[6] == 1);
  CHECK(v[7] == 1);
  return 0;
}
static int test_dist_5_4() {
  Simplex S(3);
  S.Y(0);
  S.H(1);
  S.CX(1, 2);
  S.Sdg(2);
  S.CZ(0, 2);
  S.Z(0);
  S.S(1);
  S.CX(1, 0);
  std::vector<unsigned> v(8, 0);
  for (unsigned m = 0; m < 8; m++) {
    Simplex T(S);
    unsigned c = 0;
    c += T.MeasY(0, (m >> 0) & 1) << 0;
    c += T.MeasZ(1, (m >> 1) & 1) << 1;
    c += T.MeasX(2, (m >> 2) & 1) << 2;
    v[c]++;
  }
  CHECK(v[0] == 1);
  CHECK(v[1] == 1);
  CHECK(v[2] == 1);
  CHECK(v[3] == 1);
  CHECK(v[4] == 1);
  CHECK(v[5] == 1);
  CHECK(v[6] == 1);
  CHECK(v[7] == 1);
  return 0;
}
static int test_dist_5_5() {
  Simplex S(3);
  S.Y(0);
  S.H(1);
  S.CX(1, 2);
  S.Sdg(2);
  S.CZ(0, 2);
  S.Z(0);
  S.S(1);
  S.CX(1, 0);
  std::vector<unsigned> v(8, 0);
  for (unsigned m = 0; m < 8; m++) {
    Simplex T(S);
    unsigned c = 0;
    c += T.MeasZ(0, (m >> 0) & 1) << 0;
    c += T.MeasX(1, (m >> 1) & 1) << 1;
    c += T.MeasY(2, (m >> 2) & 1) << 2;
    v[c]++;
  }
  CHECK(v[0] == 1);
  CHECK(v[1] == 1);
  CHECK(v[2] == 1);
  CHECK(v[3] == 1);
  CHECK(v[4] == 1);
  CHECK(v[5] == 1);
  CHECK(v[6] == 1);
  CHECK(v[7] == 1);
  return 0;
}
static int test_dist_6_0() {
  Simplex S(3);
  S.H(0);
  S.H(1);
  S.H(2);
  S.Sdg(1);
  S.CZ(0, 1);
  S.S(1);
  S.CZ(1, 2);
  S.H(0);
  S.H(1);
  S.H(2);
  std::vector<unsigned> v(8, 0);
  for (unsigned m = 0; m < 8; m++) {
    Simplex T(S);
    unsigned c = 0;
    c += T.MeasX(0, (m >> 0) & 1) << 0;
    c += T.MeasX(1, (m >> 1) & 1) << 1;
    c += T.MeasX(2, (m >> 2) & 1) << 2;
    v[c]++;
  }
  CHECK(v[0] == 1);
  CHECK(v[1] == 1);
  CHECK(v[2] == 1);
  CHECK(v[3] == 1);
  CHECK(v[4] == 1);
  CHECK(v[5] == 1);
  CHECK(v[6] == 1);
  CHECK(v[7] == 1);
  return 0;
}
static int test_dist_6_1() {
  Simplex S(3);
  S.H(0);
  S.H(1);
  S.H(2);
  S.Sdg(1);
  S.CZ(0, 1);
  S.S(1);
  S.CZ(1, 2);
  S.H(0);
  S.H(1);
  S.H(2);
  std::vector<unsigned> v(8, 0);
  for (unsigned m = 0; m < 8; m++) {
    Simplex T(S);
    unsigned c = 0;
    c += T.MeasY(0, (m >> 0) & 1) << 0;
    c += T.MeasY(1, (m >> 1) & 1) << 1;
    c += T.MeasY(2, (m >> 2) & 1) << 2;
    v[c]++;
  }
  CHECK(v[0] == 1);
  CHECK(v[1] == 1);
  CHECK(v[2] == 1);
  CHECK(v[3] == 1);
  CHECK(v[4] == 1);
  CHECK(v[5] == 1);
  CHECK(v[6] == 1);
  CHECK(v[7] == 1);
  return 0;
}
static int test_dist_6_2() {
  Simplex S(3);
  S.H(0);
  S.H(1);
  S.H(2);
  S.Sdg(1);
  S.CZ(0, 1);
  S.S(1);
  S.CZ(1, 2);
  S.H(0);
  S.H(1);
  S.H(2);
  std::vector<unsigned> v(8, 0);
  for (unsigned m = 0; m < 8; m++) {
    Simplex T(S);
    unsigned c = 0;
    c += T.MeasZ(0, (m >> 0) & 1) << 0;
    c += T.MeasZ(1, (m >> 1) & 1) << 1;
    c += T.MeasZ(2, (m >> 2) & 1) << 2;
    v[c]++;
  }
  CHECK(v[0] == 2);
  CHECK(v[1] == 0);
  CHECK(v[2] == 2);
  CHECK(v[3] == 0);
  CHECK(v[4] == 0);
  CHECK(v[5] == 2);
  CHECK(v[6] == 0);
  CHECK(v[7] == 2);
  return 0;
}
static int test_dist_6_3() {
  Simplex S(3);
  S.H(0);
  S.H(1);
  S.H(2);
  S.Sdg(1);
  S.CZ(0, 1);
  S.S(1);
  S.CZ(1, 2);
  S.H(0);
  S.H(1);
  S.H(2);
  std::vector<unsigned> v(8, 0);
  for (unsigned m = 0; m < 8; m++) {
    Simplex T(S);
    unsigned c = 0;
    c += T.MeasX(0, (m >> 0) & 1) << 0;
    c += T.MeasY(1, (m >> 1) & 1) << 1;
    c += T.MeasZ(2, (m >> 2) & 1) << 2;
    v[c]++;
  }
  CHECK(v[0] == 1);
  CHECK(v[1] == 1);
  CHECK(v[2] == 1);
  CHECK(v[3] == 1);
  CHECK(v[4] == 1);
  CHECK(v[5] == 1);
  CHECK(v[6] == 1);
  CHECK(v[7] == 1);
  return 0;
}
static int test_dist_6_4() {
  Simplex S(3);
  S.H(0);
  S.H(1);
  S.H(2);
  S.Sdg(1);
  S.CZ(0, 1);
  S.S(1);
  S.CZ(1, 2);
  S.H(0);
  S.H(1);
  S.H(2);
  std::vector<unsigned> v(8, 0);
  for (unsigned m = 0; m < 8; m++) {
    Simplex T(S);
    unsigned c = 0;
    c += T.MeasY(0, (m >> 0) & 1) << 0;
    c += T.MeasZ(1, (m >> 1) & 1) << 1;
    c += T.MeasX(2, (m >> 2) & 1) << 2;
    v[c]++;
  }
  CHECK(v[0] == 1);
  CHECK(v[1] == 1);
  CHECK(v[2] == 1);
  CHECK(v[3] == 1);
  CHECK(v[4] == 1);
  CHECK(v[5] == 1);
  CHECK(v[6] == 1);
  CHECK(v[7] == 1);
  return 0;
}
static int test_dist_6_5() {
  Simplex S(3);
  S.H(0);
  S.H(1);
  S.H(2);
  S.Sdg(1);
  S.CZ(0, 1);
  S.S(1);
  S.CZ(1, 2);
  S.H(0);
  S.H(1);
  S.H(2);
  std::vector<unsigned> v(8, 0);
  for (unsigned m = 0; m < 8; m++) {
    Simplex T(S);
    unsigned c = 0;
    c += T.MeasZ(0, (m >> 0) & 1) << 0;
    c += T.MeasX(1, (m >> 1) & 1) << 1;
    c += T.MeasY(2, (m >> 2) & 1) << 2;
    v[c]++;
  }
  CHECK(v[0] == 2);
  CHECK(v[1] == 0);
  CHECK(v[2] == 0);
  CHECK(v[3] == 2);
  CHECK(v[4] == 2);
  CHECK(v[5] == 0);
  CHECK(v[6] == 0);
  CHECK(v[7] == 2);
  return 0;
}
static int test_dist_7_0() {
  Simplex S(3);
  S.X(0);
  S.H(0);
  S.CX(0, 1);
  S.S(1);
  S.Y(2);
  S.CZ(1, 2);
  S.H(2);
  S.Sdg(0);
  S.Z(2);
  S.Y(2);
  S.CX(2, 0);
  std::vector<unsigned> v(8, 0);
  for (unsigned m = 0; m < 8; m++) {
    Simplex T(S);
    unsigned c = 0;
    c += T.MeasX(0, (m >> 0) & 1) << 0;
    c += T.MeasX(1, (m >> 1) & 1) << 1;
    c += T.MeasX(2, (m >> 2) & 1) << 2;
    v[c]++;
  }
  CHECK(v[0] == 0);
  CHECK(v[1] == 0);
  CHECK(v[2] == 0);
  CHECK(v[3] == 4);
  CHECK(v[4] == 4);
  CHECK(v[5] == 0);
  CHECK(v[6] == 0);
  CHECK(v[7] == 0);
  return 0;
}
static int test_dist_7_1() {
  Simplex S(3);
  S.X(0);
  S.H(0);
  S.CX(0, 1);
  S.S(1);
  S.Y(2);
  S.CZ(1, 2);
  S.H(2);
  S.Sdg(0);
  S.Z(2);
  S.Y(2);
  S.CX(2, 0);
  std::vector<unsigned> v(8, 0);
  for (unsigned m = 0; m < 8; m++) {
    Simplex T(S);
    unsigned c = 0;
    c += T.MeasY(0, (m >> 0) & 1) << 0;
    c += T.MeasY(1, (m >> 1) & 1) << 1;
    c += T.MeasY(2, (m >> 2) & 1) << 2;
    v[c]++;
  }
  CHECK(v[0] == 1);
  CHECK(v[1] == 1);
  CHECK(v[2] == 1);
  CHECK(v[3] == 1);
  CHECK(v[4] == 1);
  CHECK(v[5] == 1);
  CHECK(v[6] == 1);
  CHECK(v[7] == 1);
  return 0;
}
static int test_dist_7_2() {
  Simplex S(3);
  S.X(0);
  S.H(0);
  S.CX(0, 1);
  S.S(1);
  S.Y(2);
  S.CZ(1, 2);
  S.H(2);
  S.Sdg(0);
  S.Z(2);
  S.Y(2);
  S.CX(2, 0);
  std::vector<unsigned> v(8, 0);
  for (unsigned m = 0; m < 8; m++) {
    Simplex T(S);
    unsigned c = 0;
    c += T.MeasZ(0, (m >> 0) & 1) << 0;
    c += T.MeasZ(1, (m >> 1) & 1) << 1;
    c += T.MeasZ(2, (m >> 2) & 1) << 2;
    v[c]++;
  }
  CHECK(v[0] == 2);
  CHECK(v[1] == 0);
  CHECK(v[2] == 0);
  CHECK(v[3] == 2);
  CHECK(v[4] == 0);
  CHECK(v[5] == 2);
  CHECK(v[6] == 2);
  CHECK(v[7] == 0);
  return 0;
}
static int test_dist_7_3() {
  Simplex S(3);
  S.X(0);
  S.H(0);
  S.CX(0, 1);
  S.S(1);
  S.Y(2);
  S.CZ(1, 2);
  S.H(2);
  S.Sdg(0);
  S.Z(2);
  S.Y(2);
  S.CX(2, 0);
  std::vector<unsigned> v(8, 0);
  for (unsigned m = 0; m < 8; m++) {
    Simplex T(S);
    unsigned c = 0;
    c += T.MeasX(0, (m >> 0) & 1) << 0;
    c += T.MeasY(1, (m >> 1) & 1) << 1;
    c += T.MeasZ(2, (m >> 2) & 1) << 2;
    v[c]++;
  }
  CHECK(v[0] == 1);
  CHECK(v[1] == 1);
  CHECK(v[2] == 1);
  CHECK(v[3] == 1);
  CHECK(v[4] == 1);
  CHECK(v[5] == 1);
  CHECK(v[6] == 1);
  CHECK(v[7] == 1);
  return 0;
}
static int test_dist_7_4() {
  Simplex S(3);
  S.X(0);
  S.H(0);
  S.CX(0, 1);
  S.S(1);
  S.Y(2);
  S.CZ(1, 2);
  S.H(2);
  S.Sdg(0);
  S.Z(2);
  S.Y(2);
  S.CX(2, 0);
  std::vector<unsigned> v(8, 0);
  for (unsigned m = 0; m < 8; m++) {
    Simplex T(S);
    unsigned c = 0;
    c += T.MeasY(0, (m >> 0) & 1) << 0;
    c += T.MeasZ(1, (m >> 1) & 1) << 1;
    c += T.MeasX(2, (m >> 2) & 1) << 2;
    v[c]++;
  }
  CHECK(v[0] == 1);
  CHECK(v[1] == 1);
  CHECK(v[2] == 1);
  CHECK(v[3] == 1);
  CHECK(v[4] == 1);
  CHECK(v[5] == 1);
  CHECK(v[6] == 1);
  CHECK(v[7] == 1);
  return 0;
}
static int test_dist_7_5() {
  Simplex S(3);
  S.X(0);
  S.H(0);
  S.CX(0, 1);
  S.S(1);
  S.Y(2);
  S.CZ(1, 2);
  S.H(2);
  S.Sdg(0);
  S.Z(2);
  S.Y(2);
  S.CX(2, 0);
  std::vector<unsigned> v(8, 0);
  for (unsigned m = 0; m < 8; m++) {
    Simplex T(S);
    unsigned c = 0;
    c += T.MeasZ(0, (m >> 0) & 1) << 0;
    c += T.MeasX(1, (m >> 1) & 1) << 1;
    c += T.MeasY(2, (m >> 2) & 1) << 2;
    v[c]++;
  }
  CHECK(v[0] == 1);
  CHECK(v[1] == 1);
  CHECK(v[2] == 1);
  CHECK(v[3] == 1);
  CHECK(v[4] == 1);
  CHECK(v[5] == 1);
  CHECK(v[6] == 1);
  CHECK(v[7] == 1);
  return 0;
}
static int test_dist_8_0() {
  Simplex S(4);
  S.H(0);
  S.H(1);
  S.H(2);
  S.H(3);
  S.CZ(0, 1);
  S.CZ(2, 3);
  S.H(0);
  S.H(1);
  S.H(2);
  S.H(3);
  S.CZ(0, 2);
  S.CZ(1, 3);
  S.H(0);
  S.H(1);
  S.H(2);
  S.H(3);
  S.CZ(0, 3);
  S.CZ(1, 2);
  S.H(0);
  S.H(1);
  S.H(2);
  S.H(3);
  S.CZ(0, 1);
  S.CZ(1, 2);
  S.CZ(1, 3);
  S.H(0);
  S.H(1);
  S.H(2);
  S.H(3);
  std::vector<unsigned> v(16, 0);
  for (unsigned m = 0; m < 16; m++) {
    Simplex T(S);
    unsigned c = 0;
    c += T.MeasX(0, (m >> 0) & 1) << 0;
    c += T.MeasX(1, (m >> 1) & 1) << 1;
    c += T.MeasX(2, (m >> 2) & 1) << 2;
    c += T.MeasX(3, (m >> 3) & 1) << 3;
    v[c]++;
  }
  CHECK(v[0] == 1);
  CHECK(v[1] == 1);
  CHECK(v[2] == 1);
  CHECK(v[3] == 1);
  CHECK(v[4] == 1);
  CHECK(v[5] == 1);
  CHECK(v[6] == 1);
  CHECK(v[7] == 1);
  CHECK(v[8] == 1);
  CHECK(v[9] == 1);
  CHECK(v[10] == 1);
  CHECK(v[11] == 1);
  CHECK(v[12] == 1);
  CHECK(v[13] == 1);
  CHECK(v[14] == 1);
  CHECK(v[15] == 1);
  return 0;
}
static int test_dist_8_1() {
  Simplex S(4);
  S.H(0);
  S.H(1);
  S.H(2);
  S.H(3);
  S.CZ(0, 1);
  S.CZ(2, 3);
  S.H(0);
  S.H(1);
  S.H(2);
  S.H(3);
  S.CZ(0, 2);
  S.CZ(1, 3);
  S.H(0);
  S.H(1);
  S.H(2);
  S.H(3);
  S.CZ(0, 3);
  S.CZ(1, 2);
  S.H(0);
  S.H(1);
  S.H(2);
  S.H(3);
  S.CZ(0, 1);
  S.CZ(1, 2);
  S.CZ(1, 3);
  S.H(0);
  S.H(1);
  S.H(2);
  S.H(3);
  std::vector<unsigned> v(16, 0);
  for (unsigned m = 0; m < 16; m++) {
    Simplex T(S);
    unsigned c = 0;
    c += T.MeasY(0, (m >> 0) & 1) << 0;
    c += T.MeasY(1, (m >> 1) & 1) << 1;
    c += T.MeasY(2, (m >> 2) & 1) << 2;
    c += T.MeasY(3, (m >> 3) & 1) << 3;
    v[c]++;
  }
  CHECK(v[0] == 0);
  CHECK(v[1] == 2);
  CHECK(v[2] == 2);
  CHECK(v[3] == 0);
  CHECK(v[4] == 2);
  CHECK(v[5] == 0);
  CHECK(v[6] == 0);
  CHECK(v[7] == 2);
  CHECK(v[8] == 2);
  CHECK(v[9] == 0);
  CHECK(v[10] == 0);
  CHECK(v[11] == 2);
  CHECK(v[12] == 0);
  CHECK(v[13] == 2);
  CHECK(v[14] == 2);
  CHECK(v[15] == 0);
  return 0;
}
static int test_dist_8_2() {
  Simplex S(4);
  S.H(0);
  S.H(1);
  S.H(2);
  S.H(3);
  S.CZ(0, 1);
  S.CZ(2, 3);
  S.H(0);
  S.H(1);
  S.H(2);
  S.H(3);
  S.CZ(0, 2);
  S.CZ(1, 3);
  S.H(0);
  S.H(1);
  S.H(2);
  S.H(3);
  S.CZ(0, 3);
  S.CZ(1, 2);
  S.H(0);
  S.H(1);
  S.H(2);
  S.H(3);
  S.CZ(0, 1);
  S.CZ(1, 2);
  S.CZ(1, 3);
  S.H(0);
  S.H(1);
  S.H(2);
  S.H(3);
  std::vector<unsigned> v(16, 0);
  for (unsigned m = 0; m < 16; m++) {
    Simplex T(S);
    unsigned c = 0;
    c += T.MeasZ(0, (m >> 0) & 1) << 0;
    c += T.MeasZ(1, (m >> 1) & 1) << 1;
    c += T.MeasZ(2, (m >> 2) & 1) << 2;
    c += T.MeasZ(3, (m >> 3) & 1) << 3;
    v[c]++;
  }
  CHECK(v[0] == 4);
  CHECK(v[1] == 0);
  CHECK(v[2] == 0);
  CHECK(v[3] == 0);
  CHECK(v[4] == 4);
  CHECK(v[5] == 0);
  CHECK(v[6] == 0);
  CHECK(v[7] == 0);
  CHECK(v[8] == 0);
  CHECK(v[9] == 0);
  CHECK(v[10] == 0);
  CHECK(v[11] == 4);
  CHECK(v[12] == 0);
  CHECK(v[13] == 0);
  CHECK(v[14] == 0);
  CHECK(v[15] == 4);
  return 0;
}
static int test_dist_8_3() {
  Simplex S(4);
  S.H(0);
  S.H(1);
  S.H(2);
  S.H(3);
  S.CZ(0, 1);
  S.CZ(2, 3);
  S.H(0);
  S.H(1);
  S.H(2);
  S.H(3);
  S.CZ(0, 2);
  S.CZ(1, 3);
  S.H(0);
  S.H(1);
  S.H(2);
  S.H(3);
  S.CZ(0, 3);
  S.CZ(1, 2);
  S.H(0);
  S.H(1);
  S.H(2);
  S.H(3);
  S.CZ(0, 1);
  S.CZ(1, 2);
  S.CZ(1, 3);
  S.H(0);
  S.H(1);
  S.H(2);
  S.H(3);
  std::vector<unsigned> v(16, 0);
  for (unsigned m = 0; m < 16; m++) {
    Simplex T(S);
    unsigned c = 0;
    c += T.MeasX(0, (m >> 0) & 1) << 0;
    c += T.MeasY(1, (m >> 1) & 1) << 1;
    c += T.MeasZ(2, (m >> 2) & 1) << 2;
    c += T.MeasX(3, (m >> 3) & 1) << 3;
    v[c]++;
  }
  CHECK(v[0] == 1);
  CHECK(v[1] == 1);
  CHECK(v[2] == 1);
  CHECK(v[3] == 1);
  CHECK(v[4] == 1);
  CHECK(v[5] == 1);
  CHECK(v[6] == 1);
  CHECK(v[7] == 1);
  CHECK(v[8] == 1);
  CHECK(v[9] == 1);
  CHECK(v[10] == 1);
  CHECK(v[11] == 1);
  CHECK(v[12] == 1);
  CHECK(v[13] == 1);
  CHECK(v[14] == 1);
  CHECK(v[15] == 1);
  return 0;
}
static int test_dist_8_4() {
  Simplex S(4);
  S.H(0);
  S.H(1);
  S.H(2);
  S.H(3);
  S.CZ(0, 1);
  S.CZ(2, 3);
  S.H(0);
  S.H(1);
  S.H(2);
  S.H(3);
  S.CZ(0, 2);
  S.CZ(1, 3);
  S.H(0);
  S.H(1);
  S.H(2);
  S.H(3);
  S.CZ(0, 3);
  S.CZ(1, 2);
  S.H(0);
  S.H(1);
  S.H(2);
  S.H(3);
  S.CZ(0, 1);
  S.CZ(1, 2);
  S.CZ(1, 3);
  S.H(0);
  S.H(1);
  S.H(2);
  S.H(3);
  std::vector<unsigned> v(16, 0);
  for (unsigned m = 0; m < 16; m++) {
    Simplex T(S);
    unsigned c = 0;
    c += T.MeasY(0, (m >> 0) & 1) << 0;
    c += T.MeasZ(1, (m >> 1) & 1) << 1;
    c += T.MeasX(2, (m >> 2) & 1) << 2;
    c += T.MeasY(3, (m >> 3) & 1) << 3;
    v[c]++;
  }
  CHECK(v[0] == 0);
  CHECK(v[1] == 0);
  CHECK(v[2] == 2);
  CHECK(v[3] == 2);
  CHECK(v[4] == 2);
  CHECK(v[5] == 2);
  CHECK(v[6] == 0);
  CHECK(v[7] == 0);
  CHECK(v[8] == 0);
  CHECK(v[9] == 0);
  CHECK(v[10] == 2);
  CHECK(v[11] == 2);
  CHECK(v[12] == 2);
  CHECK(v[13] == 2);
  CHECK(v[14] == 0);
  CHECK(v[15] == 0);
  return 0;
}
static int test_dist_8_5() {
  Simplex S(4);
  S.H(0);
  S.H(1);
  S.H(2);
  S.H(3);
  S.CZ(0, 1);
  S.CZ(2, 3);
  S.H(0);
  S.H(1);
  S.H(2);
  S.H(3);
  S.CZ(0, 2);
  S.CZ(1, 3);
  S.H(0);
  S.H(1);
  S.H(2);
  S.H(3);
  S.CZ(0, 3);
  S.CZ(1, 2);
  S.H(0);
  S.H(1);
  S.H(2);
  S.H(3);
  S.CZ(0, 1);
  S.CZ(1, 2);
  S.CZ(1, 3);
  S.H(0);
  S.H(1);
  S.H(2);
  S.H(3);
  std::vector<unsigned> v(16, 0);
  for (unsigned m = 0; m < 16; m++) {
    Simplex T(S);
    unsigned c = 0;
    c += T.MeasZ(0, (m >> 0) & 1) << 0;
    c += T.MeasX(1, (m >> 1) & 1) << 1;
    c += T.MeasY(2, (m >> 2) & 1) << 2;
    c += T.MeasZ(3, (m >> 3) & 1) << 3;
    v[c]++;
  }
  CHECK(v[0] == 2);
  CHECK(v[1] == 0);
  CHECK(v[2] == 2);
  CHECK(v[3] == 0);
  CHECK(v[4] == 2);
  CHECK(v[5] == 0);
  CHECK(v[6] == 2);
  CHECK(v[7] == 0);
  CHECK(v[8] == 0);
  CHECK(v[9] == 2);
  CHECK(v[10] == 0);
  CHECK(v[11] == 2);
  CHECK(v[12] == 0);
  CHECK(v[13] == 2);
  CHECK(v[14] == 0);
  CHECK(v[15] == 2);
  return 0;
}
static int test_dist_9_0() {
  Simplex S(2);
  S.H(0);
  S.CX(0, 1);
  S.CZ(0, 1);
  S.H(0);
  S.H(1);
  std::vector<unsigned> v(4, 0);
  for (unsigned m = 0; m < 4; m++) {
    Simplex T(S);
    unsigned c = 0;
    c += T.MeasX(0, (m >> 0) & 1) << 0;
    c += T.MeasX(1, (m >> 1) & 1) << 1;
    v[c]++;
  }
  CHECK(v[0] == 2);
  CHECK(v[1] == 0);
  CHECK(v[2] == 0);
  CHECK(v[3] == 2);
  return 0;
}
static int test_dist_9_1() {
  Simplex S(2);
  S.H(0);
  S.CX(0, 1);
  S.CZ(0, 1);
  S.H(0);
  S.H(1);
  std::vector<unsigned> v(4, 0);
  for (unsigned m = 0; m < 4; m++) {
    Simplex T(S);
    unsigned c = 0;
    c += T.MeasY(0, (m >> 0) & 1) << 0;
    c += T.MeasY(1, (m >> 1) & 1) << 1;
    v[c]++;
  }
  CHECK(v[0] == 2);
  CHECK(v[1] == 0);
  CHECK(v[2] == 0);
  CHECK(v[3] == 2);
  return 0;
}
static int test_dist_9_2() {
  Simplex S(2);
  S.H(0);
  S.CX(0, 1);
  S.CZ(0, 1);
  S.H(0);
  S.H(1);
  std::vector<unsigned> v(4, 0);
  for (unsigned m = 0; m < 4; m++) {
    Simplex T(S);
    unsigned c = 0;
    c += T.MeasZ(0, (m >> 0) & 1) << 0;
    c += T.MeasZ(1, (m >> 1) & 1) << 1;
    v[c]++;
  }
  CHECK(v[0] == 0);
  CHECK(v[1] == 2);
  CHECK(v[2] == 2);
  CHECK(v[3] == 0);
  return 0;
}
static int test_dist_9_3() {
  Simplex S(2);
  S.H(0);
  S.CX(0, 1);
  S.CZ(0, 1);
  S.H(0);
  S.H(1);
  std::vector<unsigned> v(4, 0);
  for (unsigned m = 0; m < 4; m++) {
    Simplex T(S);
    unsigned c = 0;
    c += T.MeasX(0, (m >> 0) & 1) << 0;
    c += T.MeasY(1, (m >> 1) & 1) << 1;
    v[c]++;
  }
  CHECK(v[0] == 1);
  CHECK(v[1] == 1);
  CHECK(v[2] == 1);
  CHECK(v[3] == 1);
  return 0;
}
static int test_dist_9_4() {
  Simplex S(2);
  S.H(0);
  S.CX(0, 1);
  S.CZ(0, 1);
  S.H(0);
  S.H(1);
  std::vector<unsigned> v(4, 0);
  for (unsigned m = 0; m < 4; m++) {
    Simplex T(S);
    unsigned c = 0;
    c += T.MeasY(0, (m >> 0) & 1) << 0;
    c += T.MeasZ(1, (m >> 1) & 1) << 1;
    v[c]++;
  }
  CHECK(v[0] == 1);
  CHECK(v[1] == 1);
  CHECK(v[2] == 1);
  CHECK(v[3] == 1);
  return 0;
}
static int test_dist_9_5() {
  Simplex S(2);
  S.H(0);
  S.CX(0, 1);
  S.CZ(0, 1);
  S.H(0);
  S.H(1);
  std::vector<unsigned> v(4, 0);
  for (unsigned m = 0; m < 4; m++) {
    Simplex T(S);
    unsigned c = 0;
    c += T.MeasZ(0, (m >> 0) & 1) << 0;
    c += T.MeasX(1, (m >> 1) & 1) << 1;
    v[c]++;
  }
  CHECK(v[0] == 1);
  CHECK(v[1] == 1);
  CHECK(v[2] == 1);
  CHECK(v[3] == 1);
  return 0;
}

// [END] Autogenerated by `gen-test-cases` script

static int test_phase() {
  { // Phase difference of pi/2 between X;Y and Z
    Simplex S(1), T(1);
    S.X(0); S.Y(0);
    T.Z(0);
    CHECK((T.phase() - S.phase() - 2) % 8 == 0);
  }
  { // Phase difference of pi between H;Y;H and Y
    Simplex S(1), T(1);
    S.H(0); S.Y(0); S.H(0);
    T.Y(0);
    CHECK((T.phase() - S.phase() - 4) % 8 == 0);
  }
  return 0;
}

static int test_circ1() {
  // https://github.com/CQCL/simplex/issues/2
  Simplex S(8);
  S.H(2);
  S.H(4);
  S.H(7);
  S.CX(2, 1);
  S.CX(4, 6);
  S.CX(7, 5);
  S.H(4);
  S.H(6);
  S.CX(4, 1);
  S.CX(2, 3);
  S.CX(5, 6);
  S.H(4);
  S.H(7);
  S.CX(3, 0);
  S.CX(6, 4);
  S.CX(0, 5);
  CHECK(S.is_deterministic());
  return 0;
}

static int test_reset() {
  Simplex S(3);
  S.H(0);
  S.CX(0, 1);
  S.S(2);
  S.H(2);
  S.CX(2, 1);
  S.ResetX(0);
  S.ResetY(1);
  S.ResetZ(2);
  int b0 = S.MeasX(0);
  int b1 = S.MeasY(1);
  int b2 = S.MeasZ(2);
  CHECK(b0 == 0);
  CHECK(b1 == 0);
  CHECK(b2 == 0);
  return 0;
}

int main() {
  CHECK_OK(test_X());
  CHECK_OK(test_Y());
  CHECK_OK(test_Z());
  CHECK_OK(test_H());
  CHECK_OK(test_S());
  CHECK_OK(test_Sdg());
  CHECK_OK(test_CX());
  CHECK_OK(test_CZ());
  CHECK_OK(test_MeasX());
  CHECK_OK(test_MeasY());
  CHECK_OK(test_MeasZ());
  CHECK_OK(test_many_qubits());
  CHECK_OK(test_complicated_circuit_1());
  CHECK_OK(test_complicated_circuit_2());
  CHECK_OK(test_complicated_x_cx_circuit());
  CHECK_OK(test_invert());
  CHECK_OK(test_measurements());
  CHECK_OK(test_mid_circ_meas());
  // [BEGIN] Autogenerated by `gen-test-cases` script
  CHECK_OK(test_dist_0_0());
  CHECK_OK(test_dist_0_1());
  CHECK_OK(test_dist_0_2());
  CHECK_OK(test_dist_0_3());
  CHECK_OK(test_dist_0_4());
  CHECK_OK(test_dist_0_5());
  CHECK_OK(test_dist_1_0());
  CHECK_OK(test_dist_1_1());
  CHECK_OK(test_dist_1_2());
  CHECK_OK(test_dist_1_3());
  CHECK_OK(test_dist_1_4());
  CHECK_OK(test_dist_1_5());
  CHECK_OK(test_dist_2_0());
  CHECK_OK(test_dist_2_1());
  CHECK_OK(test_dist_2_2());
  CHECK_OK(test_dist_2_3());
  CHECK_OK(test_dist_2_4());
  CHECK_OK(test_dist_2_5());
  CHECK_OK(test_dist_3_0());
  CHECK_OK(test_dist_3_1());
  CHECK_OK(test_dist_3_2());
  CHECK_OK(test_dist_3_3());
  CHECK_OK(test_dist_3_4());
  CHECK_OK(test_dist_3_5());
  CHECK_OK(test_dist_4_0());
  CHECK_OK(test_dist_4_1());
  CHECK_OK(test_dist_4_2());
  CHECK_OK(test_dist_4_3());
  CHECK_OK(test_dist_4_4());
  CHECK_OK(test_dist_4_5());
  CHECK_OK(test_dist_5_0());
  CHECK_OK(test_dist_5_1());
  CHECK_OK(test_dist_5_2());
  CHECK_OK(test_dist_5_3());
  CHECK_OK(test_dist_5_4());
  CHECK_OK(test_dist_5_5());
  CHECK_OK(test_dist_6_0());
  CHECK_OK(test_dist_6_1());
  CHECK_OK(test_dist_6_2());
  CHECK_OK(test_dist_6_3());
  CHECK_OK(test_dist_6_4());
  CHECK_OK(test_dist_6_5());
  CHECK_OK(test_dist_7_0());
  CHECK_OK(test_dist_7_1());
  CHECK_OK(test_dist_7_2());
  CHECK_OK(test_dist_7_3());
  CHECK_OK(test_dist_7_4());
  CHECK_OK(test_dist_7_5());
  CHECK_OK(test_dist_8_0());
  CHECK_OK(test_dist_8_1());
  CHECK_OK(test_dist_8_2());
  CHECK_OK(test_dist_8_3());
  CHECK_OK(test_dist_8_4());
  CHECK_OK(test_dist_8_5());
  CHECK_OK(test_dist_9_0());
  CHECK_OK(test_dist_9_1());
  CHECK_OK(test_dist_9_2());
  CHECK_OK(test_dist_9_3());
  CHECK_OK(test_dist_9_4());
  CHECK_OK(test_dist_9_5());
  // [END] Autogenerated by `gen-test-cases` script
  CHECK_OK(test_phase());
  CHECK_OK(test_circ1());
  CHECK_OK(test_reset());
  return 0;
}
