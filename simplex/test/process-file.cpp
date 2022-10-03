#include <ostream>
#include <simplex.hpp>
#include <cstdlib>
#include <iostream>

int main(int argc, char *argv[]) {
  if (argc < 2 || argc > 3) {
    std::cerr << "Usage: " << argv[0] << "FILE [SEED]" << std::endl;
  }
  int seed = 0;
  if (argc == 3) {
    seed = atol(argv[2]);
  }
  Simplex S(argv[1], seed);
  unsigned n = S.n();
  for (unsigned i = 0; i < n; i++) {
    std::cout << S.MeasZ(i);
  }
  std::cout << std::endl;
  return 0;
}
