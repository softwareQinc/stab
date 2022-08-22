#include <iostream>

#include <Eigen/Dense>

#include "AffineState.h"
#include "random.h"

int main() {
    using namespace stab;

    AffineState psi(3);
    psi.H(0);
    psi.H(1);
    psi.H(2);
    psi.CZ(0, 1);
    psi.CZ(0, 2);
    psi.CZ(1, 2);
    psi.H(0);
    psi.H(1);
    psi.H(2);
    std::cout << psi;

    /*std::cout << '\n';
    for (int i = 0; i < 10; ++i)
        std::cout << random_bit(0.3) << ' ';

    std::cout << '\n';
    for (int i = 0; i < 10; ++i)
        std::cout << random_bit(0.7) << ' ';

    std::cout << '\n';
    for (int i = 0; i < 10; ++i)
        std::cout << random_bit() << ' ';

    std::cout << '\n';
    for (int i = 0; i < 10; ++i)
        std::cout << random_integer(-10, 10) << ' ';

    std::cout << '\n';
    for (int i = 0; i < 10; ++i)
        std::cout << random_real(-10.0, 10.0) << ' ';*/

    std::cout << std::endl;
}
