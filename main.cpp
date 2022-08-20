#include <iostream>

#include <Eigen/Dense>

#include "AffineState.h"
#include "random.h"

void aReduceGramRowCol(int c, Eigen::MatrixXi &Q) {
    Q(c, c) = (4 + Q(c, c) % 4) % 4;
    for (int i = 0; i < Q.cols(); i++) {
        if (i != c) {
            Q(c, i) = (2 + Q(c, i) % 2) % 2;
            Q(i, c) = (2 + Q(i, c) % 2) % 2;
        }
    }
}

int main() {
    /*MatrixXi Q_;
    Q_.setRandom(4, 4);
    std::cout << Q_ << std::endl;
    aReduceGramRowCol(3, Q_);
    std::cout << Q_ << std::endl;*/

    AffineState psi(3);
    psi.X(0);
    psi.X(2);
    psi.H(0);
    psi.H(1);
    psi.CZ(0, 2);
    psi.H(0);
    psi.Z(2);
    psi.CZ(0, 1);
    /*psi.H(1);
    psi.H(2);
    psi.CZ(0, 1);
    psi.CZ(0, 2);
    psi.H(0);
    psi.H(1);
    psi.H(2);*/
    std::cout << psi;

    // for (int i = 1; i < psi.n_; i++) {
    //	psi.CZ(0, i);
    // }

    // for (int i = 0; i < psi.n_; i++) {
    //	psi.H(i);
    // }
    // for (int i = 1; i < psi.n_; i++) {
    //	psi.CZ(0, i);
    // }
    // std::cout << psi;

    std::cout << '\n';
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
        std::cout << random_real(-10.0, 10.0) << ' ';
}
