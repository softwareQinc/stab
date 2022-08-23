#include <iostream>

#include <Eigen/Dense>

#include "AffineState.h"
#include "random.h"

int main() {
    using namespace stab;

    AffineState psi(3);
    psi.X(0);
    psi.X(1);
    psi.X(2);
    psi.H(0);
    psi.H(1);
    psi.H(2);
    psi.CZ(0, 1);
    psi.CZ(0, 2);
    psi.CZ(1, 2);
    psi.H(0);
    psi.H(1);
    std::cout << psi;



    // Initialize state directly so that we don't get flooded with console output
    //AffineState psi(3);
    //psi.A_.conservativeResize(3, 2);
    //psi.A_ << 1, 1, 0, 1, 1, 0;
    //psi.Q_.conservativeResize(2, 2);
    //psi.Q_ << 2, 1, 1, 2;
    //psi.r_ = 2;
    //psi.b_ << 1, 0, 0;
    //psi.pivots_[0] = 2;
    //psi.pivots_[1] = 1;

    //psi.H(1);
    ////psi.H(2);
    //std::cout << psi;



    std::cout << std::endl;
}
