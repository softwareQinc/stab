#include <iostream>

#include <Eigen/Dense>

#include "AffineState.h"
#include "random.h"
#include "qasm/qasm.hpp"

int main() {
    using namespace stab;

    AffineState psi(2);
    psi.A_.setZero(2,2);
    psi.A_ << 0, 1, 1, 0;
    psi.Q_.setZero(2,2);
    psi.Q_ << 3, 0, 0, 3;
    psi.b_ << 0, 1;
    psi.phase_ = 4;
    psi.pivots_[0] = 1;
    psi.pivots_[1] = 0;
    psi.r_ = 2;

    std::cout << "Initial state\n\n" << psi;

    psi.H(1);

    std::cout << psi;

    //std::string prog( // Random Clifford circuit:
    //    "OPENQASM 2.0;\ninclude \"qelib1.inc\";\nqreg q[2];\ny q[0];\ns q[0];\ny q[1];\nh q[1];\ncx q[0],q[1];\n h q[0];\n s q[0];\nsdg q[1];\n"//h q[1];\n"//cx q[0],q[1];\n"
    //    );
    //std::istringstream prog_stream(prog);
    //qasm::simulate(prog_stream);

    std::cout << std::endl;
}