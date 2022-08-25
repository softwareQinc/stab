#include <iostream>

#include <Eigen/Dense>

#include "AffineState.h"
#include "random.h"
#include "qasm/qasm.hpp"

int main() {
    using namespace stab;

    Eigen::MatrixXi A;
    A.setRandom(4, 4);
    std::cout << "A = " << A << "\n";
    Eigen::VectorXi r1 = A.row(1);
    r1(0) = 0;
    std::cout << "A = " << A << "\n";
    std::cout << "r1 = " << r1 << "\n";

    //std::string prog(
    //    "OPENQASM 2.0;\ninclude \"qelib1.inc\";\nqreg q[4];\ncreg c[4];\nh q[0];\nh q[1];\nh q[2];\nh q[3];\nz "
    //    "q[0];\nz q[2];\ns q[1];\ncx q[0], q[1];\ns q[1];\ncx q[1], q[2];\ncx "
    //    "q[2], q[3];\ncx q[0], q[1];\ns q[2];\nx q[3];\nx q[0];\ncx q[1], "
    //    "q[2];\nx q[0];\ncx q[2], q[3];\ncx q[0], q[1];\ncx q[1], q[2];\nreset "
    //    "q[1];\nmeasure q[2] -> c[2];\nmeasure q[1] -> c[1];\nmeasure q[0] -> "
    //    "c[0];\nmeasure q[3] -> c[3];\n");
    //std::istringstream prog_stream(prog);
    //qasm::simulate(prog_stream);

    std::cout << std::endl;
}