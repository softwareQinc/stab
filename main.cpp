#include <iostream>

#include <Eigen/Dense>

#include "AffineState.h"
#include "random.h"
#include "qasm/qasm.hpp"

int main() {
    using namespace stab;

    std::string prog("OPENQASM 2.0;\ninclude \"qelib1.inc\";\nqreg q[2];\ncreg c[2];\nh q[0];\ncx q[0], q[1];\nmeasure q[0] -> c[0];\nmeasure q[1] -> c[1];\n");
    std::istringstream prog_stream(prog);
    qasm::simulate(prog_stream);

    std::cout << std::endl;
}
