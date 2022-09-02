#include <iostream>

#include <Eigen/Dense>
#include <qpp/qpp.h>

// include own library **always** last
#include "AffineState.h"
#include "random.h"
#include "qasm/qasm.hpp"


int main() {

    using namespace stab;

    std::string prog =
        "OPENQASM 2.0;\n include \"qelib1.inc\";\n qreg q[1];\n h q[0];\n";
    std::istringstream prog_stream(prog);
    AffineState psi = stab::qasm_simulator::simulate_and_return(prog_stream);
    /*qpp::QCircuit qc = qpp::qasm::read(prog_stream);
    qpp::QEngine qe(qc);
    qpp::ket psi = qe.execute().get_psi();*/
    std::cout << psi;

    std::cout << std::endl;
}