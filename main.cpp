#include <iostream>

#include <Eigen/Dense>

#include "AffineState.h"
#include "random.h"
#include "qasm/qasm.hpp"

int main() {
    using namespace stab;

    std::string prog( // Random Clifford circuit:
        "OPENQASM 2.0;\ninclude \"qelib1.inc\";\nqreg q[14];\nh q[0];\ns q[1];\nh q[2];\nh q[3];\ns q[4];\ns "
        "q[5];\ns q[6];\ns q[7];\nh q[8];\ns q[8];\ns q[9];\nh q[9];\ns "
        "q[10];\ns q[11];\nh q[12];\ns q[12];\ns q[13];\nh q[13];\nswap "
        "q[2],q[0];\ncx q[2],q[1];\ncx q[2],q[3];\ncx q[2],q[5];\ncx "
        "q[2],q[6];\ncx q[2],q[0];\ncx q[9],q[2];\ncx q[13],q[2];\ncx "
        "q[4],q[7];\ncx q[2],q[4];\nh q[4];\ncx q[4],q[2];\ncx q[10],q[8];\ncx "
        "q[8],q[2];\ncx q[2],q[10];\ncx q[12],q[11];\ncx q[11],q[2];\ncx "
        "q[2],q[12];\nh q[0];\ns q[3];\ns q[4];\nh q[4];\ns q[4];\ns q[7];\ns "
        "q[8];\nh q[8];\ns q[8];\ns q[9];\nh q[9];\nh q[10];\ns q[12];\nh "
        "q[13];\nswap q[7],q[0];\ncx q[7],q[10];\ncx q[5],q[7];\ncx "
        "q[9],q[7];\ncx q[13],q[7];\ncx q[3],q[12];\ncx q[3],q[0];\ncx "
        "q[7],q[3];\nh q[3];\ncx q[3],q[7];\ncx q[4],q[1];\ncx q[1],q[7];\ncx "
        "q[7],q[4];\ncx q[11],q[8];\ncx q[8],q[7];\ncx q[7],q[11];\nh q[0];\ns "
        "q[3];\nh q[3];\ns q[4];\ns q[12];\ns q[13];\nswap q[1],q[3];\ncx "
        "q[1],q[0];\ncx q[1],q[4];\ncx q[1],q[6];\ncx q[1],q[13];\ncx "
        "q[8],q[1];\ncx q[10],q[1];\ncx q[3],q[1];\ncx q[9],q[12];\ncx "
        "q[1],q[9];\nh q[9];\ncx q[9],q[1];\nh q[0];\ns q[0];\nh q[3];\ns "
        "q[6];\nh q[6];\ns q[6];\nh q[8];\ns q[10];\nh q[10];\nh q[12];\ns "
        "q[13];\nh q[13];\nswap q[12],q[0];\ncx q[12],q[5];\ncx "
        "q[12],q[0];\ncx q[4],q[12];\ncx q[8],q[12];\ncx q[13],q[12];\ncx "
        "q[6],q[3];\ncx q[3],q[12];\ncx q[12],q[6];\ncx q[11],q[10];\ncx "
        "q[10],q[12];\ncx q[12],q[11];\ns q[0];\nh q[6];\nh q[8];\ns q[9];\nh "
        "q[11];\nswap q[3],q[0];\ncx q[3],q[13];\ncx q[4],q[3];\ncx "
        "q[6],q[3];\ncx q[3],q[0];\nh q[0];\ncx q[0],q[3];\ncx q[9],q[8];\ncx "
        "q[8],q[3];\ncx q[3],q[9];\ncx q[11],q[10];\ncx q[10],q[3];\ncx "
        "q[3],q[11];\ns q[5];\nh q[6];\ns q[8];\nh q[8];\ns q[10];\nh "
        "q[10];\ns q[10];\nh q[11];\nh q[13];\nswap q[5],q[0];\ncx "
        "q[5],q[6];\ncx q[5],q[13];\ncx q[5],q[0];\ncx q[8],q[5];\ncx "
        "q[9],q[5];\ncx q[11],q[10];\ncx q[10],q[5];\ncx q[5],q[11];\nh "
        "q[8];\ns q[9];\nh q[9];\nh q[10];\ns q[10];\nh q[11];\ns q[11];\ncx "
        "q[13],q[8];\ncx q[9],q[0];\ncx q[0],q[13];\ncx q[13],q[9];\ncx "
        "q[11],q[10];\ncx q[10],q[13];\ncx q[13],q[11];\ns q[0];\ns q[4];\ns "
        "q[6];\nh q[6];\nh q[8];\ns q[9];\nh q[10];\nh q[11];\ncx "
        "q[6],q[11];\ncx q[8],q[11];\ncx q[10],q[11];\ncx q[11],q[4];\nh "
        "q[4];\ncx q[4],q[11];\ncx q[9],q[0];\ncx q[0],q[11];\ncx "
        "q[11],q[9];\ns q[0];\nh q[0];\ns q[0];\ns q[4];\nh q[4];\ns q[6];\nh "
        "q[6];\ns q[6];\nh q[9];\ns q[10];\nh q[10];\ncx q[4],q[9];\ncx "
        "q[10],q[4];\ncx q[6],q[0];\ncx q[0],q[4];\ncx q[4],q[6];\nh q[0];\nh "
        "q[6];\ns q[8];\nswap q[6],q[0];\ncx q[6],q[9];\ncx q[6],q[0];\ncx "
        "q[8],q[10];\ncx q[6],q[8];\nh q[8];\ncx q[8],q[6];\nh q[0];\ns "
        "q[9];\nh q[9];\ncx q[8],q[10];\ncx q[9],q[0];\ncx q[0],q[8];\ncx "
        "q[8],q[9];\ns q[10];\ncx q[10],q[0];\ncx q[0],q[9];\ncx "
        "q[9],q[10];\ns q[10];\nh q[10];\nswap q[10],q[0];\ncx q[0],q[10];\ns "
        "q[0];\nx q[0];\ny q[1];\ny q[2];\nx q[3];\nx q[4];\nz q[5];\nz "
        "q[7];\nz q[8];\ny q[9];\ny q[10];\nx q[11];\nz q[13];\n"
        );
    std::istringstream prog_stream(prog);
    qasm::simulate(prog_stream);

    std::cout << std::endl;
}