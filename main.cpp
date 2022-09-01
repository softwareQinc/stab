#include <iostream>

#include <Eigen/Dense>
#include <qpp/qpp.h>

// include own library **always** last
#include "AffineState.h"
#include "random.h"
#include "qasm/qasm.hpp"


int main() {

    using namespace stab;

    int ncols = 3;
    std::cout << pow(2, ncols/2.0);

    //std::string prog( // Random Clifford circuit:
    //        "OPENQASM 2.0;\ninclude \"qelib1.inc\";\nqreg q[10];\ns q[0];\ns q[1];\nh q[4];\nh q[5];\ns q[6];\ncx "
    //        "q[0],q[1];\ncx q[0],q[4];\ncx q[0],q[6];\ncx q[0],q[9];\ncx "
    //        "q[2],q[0];\ncx q[3],q[0];\ncx q[5],q[0];\ncx q[7],q[0];\ns q[1];\ns "
    //        "q[4];\nh q[4];\nh q[5];\nh q[6];\ns q[7];\nh q[8];\ns q[9];\nswap "
    //        "q[3],q[1];\ncx q[3],q[1];\ncx q[5],q[6];\ncx q[5],q[7];\ncx "
    //        "q[5],q[9];\ncx q[3],q[5];\nh q[5];\ncx q[5],q[3];\ncx q[8],q[4];\ncx "
    //        "q[4],q[3];\ncx q[3],q[8];\ns q[1];\ns q[2];\nh q[5];\ns q[6];\nh "
    //        "q[8];\ns q[8];\nh q[9];\ncx q[5],q[2];\ncx q[5],q[6];\ncx "
    //        "q[5],q[9];\ncx q[4],q[5];\ncx q[8],q[1];\ncx q[1],q[5];\ncx "
    //        "q[5],q[8];\ns q[1];\nh q[1];\ns q[1];\ns q[2];\nh q[4];\nh q[6];\ns "
    //        "q[6];\nh q[7];\nh q[8];\nh q[9];\ns q[9];\nswap q[2],q[1];\ncx "
    //        "q[2],q[4];\ncx q[2],q[7];\ncx q[2],q[1];\ncx q[8],q[2];\ncx "
    //        "q[9],q[6];\ncx q[6],q[2];\ncx q[2],q[9];\ns q[1];\nh q[6];\nswap "
    //        "q[6],q[1];\ncx q[6],q[1];\ncx q[9],q[6];\ncx q[6],q[8];\nh q[8];\ncx "
    //        "q[8],q[6];\ns q[1];\nh q[1];\ns q[1];\ns q[4];\ns q[7];\nh q[7];\nh "
    //        "q[8];\nh q[9];\nswap q[4],q[1];\ncx q[7],q[4];\ncx q[8],q[4];\ncx "
    //        "q[9],q[1];\ncx q[4],q[9];\nh q[9];\ncx q[9],q[4];\ns q[1];\ns "
    //        "q[7];\nh q[8];\ns q[9];\nswap q[7],q[1];\ncx q[7],q[1];\nh q[1];\ncx "
    //        "q[1],q[7];\ncx q[9],q[8];\ncx q[8],q[7];\ncx q[7],q[9];\nh q[1];\nh "
    //        "q[9];\nswap q[8],q[1];\ncx q[9],q[1];\ncx q[8],q[9];\nh q[9];\ncx "
    //        "q[9],q[8];\nh q[1];\nswap q[9],q[1];\ncx q[1],q[9];\ns q[1];\nh "
    //        "q[1];\ny q[0];\nz q[2];\nx q[3];\nz q[5];\nx q[6];\nz q[7];\ny "
    //        "q[8];\ny q[9];\n"
    //);
    //std::istringstream prog_stream(prog);
    //qasm::simulate(prog_stream);

    std::cout << std::endl;
}