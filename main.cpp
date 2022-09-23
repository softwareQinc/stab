#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <random>
#include <chrono>
#include <vector>

#include <Eigen/Dense>
#include <qpp/qpp.h>

// include own library **always** last
#include "AffineState.h"
#include "random.h"
#include "qasm/qasm.hpp"

std::string random_qasm(int nq, bool measure) {
    // Generate random OPENQASM 2.0 string with or without measurements
    std::vector<std::string> gates = {"x",   "y",  "z",  "s",   "h",
                                      "sdg", "cx", "cz", "swap"};

    std::string qasm = "OPENQASM 2.0;\ninclude \"qelib1.inc\";\nqreg q[" +
                       std::to_string(nq) + "];\ncreg c[" + std::to_string(nq) +
                       "];\n";

    // Now add some gates
    for (int gatenumber = 0; gatenumber < pow(nq, 2); ++gatenumber) {
        int randno = stab::random_integer(0, 8);
        std::string gate = gates[randno];

        int q1 = stab::random_integer(0, nq - 1);
        int q2; // Only needed if randno > 5:
        if (randno > 5) {
            if (nq == 1) { // Can't have two-qubit gates in this case
                continue;
            }
            q2 = q1;
            while (q2 == q1) {
                q2 = stab::random_integer(0, nq - 1);
            }
        }

        qasm += gate + " q[" + std::to_string(q1) + "]";
        if (randno > 5) {
            qasm += ",q[" + std::to_string(q2) + "]";
        }
        qasm += ";\n";
    }

    if (measure) {
        for (int i = 0; i < nq; ++i) {
            qasm += "measure q[" + std::to_string(i) + "] -> c[" +
                    std::to_string(i) + "];\n";
        }
    }

    return qasm;
}

int main() {
    using namespace stab;
}