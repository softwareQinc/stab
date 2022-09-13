#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <random>
#include <chrono>

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

    AffineState psi(64);

    std::fstream infile("random_64_stim.stim");
    std::string line;
    while (std::getline(infile, line)) {
        std::istringstream iss(line);
        std::string op;
        int a, b;
        iss >> op;
        if (op == "tick") {
            continue;
        } else if (op == "S") {
            iss >> a;
            psi.S(a);
        } else if (op == "H") {
            iss >> a;
            psi.H(a);
        } else if (op == "M") {
            iss >> a;
            psi.MeasureZ(a);
        } else if (op == "CNOT") {
            iss >> a >> b;
            psi.CX(a, b);
        }
    }

    std::cout << psi.b();




    /*std::vector<std::pair<int, double>> times;

    for (int nq = 200; nq <=400; nq += 25) {
        std::string prog = random_qasm(nq, true);
        std::istringstream prog_stream(prog);

        std::string nqstring = "circuit" + std::to_string(nq) + ".txt";
        std::fstream qasm_file(nqstring.c_str(), std::fstream::out);
        qasm_file << prog;
        qasm_file.close();

        auto start = std::chrono::steady_clock::now();
        AffineState psi = stab::qasm_simulator::simulate_and_return(prog_stream);
        auto end = std::chrono::steady_clock::now();
        std::chrono::duration<double> diff = end - start;

        times.push_back(std::make_pair(nq, diff.count()));
    }

    std::fstream myfile("all_times.csv", std::fstream::out);
    for (auto p : times) {
        myfile << p.first << "," << p.second << "\n";
    }
    myfile.close();*/

    std::cout << std::endl;
}