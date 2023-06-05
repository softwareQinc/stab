#include <chrono>
#include <fstream>
#include <iostream>
#include <memory>
#include <random>
#include <sstream>
#include <string>
#include <vector>

#include <Eigen/Dense>

#ifdef USE_QPP
#include <qpp/qpp.h>
#endif // USE_QPP

// include own library **always** last
#include "AffineState.h"
#include "qasm/qasm.hpp"
#include "random.h"

stab::AffineState run_stim(std::fstream& infile, const int& nq) {
    // Naive little function to run stim files
    stab::AffineState psi(nq);
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
        } else if (op == "S_DAG") {
            iss >> a;
            psi.SDG(a);
        } else if (op == "X") {
            iss >> a;
            psi.X(a);
        } else if (op == "Y") {
            iss >> a;
            psi.Y(a);
        } else if (op == "Z") {
            iss >> a;
            psi.Z(a);
        } else if (op == "I") {
            continue;
        } else if (op == "M") {
            iss >> a;
            psi.MeasureZ(a);
        } else if (op == "CNOT") {
            iss >> a >> b;
            psi.CX(a, b);
        } else if (op == "CZ") {
            iss >> a >> b;
            psi.CZ(a, b);
        } else if (op == "SWAP") {
            iss >> a >> b;
            psi.SWAP(a, b);
        } else { // Some error just for now
            assert(false);
        }
    }
    return psi;
}

int main() {
    using namespace stab;

    // Simple example
    AffineState psi(3);
    psi.H(0);
    psi.CX(0, 1);
    psi.CX(0, 2);
    std::cout << "Statevector representation:\n" << psi.to_ket();
    int result = psi.MeasureZ(0);
    std::cout << "\nMeasured qubit 0 and observed result " << result;
    psi.MeasureZ(1, true, result);
    psi.MeasureZ(2, true, result);
    std::cout << "\nAfter postselecting, new state is:\n" << psi.to_ket() << "\n";
}