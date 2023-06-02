#include <iostream>
#include <memory>
#include <sstream>
#include <fstream>
#include <string>
#include <random>
#include <chrono>
#include <vector>
#include <memory>

#include <Eigen/Dense>
#include <qpp/qpp.h>

// include own library **always** last
#include "AffineState.h"
#include "random.h"
#include "qasm/qasm.hpp"

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

    // std::cout << "Hello";

    // std::ifstream ifs("/Users/alexkerzner/Documents/GitHub/stab/build/unit_tests/qasm_error_example.qasm");
    // std::string content( (std::istreambuf_iterator<char>(ifs) ),
    //                    (std::istreambuf_iterator<char>()    ) );
    // std::cout << content;

    // auto psi = stab::qasm_simulator::simulate_and_return(ifs);
    // std::cout << psi;
    // // // stab::qasm_simulator::simulate_file(fname);
}