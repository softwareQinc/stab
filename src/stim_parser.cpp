#include <fstream>
#include <iostream>
#include <numeric>
#include <sstream>

#include "AffineState.h"

stab::AffineState run_stim(std::fstream& infile, int nq) {
    using namespace stab;
    AffineState psi(nq);
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
        } else { // Some error just for now
            assert(false);
        }
    }
    return psi;
}