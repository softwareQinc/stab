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
    stab::AffineState psi(nq);
    for (int i = 0; i < nq; ++i) {  //
        psi.H(i);
    }  //
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
            //psi.H(a);
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
    Eigen::MatrixXi m(4, 4);
    m << 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16;
    std::cout << m << "\n\n";

    using block_t = Eigen::Block<Eigen::MatrixXi>; // this is our type
    // this will be our smart pointer to matrix that'll be shared by all member functions;
    // takes care of memory allocations/de-allocations automatically with zero overhead
    std::unique_ptr<block_t> Q;

    Q = std::make_unique<block_t>(m, 0, 0, 2, 2);
    std::cout << Q->rows() << "x" << Q->cols() << "\n";
    std::cout << *Q << "\n\n"; // use *Q to dereference the pointer and refer to the matrix
    (*Q)(0, 0) = 100;

    Q = std::make_unique<block_t>(m, 0, 0, 3, 3);
    std::cout << Q->rows() << "x" << Q->cols() << "\n";
    std::cout << *Q << "\n\n";

    std::cout << (*Q) * (*Q) << std::endl;
}
