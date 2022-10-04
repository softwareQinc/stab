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
