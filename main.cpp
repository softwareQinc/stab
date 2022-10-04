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

int main() {
    Eigen::MatrixXd data(4, 4);
    data << 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16;
    Eigen::Ref<Eigen::MatrixXd> ref{data};
    std::cout << ref << '\n';
}