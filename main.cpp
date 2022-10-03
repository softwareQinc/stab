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
    Eigen::MatrixXi master(3, 3);
    master << 1, 2, 3, 4, 5, 6, 7, 8, 9;

    Eigen::Ref<Eigen::MatrixXi> submat = master.topLeftCorner(2, 2);
    std::cout << submat << "\n\n";

    submat.row(0) += submat.row(1);
    std::cout << submat << "\n\n";
    submat.row(0).swap(submat.row(1));
    std::cout << submat;


}

//int main() {
//    using namespace stab;
//
//    std::cout << "Beginning tests \n";
//
//    std::vector<std::pair<int, double>> times;
//    int nmin = 200;
//    int nmax = 1000;
//    int step = 25;
//    int copies_per_n = 3;
//
//
//    for (int n = nmin; n <= nmax; n += step) {
//        std::cout << n << "\n";
//        //for (int j = 1; j <= copies_per_n; ++j) {
//            /*std::string fname = R"(C:\Users\Alex\Desktop\random_shallow_nonuniform_stims\)";
//            fname += "random_clifford_" + std::to_string(n) + "_" +
//                     std::to_string(j) + ".stim";
//            std::fstream infile(fname);*/
//
//            AffineState psi(n);
//            auto start = std::chrono::steady_clock::now();
//            //run_stim(infile, n);
//            for (int q = 0; q < n; ++q) {
//                psi.H(q);
//            }
//            auto end = std::chrono::steady_clock::now();
//            std::chrono::duration<double> diff = end - start;
//
//            times.push_back(std::make_pair(n, diff.count()));
//        //}
//    }
//
//    std::fstream myfile("just_hadamards_skip_initialization.csv", std::fstream::out);
//    for (auto p : times) {
//        myfile << p.first << "," << p.second << "\n";
//    }
//    myfile.close();
//}