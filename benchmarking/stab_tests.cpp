#include <chrono>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "AffineState.h"

stab::AffineState run_stim(std::fstream& infile, const int& nq) {
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

int main() {
    using namespace stab;

    std::vector<std::pair<int, double>> times;

    int nmin = 200;
    int nmax = 1000;
    int step = 25;
    int copies_per_n = 3;

    for (int n = nmin; n <= nmax; n += step) {
        for (int j = 1; j <= copies_per_n; ++j) {
            std::string fname = R"(random_stims\)";
            fname += "random_clifford_" + std::to_string(n) + "_" +
                     std::to_string(j) + ".stim";
            std::fstream infile(fname);

            auto start = std::chrono::steady_clock::now();
            run_stim(infile, n);
            auto end = std::chrono::steady_clock::now();
            std::chrono::duration<double> diff = end - start;

            times.emplace_back(n, diff.count());
        }
    }

    std::fstream myfile("all_times.csv", std::fstream::out);
    for (auto p : times) {
        myfile << p.first << "," << p.second << "\n";
    }
    myfile.close();
}