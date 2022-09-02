//#include <iostream> // don't include unless really needed, it's a huge header
#include <string>
#include <vector>

#include <qpp/qpp.h>

#include "gtest/gtest.h"

// include own library **always** last
#include "AffineState.h"
#include "random.h"
#include "qasm/qasm.hpp"

using namespace stab;

// anonymous namespace, now this function is local to this translation unit
namespace {
    std::string random_qasm(int nq, bool measure) {
        std::vector<std::string> gates = {"x", "y", "z", "s", "h",
                                          "sdg", "cx", "cz", "swap"};
        // Boilerplate stuff at top of qasm code:
        std::string qasm = "OPENQASM 2.0;\ninclude \"qelib1.inc\";\nqreg q[" +
                           std::to_string(nq) + "];\ncreg c[" + std::to_string(nq) +
                           "];\n";

        // Now add some gates

        for (int gatenumber = 0; gatenumber < pow(nq, 3); ++gatenumber) {
            int randno = random_integer(0, 8);
            std::string gate = gates[randno];

            int q1 = random_integer(0, nq - 1);
            int q2; // Only needed if randno > 5:
            if (randno > 5) {
                if (nq == 1) { // Can't have two-qubit gates in this case
                    continue;
                }
                q2 = q1;
                while (q2 == q1) {
                    q2 = random_integer(0, nq - 1);
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

    std::pair<std::vector<std::string>, std::vector<std::string>>
    get_random_circuits(int nmax) {
        std::pair<std::vector<std::string>, std::vector<std::string>> all_circs;
        for (int n = 1; n <= nmax; ++n) {
            std::string s = random_qasm(n, false);
            all_circs.first.push_back(s);
            all_circs.second.push_back(s + "measure q -> c;\n");
        }
        return all_circs;
    }

    const std::pair<std::vector<std::string>, std::vector<std::string>> circs =
            get_random_circuits(12);
    const std::vector<std::string> circs_without = circs.first;
    const std::vector<std::string> circs_with = circs.second;
}


TEST(QppWorking, AllZeroState) {
    // test that qpp is working
    using namespace qpp;
    ket psi = 0_ket;
    EXPECT_EQ(psi, 0_ket);
}

TEST(RandomQASM, Generation) {
    for (int n = 1; n < 6; ++n) {
        std::string s = random_qasm(n, false);
    }
    EXPECT_EQ(0, 0);
}

TEST(GenerateRandomStates, CheckNorms) {
    bool success = true;
    for (auto const &s: circs_without) {
        std::istringstream prog_stream(s);
        AffineState psi =
                stab::qasm_simulator::simulate_and_return(prog_stream);
        Eigen::VectorXcd vec = psi.to_vec();
        success = (abs(vec.norm() - 1) < 1e-8);
        if (!success) {
            break;
        }
    }
    EXPECT_TRUE(success);
}

TEST(RunRandomQASM, PerformMeasurements) {
    bool success = true;
    for (auto const &s: circs_with) {
        std::istringstream prog_stream(s);
        AffineState psi =
                stab::qasm_simulator::simulate_and_return(prog_stream);
        if (psi.r_ > 0 || psi.A_.size() > 0 || psi.Q_.size() > 0) {
            success = false;
            break;
        }
    }

    EXPECT_TRUE(success);
}

TEST(CompareWithQPP, NoMeasurements) {
    bool success = true;
    for (auto const &s: circs_without) {
        std::istringstream prog_stream(s);
        AffineState psi1 =
                stab::qasm_simulator::simulate_and_return(prog_stream);
        Eigen::VectorXcd vec1 = psi1.to_vec();
        prog_stream.str(s); // Reset prog stream
        prog_stream.clear();  // Reset EOF bit

        qpp::QCircuit qc = qpp::qasm::read(prog_stream);
        qpp::QEngine qe(qc);
        qpp::ket vec2 = qe.execute().get_psi();
        // Recall qpp::ket is just Eigen::VectorXcd

        if ((vec1 - vec2).norm() > 1e-12) {
            success = false;
            break;
        }
    }

    EXPECT_TRUE(success);
}
