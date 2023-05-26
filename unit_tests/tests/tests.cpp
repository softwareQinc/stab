#include <string>
#include <vector>

#include <qpp/qpp.h>

#include "gtest/gtest.h"

#include "AffineState.h"
#include "qasm/qasm.hpp"
#include "random.h"

using namespace stab;

// anonymous namespace, now those functions are local to this translation unit
namespace {
std::string random_qasm(int nq, bool measure) {
    // Generate random OPENQASM 2.0 string with or without measurements
    std::vector<std::string> gates = {"x",   "y",  "z",  "s",   "h",
                                      "sdg", "cx", "cz", "swap"};

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

std::string random_identity(int nq) {
    // Generate random OPENQASM 2.0 string for an identity circuit. We construct
    // this in the form UU^\dagger
    std::vector<std::string> gates = {"x",   "y",  "z",  "s",   "h",
                                      "sdg", "cx", "cz", "swap"};

    std::string qasm = "OPENQASM 2.0;\ninclude \"qelib1.inc\";\nqreg q[" +
                       std::to_string(nq) + "];\n";
    std::string inverse = "";

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

        std::string next_line = gate + " q[" + std::to_string(q1) + "]";
        if (randno > 5) {
            next_line += ",q[" + std::to_string(q2) + "]";
        }
        next_line += ";\n";
        qasm += next_line;

        inverse = next_line + inverse;
        if (gate == "s" ||
            gate == "sdg") { // Lazy way to handle inverses: s^{-1} = s^{3},
                             // sdg^{-1} = sdg^{3}
            inverse = next_line + inverse;
            inverse = next_line + inverse;
        }
    }

    return qasm + inverse;
}

std::pair<std::vector<std::string>, std::vector<std::string>>
get_random_circuits(int nmax) {
    // Get a bunch of circuits, both with an without measurements
    std::pair<std::vector<std::string>, std::vector<std::string>> all_circs;
    for (int n = 1; n <= nmax; ++n) {
        std::string s = random_qasm(n, false);
        all_circs.first.push_back(s);
        all_circs.second.push_back(s + "measure q -> c;\n");
    }
    return all_circs;
}

// Save time by only generating the random circuits once
const std::pair<std::vector<std::string>, std::vector<std::string>> circs =
    get_random_circuits(10);
const std::vector<std::string> circs_without = circs.first;
const std::vector<std::string> circs_with = circs.second;
} // namespace

TEST(QppWorking, AllZeroState) {
    // test that qpp is working
    using namespace qpp;
    ket psi = 0_ket;
    EXPECT_EQ(psi, 0_ket);
}

TEST(RandomQASM, Generation) {
    // test that qasm string generation works
    for (int n = 1; n < 6; ++n) {
        std::string s = random_qasm(n, false);
    }
    EXPECT_TRUE(true);
}

TEST(CompareWithQpp, Measurements) {
    // Measure with qpp, then check whether measurement outcome is possible with
    // stab
    bool success = true;
    for (auto const& s : circs_without) {
        std::istringstream prog_stream(s);
        AffineState psi1 =
            stab::qasm_simulator::simulate_and_return(prog_stream);

        prog_stream.str(s);  // Reset prog stream
        prog_stream.clear(); // Reset EOF bit

        qpp::QCircuit qc = qpp::qasm::read(prog_stream);
        for (int i = 0; i < psi1.n(); ++i) {
            qc.measure(i, i);
        }

        qpp::QEngine qe(qc);
        qe.execute();
        std::vector<size_t> results = qe.get_dits();

        /*try {*/
        for (int i = 0; i < psi1.n(); ++i) {
            psi1.MeasureZ(i, true, results[i]);
        }
        /*} catch (std::logic_error& e) {
            std::cout << s << "\n\n";
            for (auto i : results) std::cout << i << " ";
        }*/

        for (int i = 0; i < psi1.n(); ++i) {
            assert(psi1.b()(i) == results[i]);
        }
    }
    EXPECT_TRUE(success);
}

TEST(GenerateRandomStates, CheckNorms) {
    // Test that the AffineState::to_ket() function gives unit vector
    bool success = true;
    for (auto const& s : circs_without) {
        std::istringstream prog_stream(s);
        AffineState psi =
            stab::qasm_simulator::simulate_and_return(prog_stream);

        Eigen::VectorXcd vec = psi.to_ket();
        success = (abs(vec.norm() - 1) < 1e-12);
        if (!success) {
            break;
        }
    }
    EXPECT_TRUE(success);
}

TEST(RunRandomQASM, PerformMeasurements) {
    // Check that simulator runs qasm code and gives a sensible output state
    bool success = true;
    for (auto const& s : circs_with) {
        std::istringstream prog_stream(s);
        AffineState psi =
            stab::qasm_simulator::simulate_and_return(prog_stream);
        if (psi.r() > 0 || psi.A().size() > 0 || psi.Q().size() > 0) {
            success = false;
            break;
        }
    }
    EXPECT_TRUE(success);
}

// TEST(RunRandomQASM, SamplingCorrectness) {
//     // Check that measurement outcomes are indeed possible
//     for (int nq = 25; nq < 51; nq += 5) {
//         std::string s = random_qasm(nq, true);
//         std::istringstream prog_stream(s);
//         AffineState psi =
//             stab::qasm_simulator::simulate_and_return(prog_stream);
//         mat_u_t A = psi.A();
//         vec_u_t b = psi.b();
//         int r = psi.r();
//
//         std::map<vec_u_t, int> results = psi.Sample(100);
//     }
//
//     EXPECT_TRUE(true);
// }

TEST(RunRandomQASM, Identity) {
    // Make sure that random identity circuits indeed evaluate to |0^n>
    // |0^n> is encoded by psi.b() being equal to zero, and A and Q having size
    // zero.
    bool success = true;
    for (int nq = 20; nq <= 50; nq += 10) {
        std::string s = random_identity(nq);
        std::istringstream prog_stream(s);
        AffineState psi =
            stab::qasm_simulator::simulate_and_return(prog_stream);
        if (psi.r() != 0 || psi.A().size() != 0 || psi.Q().size() != 0 ||
            psi.phase() != 0 || !psi.b().isZero()) {
            success = false;
            break;
        }
    }
    EXPECT_TRUE(success);
}

// TODO check this please, fails once in a blue moon
// To run multiple times and stop when the test fails, execute (from inside
// ./build): ./unit_tests/unit_tests --gtest_filter=CompareWithQpp.*
// --gtest_repeat=1000 --gtest_throw_on_failure
TEST(CompareWithQpp, NoMeasurements) {
    // Run the same circuit with qpp and stab and compare results
    bool success = true;
    for (auto const& s : circs_without) {
        std::istringstream prog_stream(s);
        AffineState psi1 =
            stab::qasm_simulator::simulate_and_return(prog_stream);
        Eigen::VectorXcd vec1 = psi1.to_ket();
        prog_stream.str(s);  // Reset prog stream
        prog_stream.clear(); // Reset EOF bit

        qpp::QCircuit qc = qpp::qasm::read(prog_stream);
        qpp::QEngine qe(qc);
        qpp::ket vec2 = qe.execute().get_psi();
        // Recall qpp::ket is just Eigen::VectorXcd so vec1 and vec2 are
        // compatible

        /*auto diff = vec1 - vec2;
        for (int i = 0; i < pow(2, psi1.n()); ++i) {
            if (abs(diff(i)) > 1e-14) {
                success = false;
                break;
            }
        }

        if (!success) {
            std::cout << s << "\n\n";
        }*/

        if ((vec1 - vec2).norm() > 1e-12) {
            success = false;
            break;
        }
    }
    EXPECT_TRUE(success);
}
