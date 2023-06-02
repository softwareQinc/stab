#include <string>
#include <vector>

#include <qpp/qpp.h>

#include "gtest/gtest.h"

#include "AffineState.h"
#include "random.h"
#include "qasm/qasm.hpp"

using namespace stab;

//anonymous namespace, now those functions are local to this translation unit
namespace {
std::string random_qasm(int nq, bool measure) {
    // We will need to generate a bunch of OPENQASM 2.0 strings, both with and without measurements
    std::vector<std::string> gates = {"x",   "y",  "z",  "s",   "h",
                                      "sdg", "cx", "cz", "swap"};

    // We now build a string called qasm by appending a bunch of random gates
    std::string qasm = "OPENQASM 2.0;\ninclude \"qelib1.inc\";\nqreg q[" +
                       std::to_string(nq) + "];\ncreg c[" + std::to_string(nq) +
                       "];\n";

    // Start with layer of Hadamards to make things interesting:

    for (int i=0;i<nq;++i){
        qasm += "h q[" + std::to_string(i) + "];\n";
    }

    // Now add n^3 random gates
    // for (int gatenumber = 0; gatenumber < pow(nq, 4); ++gatenumber) {
    for (int gatenumber = 0; gatenumber < 10000; ++gatenumber) { // TODO: Change back to line above
        int randno = random_integer(0, 8);
        std::string gate = gates[randno]; // Select random number

        // Now we choose random qubits to apply that gate to
        int q1 = random_integer(0, nq - 1);
        int q2; // Only needed if randno > 5, since otherwise it's a one-qubit gate:
        if (randno > 5) {
            if (nq == 1) { // Can't have two-qubit gates in this case
                continue;
            }
            q2 = q1;
            while (q2 == q1) { // Hack-y way to choose a distinct random number
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
    // Generate random OPENQASM 2.0 string for an IDENTITY circuit. We construct
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
get_random_circuits(int nmin, int nmax) {
    // Get random circuits with 1, ..., nmax qubits, both with and without measurements
    std::pair<std::vector<std::string>, std::vector<std::string>> all_circs;
    for (int n = nmin; n <= nmax; ++n) { // TODO: Change back to int n = 1; ...
        std::string s = random_qasm(n, false);
        all_circs.first.push_back(s);
        all_circs.second.push_back(s + "measure q -> c;\n");
    }
    return all_circs;
}

// We will use a bunch of random circuits. We can save time by only generating them once
const std::pair<std::vector<std::string>, std::vector<std::string>> circs =
    get_random_circuits(4, 5); // TODO: Change back to 10ish
const std::vector<std::string> circs_without = circs.first;
const std::vector<std::string> circs_with = circs.second;
} // namespace

// TEST(QppWorking, AllZeroState) {
//     // test that qpp is working
//     using namespace qpp;
//     ket psi = 0_ket;
//     EXPECT_EQ(psi, 0_ket);
// }

// TEST(RandomQASM, Generation) {
//     // test that qasm string generation works
//     for (int n = 1; n < 6; ++n) {
//         std::string s = random_qasm(n, false);
//     }
//     EXPECT_TRUE(true);
// }

// TEST(CompareWithQpp, Measurements) {
//     // Run the circuit with qpp, then check that the measured outcome is possible according to stab
//     bool success = true;
//     for (auto const& s : circs_without) {
//         std::istringstream prog_stream(s);
//         AffineState psi1 =
//             stab::qasm_simulator::simulate_and_return(prog_stream);

//         prog_stream.str(s);  // Reset prog stream
//         prog_stream.clear(); // Reset EOF bit

//         qpp::QCircuit qc = qpp::qasm::read(prog_stream);
//         for (int i = 0; i < psi1.n(); ++i) {
//             qc.measure(i, i);
//         }

//         qpp::QEngine qe(qc);
//         qe.execute();
//         std::vector<size_t> results = qe.get_dits();

//         /*try {*/
//         for (int i = 0; i < psi1.n(); ++i) {
//             psi1.MeasureZ(i, true, results[i]);
//         }
//         /*} catch (std::logic_error& e) {
//             std::cout << s << "\n\n";
//             for (auto i : results) std::cout << i << " ";
//         }*/

//         for (int i = 0; i < psi1.n(); ++i) {
//             assert(psi1.b()(i) == results[i]);
//         }
//     }
//     EXPECT_TRUE(success);
// }

// TEST(GenerateRandomStates, CheckNorms) {
//     // Test that the AffineState::to_ket() function gives unit vector
//     bool success = true;
//     for (auto const& s : circs_without) {
//         std::istringstream prog_stream(s);
//         AffineState psi =
//             stab::qasm_simulator::simulate_and_return(prog_stream);

//         Eigen::VectorXcd vec = psi.to_ket();
//         success = (abs(vec.norm() - 1) < 1e-12);
//         if (!success) {
//             break;
//         }
//     }
//     EXPECT_TRUE(success);
// }

// TEST(RunRandomQASM, PerformMeasurements) {
//     // Check that simulator runs qasm code and gives a sensible output state
//     bool success = true;
//     for (auto const& s : circs_with) {
//         std::istringstream prog_stream(s);
//         AffineState psi =
//             stab::qasm_simulator::simulate_and_return(prog_stream);
//         if (psi.r() > 0 || psi.A().size() > 0 || psi.Q().size() > 0) {
//             success = false;
//             break;
//         }
//     }
//     EXPECT_TRUE(success);
// }

// TEST(RunRandomQASM, Identity) {
//     // Make sure that random identity circuits indeed evaluate to |0^n>
//     // |0^n> is encoded by psi.b() being equal to zero, and A and Q having size
//     // zero.
//     bool success = true;
//     for (int nq = 20; nq <= 50; nq += 5) {
//         std::string s = random_identity(nq);
//         std::istringstream prog_stream(s);
//         AffineState psi =
//             stab::qasm_simulator::simulate_and_return(prog_stream);
//         if (psi.r() != 0 || psi.A().size() != 0 || psi.Q().size() != 0 ||
//             psi.phase() != 0 || !psi.b().isZero()) {
//             success = false;
//             break;
//         }
//     }
//     EXPECT_TRUE(success);
// }

// TODO check this please, fails once in a blue moon
// To run multiple times and stop when the test fails, execute (from inside
// ./build): ./unit_tests/unit_tests --gtest_filter=CompareWithQpp.\* --gtest_repeat=1000 --gtest_throw_on_failure
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
            std::ofstream out("qasm_code.txt");
            out << s;
            out.close();
            std::cout << "\n" << s << "\n";
            std::cout << vec1 - vec2 << "\n";
            std::cout << (vec1 - vec2).norm() << "\n";
            break;
        }
    }
    EXPECT_TRUE(success);
}

// TEST(CompareWithQpp, SpecficCircuit){
//     std::ifstream ifs("/Users/alexkerzner/test_circuit_jun_1.qasm");
//     std::string s( (std::istreambuf_iterator<char>(ifs) ),
//                        (std::istreambuf_iterator<char>()    ) );

//     bool success = true;
//     std::istringstream prog_stream(s);
//     AffineState psi1 =
//         stab::qasm_simulator::simulate_and_return(prog_stream);
//     Eigen::VectorXcd vec1 = psi1.to_ket();
//     prog_stream.str(s);  // Reset prog stream
//     prog_stream.clear(); // Reset EOF bit

//     qpp::QCircuit qc = qpp::qasm::read(prog_stream);
//     qpp::QEngine qe(qc);
//     qpp::ket vec2 = qe.execute().get_psi();

//     // std::cout << vec1;
//     // std::cout << vec2;

//     std::cout << (vec1 - vec2).norm();
//     EXPECT_TRUE(success);
// }
