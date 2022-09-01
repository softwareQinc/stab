#include "AffineState.h"
#include "random.h"
#include "qasm/qasm.hpp"

#include "gtest/gtest.h"

#include <qpp/qpp.h>

#include <iostream>

#include <string>


using namespace stab;

AffineState random_state(int nq) {
    // Generate random state (not uniformly random but probably good enough)
    // Use fact that all stabilizer states are graph states plus local Cliffords
    AffineState psi(nq);
    // Generate graph state
    for (int i = 0; i < nq; ++i) { psi.H(i); }
    for (int i = 0; i < nq; ++i) {
        for (int j = 0; j < i; ++j) {
            if (random_bit(0.5) == 0) {
                psi.CZ(i, j);
            }
        }
    }

    // Apply local Clifford using decomposition S^a H^b S^c, where 0<=a,c<=3 and 0<=b<=1
    for (int i = 0; i < nq; ++i) { // Really silly way of doing this, but it's adequate for now
        if (random_bit(0.5) == 0) {
            psi.S(i);
        }
        if (random_bit(0.5) == 0) {
            psi.S(i);
        }
        if (random_bit(0.75) == 0) {
            psi.H(i);
            if (random_bit(0.5) == 0) {
                psi.S(i);
            }
            if (random_bit(0.5) == 0) {
                psi.S(i);
            }
        }
    }
    return psi;
}


std::string random_qasm(int nq, bool measure) {
    std::vector<std::string> gates = {"x",   "y",  "z",  "s",   "h",
                                      "sdg", "cx", "cz", "swap"};
    // Boilerplate stuff at top of qasm code:
    std::string qasm = "OPENQASM 2.0;\ninclude \"qelib1.inc\";\nqreg[" +
                       std::to_string(nq) + "];\n";
    if (measure) {
        qasm += "creg[" + std::to_string(nq) + "];\n";
    }

    // Now add some gates

    for (int gatenumber = 0; gatenumber < pow(nq, 3); ++gatenumber) {
        int randno = random_integer(0, 8);
        std::string gate = gates[randno];
        
        int q1 = random_integer(0, nq - 1);
        int q2; // Only needed if randno > 5:
        if (randno > 5) {
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

TEST(QppWorking, AllZeroState) {
    std::cout << "hello first test!\n";

    // test that qpp is working
    using namespace qpp;
    ket psi = 0_ket;

    EXPECT_EQ(psi, 0_ket);
}

TEST(GenerateRandomStates, JustGenerate) {
    bool success = true;
    for (int n = 1; n < 11; ++n) {
        AffineState psi = random_state(n);
    }
    EXPECT_TRUE(success);
}

TEST(GenerateRandomStates, CheckNorms) {
    bool success = true;
    for (int n = 1; n < 10; ++n) {
        AffineState psi = random_state(n);
        Eigen::VectorXcd vec = psi.to_vec();
        success = (abs(vec.norm() - 1) < 1e-4);
    }
    EXPECT_TRUE(success);
}

TEST(CompareWithQPP, RandomQasm) {
    bool success = true;
    for (int nq = 1; nq < 10; ++nq) {
        std::string qsm = random_qasm(nq, false);
        qasm::QASMSimulator qs(nq);
        qs.run();
    }
}

// this test will pass
TEST(Test2, AllTests) { EXPECT_NEAR(0.999999, 1, 1e-3); }
