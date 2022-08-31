#include "AffineState.h"
#include "random.h"
#include "qasm/qasm.hpp"

#include "gtest/gtest.h"

#include <qpp/qpp.h>

#include <iostream>


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
        if (random_bit(0.5) == 0) {
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

TEST(ProduceRandomStates, NQubits1To10) {
    std::cout << "hello first test!\n";

    // test that qpp is working
    using namespace qpp;
    ket psi = 0_ket;

    EXPECT_EQ(psi, 0_ket);
}

TEST(ConvertToStatevector, RandomStates) {
    bool success = true;
    for (int n = 1; n < 11; ++n) {
        AffineState psi(n);
        Eigen::VectorXcd vec = psi.to_vec();
    }
    EXPECT_EQ(0, 0); 
}

// this test will fail
TEST(Test2, AllTests) { EXPECT_NEAR(0, 1, 1e-3); }
