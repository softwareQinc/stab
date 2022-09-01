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
        //std::cout << "\npsi = " << psi << "\n";
        //std::cout << "\nvec is " << vec << "\n";
        //double nrm = vec.norm();
        //std::cout << "norm is " << nrm << "\n";
        success = (abs(vec.norm() - 1) < 1e-4);
    }
    EXPECT_TRUE(success);
}

// this test will pass
TEST(Test2, AllTests) { EXPECT_NEAR(0.999999, 1, 1e-3); }
