#include "AffineState.h"
#include "random.h"
#include "qasm/qasm.hpp"

#include "gtest/gtest.h"

#include <qpp/qpp.h>

#include <iostream>


using namespace stab;

AffineState random_state() { return AffineState{3}; }

TEST(Test1, Case1) {
    std::cout << "hello first test!\n";

    // test that qpp is working
    using namespace qpp;
    ket psi = 0_ket;

    EXPECT_EQ(psi, 0_ket);
}

TEST(Test1, Case2) { EXPECT_EQ(0, 0); }

// this test will fail
TEST(Test2, AllTests) { EXPECT_NEAR(0, 1, 1e-3); }
