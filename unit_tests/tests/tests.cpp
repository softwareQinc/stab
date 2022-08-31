#include "AffineState.h"
#include "random.h"
#include "qasm/qasm.hpp"

#include "gtest/gtest.h"

#include <iostream>

using namespace stab;

TEST(Test1, Case1) {
    std::cout << "hello first test!\n";
    EXPECT_EQ(0, 0);
}

TEST(Test1, Case2) { EXPECT_EQ(0, 0); }

// this test will fail
TEST(Test2, AllTests) { EXPECT_NEAR(0, 1, 1e-3); }
