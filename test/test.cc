/**
 * @file test.cc
 * @author kyungpyo-kim (kyungpyo94@gmail.com)
 * @date 2022-08-24
 * @brief test for Hungarian algorithm
 */

#include <gtest/gtest.h>

#include <memory>

#include "hungarian_assigner.hpp"

TEST(HungarianAssigner, declare) {
  try {
    HungarianAssigner assigner;
    EXPECT_TRUE(true);
  } catch (...) {
    FAIL();
  }
}

// Run all the tests that were declared with TEST()
int main(int argc, char **argv) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
