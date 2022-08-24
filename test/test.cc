/**
 * @file test.cc
 * @author kyungpyo-kim (kyungpyo94@gmail.com)
 * @date 2022-08-24
 * @brief test for Hungarian algorithm
 */

#include <gtest/gtest.h>

#include <memory>

#include "hungarian_assigner.hpp"

TEST(HungarianAssigner, declare)
{
  try
  {
    HungarianAssigner assigner;
    EXPECT_TRUE(true);
  }
  catch (...)
  {
    FAIL();
  }
}

TEST(HungarianAssigner, solve_basic)
{
  HungarianAssigner assigner;
  HungarianAssigner::CostType cost_matrix = {
    {3, 7, 5, 11},
    {5, 4, 6, 3},
    {6, 10, 1, 1},
  };
  std::vector<size_t> assignment_index;
  assigner.solve(cost_matrix, 3, 4, 0, &assignment_index);
}

// Run all the tests that were declared with TEST()
int main(int argc, char **argv)
{
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
