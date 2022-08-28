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

TEST(HungarianAssigner, solve_basic1) {
  HungarianAssigner assigner;
  HungarianAssigner::CostType cost_matrix = {
      {3, 7, 5, 11},
      {5, 4, 6, 3},
      {6, 10, 1, 1},
  };
  HungarianAssigner::AssignmentType assignment_index;

  EXPECT_EQ(assigner.solve(cost_matrix, 3, 4, 0, &assignment_index), 7);
  EXPECT_EQ(assignment_index[0], 0);
  EXPECT_EQ(assignment_index[1], 3);
  EXPECT_EQ(assignment_index[2], 2);
}

TEST(HungarianAssigner, solve_basic2) {
  HungarianAssigner assigner;
  HungarianAssigner::CostType cost_matrix = {
      {3, 7, 5, 11},
      {5, 4, 6, 3},
      {6, 10, 1, 1},
  };
  HungarianAssigner::AssignmentType assignment_index;

  EXPECT_EQ(assigner.solve(cost_matrix, 3, 4, 1, &assignment_index), 27);
  EXPECT_EQ(assignment_index[0], 3);
  EXPECT_EQ(assignment_index[1], 2);
  EXPECT_EQ(assignment_index[2], 1);
}

TEST(HungarianAssigner, solve_complicated) {
  HungarianAssigner assigner;
  HungarianAssigner::CostType cost_matrix = {
    {300, 290, 280, 290, 210, 300, 290, 280, 290, 210},
    {250, 310, 290, 300, 200, 250, 310, 290, 300, 200},
    {180, 190, 300, 190, 180, 180, 190, 300, 190, 180},
    {320, 180, 190, 240, 170, 320, 180, 190, 240, 170},
    {270, 210, 190, 250, 160, 270, 210, 190, 250, 160},
    {190, 200, 220, 190, 140, 190, 200, 220, 190, 140},
    {220, 300, 230, 180, 160, 220, 300, 230, 180, 160},
    {260, 190, 260, 210, 180, 260, 190, 260, 210, 180}
  };
  HungarianAssigner::AssignmentType assignment_index;

  EXPECT_EQ(assigner.solve(cost_matrix, 8, 10, 0, &assignment_index), 1520);
}


// Run all the tests that were declared with TEST()
int main(int argc, char **argv) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
