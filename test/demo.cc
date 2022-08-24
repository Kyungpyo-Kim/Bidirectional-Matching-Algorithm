/**
 * @file test.cc
 * @author kyungpyo-kim (kyungpyo94@gmail.com)
 * @date 2022-08-24
 * @brief test for Hungarian algorithm
 */

#include "hungarian_assigner.hpp"

// Run all the tests that were declared with TEST()
int main(int argc, char **argv)
{
  HungarianAssigner assigner;
  HungarianAssigner::CostType cost_matrix = {
      {1, 2, 3},
      {2, 4, 6},
      {3, 6, 9},
  };
  std::vector<size_t> assignment_index;

  assigner.solve(cost_matrix, 3, 3, 0, &assignment_index);

  return 0;
}
