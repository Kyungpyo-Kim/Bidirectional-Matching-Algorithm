/**
 * @file test.cc
 * @author kyungpyo-kim (kyungpyo94@gmail.com)
 * @date 2022-08-24
 * @brief test for Hungarian algorithm
 */

#include "hungarian_assigner.hpp"
#include <iostream>

// Run all the tests that were declared with TEST()
int main(int argc, char **argv)
{
  {
    HungarianAssigner assigner;
    HungarianAssigner::CostType cost_matrix = {
        {1, 2, 3},
        {2, 4, 6},
        {3, 6, 9},
    };
    std::vector<size_t> assignment_index;

    auto cost = assigner.solve(cost_matrix, 3, 3, 0, &assignment_index);
    std::cout << "cost: " << cost << std::endl;
    for (size_t i = 0; i < assignment_index.size(); ++i)
    {
      std::cout << "assignment_index[" << i << "]: " << assignment_index[i] << std::endl;
    }
  }

  {
    // HungarianAssigner assigner;
    // HungarianAssigner::CostType cost_matrix = {
    //     {3, 7,  5, 11},
    //     {5, 4,  6, 3},
    //     {6, 10, 1, 1},
    // };
    // std::vector<size_t> assignment_index;

    // auto cost = assigner.solve(cost_matrix, 3, 4, 0, &assignment_index);
    // std::cout << "cost: " << cost << std::endl;
    // for (size_t i = 0; i < assignment_index.size(); ++i)
    // {
    //   std::cout << "assignment_index[" << i << "]: " << assignment_index[i] << std::endl;
    // }
  }

  {
    // HungarianAssigner assigner;
    // HungarianAssigner::CostType cost_matrix = {
    //     {3, 7,  5, 11},
    //     {5, 4,  6, 3},
    //     {6, 10, 1, 1},
    // };
    // std::vector<size_t> assignment_index;

    // auto cost = assigner.solve(cost_matrix, 3, 4, 1, &assignment_index);
    // std::cout << "cost: " << cost << std::endl;
    // for (size_t i = 0; i < assignment_index.size(); ++i)
    // {
    //   std::cout << "assignment_index[" << i << "]: " << assignment_index[i] << std::endl;
    // }
  }

  return 0;
}
