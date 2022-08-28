/**
 * @file test.cc
 * @author kyungpyo-kim (kyungpyo94@gmail.com)
 * @date 2022-08-24
 * @brief demo for Hungarian algorithm
 */

#include "hungarian_assigner.hpp"
#include <iostream>

// Run all the tests that were declared with TEST()
int main(int argc, char **argv) {
  {
    HungarianAssigner assigner;
    HungarianAssigner::CostType cost_matrix = {
        {1, 1, 3},
        {2, 4, 6},
        {3, 6, 9},
    };
    HungarianAssigner::AssignmentType assignment_index;

    auto cost = assigner.solve(cost_matrix, 3, 3, 0, &assignment_index);
    assigner.show("demo");
    std::cout << "cost: " << cost << std::endl;
    for (size_t i = 0; i < assignment_index.size(); ++i) {
      std::cout << "assignment_index[" << i << "]: " << assignment_index[i]
                << std::endl;
    }
    std::cout << "----------------------------------------" << std::endl;
  }

  {
    HungarianAssigner assigner;
    HungarianAssigner::CostType cost_matrix = {
        {3, 7, 5, 11},
        {5, 4, 6, 3},
        {6, 10, 1, 1},
    };
    HungarianAssigner::AssignmentType assignment_index;

    auto cost = assigner.solve(cost_matrix, 3, 4, 0, &assignment_index);
    assigner.show("demo");
    std::cout << "cost: " << cost << std::endl;
    for (size_t i = 0; i < assignment_index.size(); ++i) {
      std::cout << "assignment_index[" << i << "]: " << assignment_index[i]
                << std::endl;
    }
    std::cout << "----------------------------------------" << std::endl;
  }

  {
    HungarianAssigner assigner;
    HungarianAssigner::CostType cost_matrix = {
        {3, 7, 5, 11},
        {5, 4, 6, 3},
        {6, 10, 1, 1},
    };
    HungarianAssigner::AssignmentType assignment_index;

    auto cost = assigner.solve(cost_matrix, 3, 4, 1, &assignment_index);
    assigner.show("demo");
    std::cout << "cost: " << cost << std::endl;
    for (size_t i = 0; i < assignment_index.size(); ++i) {
      std::cout << "assignment_index[" << i << "]: " << assignment_index[i]
                << std::endl;
    }
    std::cout << "----------------------------------------" << std::endl;
  }

  {
    HungarianAssigner assigner;
    HungarianAssigner::CostType cost_matrix = {
        {3, 5, 6},
        {7, 4, 10},
        {5, 6, 1},
        {11, 3, 1},
    };

    HungarianAssigner::AssignmentType assignment_index;

    auto cost = assigner.solve(cost_matrix, 4, 3, 0, &assignment_index);
    assigner.show("demo");
    std::cout << "cost: " << cost << std::endl;
    for (size_t i = 0; i < assignment_index.size(); ++i) {
      std::cout << "assignment_index[" << i << "]: " << assignment_index[i]
                << std::endl;
    }
    std::cout << "----------------------------------------" << std::endl;
  }

  {
    HungarianAssigner assigner;
    HungarianAssigner::CostType cost_matrix = {
        {300, 290, 280, 290, 210, 300, 290, 280, 290, 210},
        {250, 310, 290, 300, 200, 250, 310, 290, 300, 200},
        {180, 190, 300, 190, 180, 180, 190, 300, 190, 180},
        {320, 180, 190, 240, 170, 320, 180, 190, 240, 170},
        {270, 210, 190, 250, 160, 270, 210, 190, 250, 160},
        {190, 200, 220, 190, 140, 190, 200, 220, 190, 140},
        {220, 300, 230, 180, 160, 220, 300, 230, 180, 160},
        {260, 190, 260, 210, 180, 260, 190, 260, 210, 180},
    };
    HungarianAssigner::AssignmentType assignment_index;

    auto cost = assigner.solve(cost_matrix, 8, 10, 0, &assignment_index);
    assigner.show("demo");

    std::cout << "cost: " << cost << std::endl;
    for (size_t i = 0; i < assignment_index.size(); ++i) {
      std::cout << "assignment_index[" << i << "]: " << assignment_index[i]
                << std::endl;
    }
    std::cout << "----------------------------------------" << std::endl;
  }

  return 0;
}
