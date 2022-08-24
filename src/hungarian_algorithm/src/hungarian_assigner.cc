/**
 * @file hungarian_algorithm.cc
 * @author kyungpyo-kim (kyungpyo94@gmail.com)
 * @date 2022-08-24
 * @brief Hungarian algorithm
 */

#include "hungarian_assigner.hpp"
#include <iostream>

float HungarianAssigner::solve(
    const CostType &cost_matrix, const size_t n,
    const size_t m, const int mode, std::vector<size_t> *assignment_index) {

  // initialize variables
  float cost = 0;
  n_ = n;
  m_ = m;
  dim_ = std::max(n, m);
  mode_ = mode;

  // build cost matrix
  build_cost_matrix(cost_matrix);

  // step 1 and step 2 to make zero elements
  

  return cost;
}

void HungarianAssigner::build_cost_matrix(const CostType &cost_matrix){
  cost_matrix_.resize(dim_, std::vector<float>(dim_, 0));
  
  for (size_t row = 0; row < dim_; ++row) {
    for (size_t col = 0; col < dim_; ++col) {
      if (row < n_ & col < m_){
        cost_matrix_[row][col] = cost_matrix[row][col];
      }
    }
  }
}


void HungarianAssigner::step1() {}
void HungarianAssigner::step2() {}
void HungarianAssigner::step3() {}
void HungarianAssigner::step4() {}
void HungarianAssigner::step5() {}
void HungarianAssigner::step6() {}
