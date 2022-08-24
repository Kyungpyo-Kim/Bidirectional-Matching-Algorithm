/**
 * @file hungarian_algorithm.cc
 * @author kyungpyo-kim (kyungpyo94@gmail.com)
 * @date 2022-08-24
 * @brief Hungarian algorithm
 */

#include "hungarian_assigner.hpp"
#include <iostream>
#include <algorithm>
#include <cassert>

float HungarianAssigner::solve(
    const CostType &cost_matrix, const size_t n,
    const size_t m, const int mode, std::vector<size_t> *assignment_index)
{

  // initialize variables
  float cost = 0;
  n_ = n;
  m_ = m;
  dim_ = std::max(n, m);
  mode_ = mode;

  // build cost matrix
  build_cost_matrix(cost_matrix); // O(dim_^2)

  // step 1 and step 2 to make zero elements
  step1(); // O(dim_^2)
  show();

  step2(); // O(dim_^2)
  show();

  // step 3 to cover all zeros with a minimum number of lines
  step3(); // O(dim_^3)
  show();

  // step 4 if all lines are covered, step 6 if not step 5
  step4();

  show();

  return cost;
}

void HungarianAssigner::build_cost_matrix(const CostType &cost_matrix)
{
  cost_matrix_ = cost_matrix;
  opt_matrix_.resize(dim_, std::vector<float>(dim_, 0));
  prime_.resize(dim_, std::vector<bool>(dim_, false));
  star_.resize(dim_, std::vector<bool>(dim_, false));

  for (size_t row = 0; row < dim_; ++row)
  {
    for (size_t col = 0; col < dim_; ++col)
    {
      if (row < n_ & col < m_)
      {
        if (mode_ == 0)
        {
          opt_matrix_[row][col] = cost_matrix[row][col];
        }
        else if (mode_ == 1)
        {
          opt_matrix_[row][col] = -cost_matrix[row][col];
        }
        else
        {
          assert(false && "Invalid mode");
        }
      }
    }
  }
}

void HungarianAssigner::show()
{
  std::cout << "opt_matrix_: " << std::endl;
  for (size_t row = 0; row < dim_; ++row)
  {
    for (size_t col = 0; col < dim_; ++col)
    {
      std::cout << opt_matrix_[row][col] << " ";
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;
}

void HungarianAssigner::step1()
{
  /** step 1: make minimum value to zero along the row */
  for (auto &row : opt_matrix_)
  {
    auto min_value = *std::min_element(row.begin(), row.end());
    for (auto &r : row)
    {
      r -= min_value;
      
    }
  }

  for (size_t row = 0; row < dim_; ++row)
  {
    for (size_t col = 0; col < dim_; ++col)
    {
      if (std::abs(opt_matrix_[row][col]) < epsilon_)
      {
        prime_[row][col] = true;
      }
    }
  }
}

void HungarianAssigner::step2()
{
  /** step 2: make minimum value to zero along the column */
  for (size_t col = 0; col < dim_; ++col)
  {
    std::vector<float> column(dim_, 0);
    for (size_t row = 0; row < dim_; ++row)
    {
      column[row] = opt_matrix_[row][col];
    }
    auto min_value = *std::min_element(column.begin(),
                                       column.end());
    for (size_t row = 0; row < dim_; ++row)
    {
      opt_matrix_[row][col] -= min_value;
    }
  }
}
void HungarianAssigner::step3()
{
  line_row_.clear();
  line_col_.clear();

  std::vector<std::vector<bool>> zero_checker(dim_, std::vector<bool>(dim_, false));

  for (size_t row = 0; row < dim_; ++row)
  {
    for (size_t col = 0; col < dim_; ++col)
    {
      if (std::abs(opt_matrix_[row][col]) <= epsilon_ && !zero_checker[row][col])
      {
        zero_checker[row][col] = true;

        size_t zero_in_row = 0;
        for (size_t col_check = 0; col_check < dim_; ++col_check)
        {
          if (std::abs(opt_matrix_[row][col_check]) <= epsilon_ && !zero_checker[row][col_check])
          {
            ++zero_in_row;
          }
        }

        size_t zero_in_col = 0;
        for (size_t row_check = 0; row_check < dim_; ++row_check)
        {
          if (std::abs(opt_matrix_[row_check][col]) <= epsilon_ && !zero_checker[row_check][col])
          {
            ++zero_in_col;
          }
        }

        if (zero_in_row < zero_in_col)
        {
          line_col_.insert(col);
          for (size_t row_check = 0; row_check < dim_; ++row_check)
          {
            if (std::abs(opt_matrix_[row_check][col]) <= epsilon_ && !zero_checker[row_check][col])
            {
              zero_checker[row_check][col] = true;
            }
          }
        }
        else
        {
          line_row_.insert(row);
          for (size_t col_check = 0; col_check < dim_; ++col_check)
          {
            if (std::abs(opt_matrix_[row][col_check]) <= epsilon_ && !zero_checker[row][col_check])
            {
              zero_checker[row][col_check] = true;
            }
          }
        }
      }
    }
  }
}
void HungarianAssigner::step4()
{

  if (line_row_.size() + line_col_.size() < dim_)
  {
    step5();
  }
  else
  {
    step6();
  }
}
void HungarianAssigner::step5()
{
  float min_value = std::numeric_limits<float>::max();

  for (size_t row = 0; row < dim_; ++row)
  {
    for (size_t col = 0; col < dim_; ++col)
    {
      if (line_row_.count(row) == 0 && line_col_.count(col) == 0)
      {
        if (opt_matrix_[row][col] < min_value)
        {
          min_value = opt_matrix_[row][col];
        }
      }
    }
  }

  for (size_t row = 0; row < dim_; ++row)
  {
    for (size_t col = 0; col < dim_; ++col)
    {
      if (line_row_.count(row) == 0 && line_col_.count(col) == 0)
      {
        opt_matrix_[row][col] -= min_value;
      }
      if (line_row_.count(row) == 1 && line_col_.count(col) == 1)
      {
        opt_matrix_[row][col] += min_value;
      }
    }
  }

  step3();
  step4();
}

void HungarianAssigner::step6()
{
  std::cout << "step6" << std::endl;

  std::vector<std::vector<size_t>> results(dim_);

  for (size_t row = 0; row < dim_; ++row)
  {
    for (size_t col = 0; col < dim_; ++col)
    {
      if (line_row_.count(row) == 1 || line_col_.count(col) == 1)
      {
        if (std::abs(opt_matrix_[row][col]) < epsilon_)
        {
          results[row].push_back(col);
        }
      }
    }
  }

  for (auto row : results)
  {
    std::cout << "row: ";
    for (auto col : row)
    {
      std::cout << col << " ";
    }
    std::cout << std::endl;
  }
}
