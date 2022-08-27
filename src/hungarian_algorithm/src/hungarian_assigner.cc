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
  cost_matrix_.clear();
  cost_matrix_ = cost_matrix;
  n_ = n;
  m_ = m;
  mode_ = mode;
  is_working_ = true;
  next_step_ = 0;

  // run assignment algorithm
  while (is_working_)
  {
    switch (next_step_)
    {
    case 0:
      // build cost matrix
      step0(); // O(dim_^2)
      show("step0");
      break;

    case 1:
      // step 1
      step1(); // O(dim_^2)
      show("step1");
      break;

    case 2:
      // step 2
      step2(); // O(dim_^2)
      show("step2");
      break;

    case 3:
      // step 3
      step3(); // O(dim_^2)
      show("step3");
      break;

    case 4:
      // step 4
      step4(); // O(dim_^2)
      show("step4");
      break;

    case 5:
      // step 5
      step5(); // O(dim_^2)
      show("step5");
      break;

    case 6:
      // step 6
      step6(); // O(dim_^2)
      show("step6");
      break;

    case 7:
      // step 7
      step7(); // O(dim_^2)
      show("step7");
      break;

    default:
      break;
    }

    if (next_step_ == -1)
      break;
  }

  // calculate assignment cost
  float cost = 0;
  assignment_index->clear();
  for (size_t row = 0; row < n_; ++row)
  {
    assignment_index->emplace_back(results_[row]);
    cost += cost_matrix_[row][results_[row]];
  }
  return cost;
}

void HungarianAssigner::show(std::string name = "")
{
  std::cout << name << std::endl;
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

  std::cout << "mask_: " << std::endl;
  for (size_t row = 0; row < dim_; ++row)
  {
    for (size_t col = 0; col < dim_; ++col)
    {
      std::cout << mask_[row][col] << " ";
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;
}

void HungarianAssigner::step0()
{
  /**
   * @brief step 0: initialize cost matrix, mask matrix and path matrix
   *
   * time complexity: O(dim_^2)
   *
   */

  dim_ = std::max(n_, m_);
  // to solve unbalanced problem, set cost matrix to be square matrix
  // and set 0 to the rest of the matrix
  opt_matrix_.clear();
  opt_matrix_.resize(dim_, std::vector<float>(dim_, 0));

  mask_.clear();
  mask_.resize(dim_, std::vector<int>(dim_, 0));
  path_.resize(dim_ + 2, std::vector<int>(2, 0));
  line_row_.clear();
  line_col_.clear();

  for (size_t row = 0; row < dim_; ++row)
  {
    for (size_t col = 0; col < dim_; ++col)
    {
      if (row < n_ & col < m_)
      {
        if (mode_ == 0)
        {
          // to minimize the cost, set cost matrix
          opt_matrix_[row][col] = cost_matrix_[row][col];
        }
        else if (mode_ == 1)
        {
          // to maximize the cost, set cost matrix to be negative
          opt_matrix_[row][col] = -cost_matrix_[row][col];
        }
        else
        {
          assert(false && "Invalid mode");
        }
      }
    }
  }

  next_step_ = 1;
}

void HungarianAssigner::step1()
{
  /**
   * @brief step 1: make minimum value to zero along the row
   *
   * time complexity: O(dim_^2)
   */

  for (auto &row : opt_matrix_)
  {
    auto min_value = *std::min_element(row.begin(), row.end());
    for (auto &r : row)
    {
      r -= min_value;
    }
  }

  next_step_ = 2;
}

void HungarianAssigner::step2()
{
  /**
   * @brief step 2: find zero and uncovered row and column, star zeros
   *
   * time complexity: O(dim_^2)
   */

  for (size_t col = 0; col < dim_; ++col)
  {
    for (size_t row = 0; row < dim_; ++row)
    {
      if (opt_matrix_[row][col] == 0 && line_row_.count(row) == 0 && line_col_.count(col) == 0)
      {
        mask_[row][col] = kStar;
        line_row_.insert(row);
        line_col_.insert(col);
      }
    }
  }

  next_step_ = 3;
}

void HungarianAssigner::step3()
{
  /**
   * @brief step 3: find starred zero in current col and cover it
   * 
   * 모든 column이 star 되어 있으면 종료 -> step 7
   * 그렇지 않으면 -> step 4
   * 
   * time complexity: O(dim_^2)
   */

  for (size_t row = 0; row < dim_; ++row)
  {
    for (size_t col = 0; col < dim_; ++col)
    {
      if (mask_[row][col] == kStar)
      {
        line_col_.insert(col);
      }
    }
  }

  size_t cnt_col = 0;
  for (size_t col = 0; col < dim_; ++col)
  {
    if (line_col_.count(col) == 1)
    {
      cnt_col++;
    }
  }

  if (cnt_col < dim_)
  {
    next_step_ = 4;
  }
  else
  {
    next_step_ = 7;
  }
}

void HungarianAssigner::step4()
{
  /**
   * @brief 
   * 
   * 0인 element 가 없으면 -> step 6
   * 0인 element 있으면 Prime 으로 체크하고,
   * 같은 row 에 star 가 있으면 row 를 cover 하고 column 은 uncover 한다.
   * prime 한 row 에 star 가 없으면 step 5
   * uncovered zero 가 없어질때까지 반복한다.
   * 
   * time complexity: O(dim_)
   *  uncovered zero 의 숫자는 최대 dim_ 개이다.
   */

  while (true)
  {
    auto [row, col] = step4_find_zero_cost_element();
    if (row == -1)
    {
      next_step_ = 6;
      return;
    }

    mask_[row][col] = kPrime;
    if (step4_star_in_row(row))
    {
      auto col = step4_find_col_star_in_row(row);
      line_row_.insert(row);
      line_col_.erase(col);
    }
    else
    {
      next_step_ = 5;
      path_row_0_ = row;
      path_col_0_ = col;
      return;
    }
  }
}

std::pair<int, int> HungarianAssigner::step4_find_zero_cost_element()
{
  /**
   * @brief find uncovered zero in cost matrix
   * 
   */

  bool found = false;
  int zero_row = -1, zero_col = -1;

  for (size_t row = 0; row < dim_; ++row)
  {
    for (size_t col = 0; col < dim_; ++col)
    {
      if (opt_matrix_[row][col] == 0 && line_row_.count(row) == 0 && line_col_.count(col) == 0)
      {
        found = true;
        zero_row = static_cast<int>(row);
        zero_col = static_cast<int>(col);
      }

      if (found)
      {
        break;
      }
    }
    if (found)
    {
      break;
    }
  }

  return {zero_row, zero_col};
}

bool HungarianAssigner::step4_star_in_row(const size_t row)
{
  for (size_t col = 0; col < dim_; ++col)
  {
    if (mask_[row][col] == kStar)
    {
      return true;
    }
  }
  return false;
}

size_t HungarianAssigner::step4_find_col_star_in_row(const size_t row)
{
  for (size_t col = 0; col < dim_; ++col)
  {
    if (mask_[row][col] == kStar)
    {
      return col;
    }
  }
  assert(false && "No star in row");
}

void HungarianAssigner::step5()
{
  cnt_path_ = 1;
  path_[cnt_path_ - 1][0] = path_row_0_;
  path_[cnt_path_ - 1][1] = path_col_0_;

  bool is_working = true;
  while (is_working)
  {
    auto row = find_row_star_in_col(path_[cnt_path_ - 1][1]);
    if (row > -1)
    {
      cnt_path_ += 1;
      path_[cnt_path_ - 1][0] = row;
      path_[cnt_path_ - 1][1] = path_[cnt_path_ - 2][1];
    }
    else
    {
      is_working = false;
    }

    if (is_working)
    {
      auto col = find_col_prime_in_row(path_[cnt_path_ - 1][0]);
      cnt_path_ += 1;
      path_[cnt_path_ - 1][0] = path_[cnt_path_ - 2][0];
      path_[cnt_path_ - 1][1] = col;
    }
  }
  augment_path();
  clear_lines();
  clear_primes();
  next_step_ = 3;
}

void HungarianAssigner::step6()
{
  float min_val = std::numeric_limits<float>::max();
  for (size_t row = 0; row < dim_; ++row)
  {
    for (size_t col = 0; col < dim_; ++col)
    {
      if (opt_matrix_[row][col] < min_val && line_row_.count(row) == 0 && line_col_.count(col) == 0)
      {
        min_val = opt_matrix_[row][col];
      }
    }
  }

  for (size_t row = 0; row < dim_; ++row)
  {
    for (size_t col = 0; col < dim_; ++col)
    {
      if (line_row_.count(row) == 1)
      {
        opt_matrix_[row][col] += min_val;
      }
      if (line_col_.count(col) == 0)
      {
        opt_matrix_[row][col] -= min_val;
      }
    }
  }

  next_step_ = 4;
}

int HungarianAssigner::find_row_star_in_col(const size_t col)
{
  for (size_t row = 0; row < dim_; ++row)
  {
    if (mask_[row][col] == kStar)
    {
      return row;
    }
  }

  return -1;
}

size_t HungarianAssigner::find_col_prime_in_row(const size_t row)
{
  for (size_t col = 0; col < dim_; ++col)
  {
    if (mask_[row][col] == kPrime)
    {
      return col;
    }
  }
  assert(false && "No prime in row");
}

void HungarianAssigner::augment_path()
{
  for (size_t path = 0; path < cnt_path_; ++path)
  {
    if (mask_[path_[path][0]][path_[path][1]] == kStar)
    {
      mask_[path_[path][0]][path_[path][1]] = 0;
    }
    else
    {
      mask_[path_[path][0]][path_[path][1]] = kStar;
    }
  }
}
void HungarianAssigner::clear_lines()
{
  line_row_.clear();
  line_col_.clear();
}
void HungarianAssigner::clear_primes()
{
  for (size_t row = 0; row < dim_; ++row)
  {
    for (size_t col = 0; col < dim_; ++col)
    {
      if (mask_[row][col] == kPrime)
      {
        mask_[row][col] = 0;
      }
    }
  }
}

void HungarianAssigner::step7()
{
  results_.resize(n_);
  for (size_t row = 0; row < n_; ++row)
  {
    for (size_t col = 0; col < dim_; ++col)
    {
      if (mask_[row][col] == kStar)
      {
        results_[row] = col;
      }
    }
  }

  next_step_ = -1;
}
