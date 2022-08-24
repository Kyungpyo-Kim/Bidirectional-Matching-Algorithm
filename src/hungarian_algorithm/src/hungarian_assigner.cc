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
  float cost = 0;
  n_ = n;
  m_ = m;
  dim_ = std::max(n, m);
  mode_ = mode;

  is_working_ = true;
  next_step_ = 0;

  // run optimization
  while (is_working_)
  {
    switch (next_step_)
    {
    case 0:
      // build cost matrix
      step0(); // O(dim_^2)
      break;

    case 1:
      // step 1
      step1(); // O(dim_^2)
      break;

    case 2:
      // step 2
      step2(); // O(dim_^2)
      break;

    case 3:
      // step 3
      step3(); // O(dim_^2)
      break;

    case 4:
      // step 4
      step4(); // O(dim_^2)
      break;

    case 5:
      // step 5
      step5(); // O(dim_^2)
      break;

    case 6:
      // step 6
      step6(); // O(dim_^2)
      break;

    case 7:
      // step 7
      step7(); // O(dim_^2)
      break;

    default:
      break;
    }

    if (next_step_ == -1)
      break;
  }

  assignment_index->clear();
  for (size_t row = 0; row < n_; ++row)
  {
    assignment_index->emplace_back(results_[row]);
    cost += cost_matrix_[row][results_[row]];
  }
  return cost;
}

void HungarianAssigner::step0()
{
  opt_matrix_.clear();
  opt_matrix_.resize(dim_, std::vector<float>(dim_, 0));
  mask_.clear();
  mask_.resize(dim_, std::vector<int>(dim_, 0));
  path_.resize(dim_ + 2, std::vector<int>(2, 0));

  for (size_t row = 0; row < dim_; ++row)
  {
    for (size_t col = 0; col < dim_; ++col)
    {
      if (row < n_ & col < m_)
      {
        if (mode_ == 0)
        {
          opt_matrix_[row][col] = cost_matrix_[row][col];
        }
        else if (mode_ == 1)
        {
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

  show("step1");
  next_step_ = 2;
}

void HungarianAssigner::step2()
{
  line_row_.clear();
  line_col_.clear();

  /** step 2: make minimum value to zero along the column */
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
  show("step2");
}

void HungarianAssigner::step3()
{
  line_row_.clear();
  line_col_.clear();

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

  show("step3");
}

void HungarianAssigner::step4()
{
  while (true)
  {
    auto [row, col] = find_zero_cost_element();
    if (row == -1)
    {
      next_step_ = 6;
      show("step4 - 6");
      return;
    }

    mask_[row][col] = kPrime;
    if (star_in_row(row))
    {
      auto col = find_col_star_in_row(row);
      line_row_.insert(row);
      line_col_.erase(col);
    }
    else
    {
      next_step_ = 5;
      path_row_0_ = row;
      path_col_0_ = col;
      show("step4 - 5");
      return;
    }
  }
}

std::pair<int, int> HungarianAssigner::find_zero_cost_element()
{

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

bool HungarianAssigner::star_in_row(const size_t row)
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

size_t HungarianAssigner::find_col_star_in_row(const size_t row)
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
  std::cout << "1\n";
  cnt_path_ = 1;
  path_[cnt_path_ - 1][0] = path_row_0_;
  path_[cnt_path_ - 1][1] = path_col_0_;
  std::cout << "2\n";

  bool is_working = true;
  while (is_working)
  {
    auto row = find_row_star_in_col(path_[cnt_path_ - 1][1]);
    if (row > -1)
    {
      std::cout << "3\n";
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
      std::cout << "4 1\n";
      auto col = find_col_prime_in_row(path_[cnt_path_ - 1][0]);
      std::cout << "4 2\n";
      cnt_path_ += 1;
      std::cout << "4 3\n";
      std::cout << "size: " << path_.size() << std::endl;
      std::cout << "cnt_path_: " << cnt_path_ << std::endl;
      path_[cnt_path_ - 1][0] = path_[cnt_path_ - 2][0];
      std::cout << "4 4\n";
      path_[cnt_path_ - 1][1] = col;
    }
  }
  std::cout << "5\n";
  augment_path();
  std::cout << "6\n";
  clear_lines();
  std::cout << "7\n";
  clear_primes();
  std::cout << "8\n";
  next_step_ = 3;
  show("step5");
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
  show("step6");
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
  show("step7");
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