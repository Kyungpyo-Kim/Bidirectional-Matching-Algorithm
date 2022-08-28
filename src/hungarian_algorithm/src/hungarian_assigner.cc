/**
 * @file hungarian_algorithm.cc
 * @author kyungpyo-kim (kyungpyo94@gmail.com)
 * @date 2022-08-24
 * @brief Hungarian algorithm
 */

#include "hungarian_assigner.hpp"
#include <algorithm>
#include <cassert>
#include <iomanip>
#include <iostream>

float HungarianAssigner::solve(const CostType &cost_matrix, const size_t n,
                               const size_t m, const int mode,
                               std::vector<int> *assignment_index) {
  // initialize variables
  cost_matrix_.clear();
  cost_matrix_ = cost_matrix;
  n_ = n;
  m_ = m;
  mode_ = mode;
  init();

  // run assignment algorithm
  next_step_ = 0;

  while (next_step_ > -1) {
    switch (next_step_) {
    case 0:
      preliminaries();
      // show("preliminaries");
      break;

    case 1:
      step1();
      // show("step1");
      break;

    case 2:
      step2();
      // show("step2");
      break;

    case 3:
      step3();
      // show("step3");
      break;

    default:
      break;
    }
  }

  return wrapUp(assignment_index);
}

void HungarianAssigner::show(std::string name = "") {
  std::cout << "[ " << name << " ]" << std::endl;
  std::cout << "cost_matrix_: " << std::endl;
  for (size_t row = 0; row < n_; ++row) {
    for (size_t col = 0; col < m_; ++col) {
      std::cout << std::setw(8) << cost_matrix_[row][col] << " ";
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;

  std::cout << "opt_matrix_: " << std::endl;
  for (size_t row = 0; row < dim_; ++row) {
    for (size_t col = 0; col < dim_; ++col) {
      std::cout << std::setw(8) << opt_matrix_[row][col] << " ";
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;

  std::cout << "mask: " << std::endl;
  for (size_t row = 0; row < dim_; ++row) {
    for (size_t col = 0; col < dim_; ++col) {
      std::cout << mask_[row][col] << " ";
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;

  std::cout << "cover: " << std::endl;
  for (size_t row = 0; row < dim_; ++row) {
    for (size_t col = 0; col < dim_; ++col) {
      if (row_cover_.count(row) > 0 || col_cover_.count(col) > 0)
        std::cout << "1 ";
      else
        std::cout << "0 ";
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;
}

void HungarianAssigner::init() {
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

  // initialize cover set
  row_cover_.clear();
  col_cover_.clear();

  // initialize mask matrix
  mask_.clear();
  mask_.resize(dim_, std::vector<int>(dim_, 0));

  // time complexity: O(dim_^2)
  for (size_t row = 0; row < dim_; ++row) {
    for (size_t col = 0; col < dim_; ++col) {
      if (row < n_ & col < m_) {
        if (mode_ == 0) {
          // to minimize the cost, set cost matrix
          opt_matrix_[row][col] = cost_matrix_[row][col];
        } else if (mode_ == 1) {
          // to maximize the cost, set cost matrix to be negative
          opt_matrix_[row][col] = -cost_matrix_[row][col];
        } else {
          assert(false && "Invalid mode");
        }
      }
    }
  }

  next_step_ = 1;
}

void HungarianAssigner::preliminaries() {
  /**
   * @brief step 1: subtract minimum from each row and mark starred zeros
   *
   * No lines are covered; no zeros are starred or primed.
   * Subtract from each element in this row the smallest element of this row. Do
   * the same for each row of A. Then consider each column of the resulting
   * matrix and subtract from each column its smallest entry.
   * Consider a zero Z of the matrix. If there is no starred zero in its row and
   * none in its column, star Z. Repeat, considering each zero in the matrix in
   * turn. Then cover every column containing a starred zero. [These starred
   * zeros are independent.]
   *
   * time complexity: O(dim_^2)
   */

  // subtract minimum from each row
  // time complexity: O(dim_^2)
  for (auto &row : opt_matrix_) {
    auto min_value = *std::min_element(row.begin(), row.end());
    for (auto &elem : row) {
      elem -= min_value;
    }
  }

  // time complexity: O(dim_^2)
  std::unordered_set<size_t> row_check;
  std::unordered_set<size_t> col_check;
  for (size_t row = 0; row < dim_; ++row) {
    for (size_t col = 0; col < dim_; ++col) {
      if (std::abs(opt_matrix_[row][col]) <= EPSILON_FLOAT &&
          row_check.count(row) == 0 && col_check.count(col) == 0) {
        // mark starred zeros
        mask_[row][col] = kStar;
        row_check.insert(row);
        col_check.insert(col);

        // cover each column containing a starred zero
        col_cover_.insert(col);
      }
    }
  }

  if (row_cover_.size() == dim_ || col_cover_.size() == dim_) {
    next_step_ = -1;
    return;
  }

  next_step_ = 1;
  return;
}

void HungarianAssigner::step1() {
  /**
   * @brief step 1: find uncovered zeros and prime them
   *
   * Choose a non-covered zero and prime it. Consider the row containing it. If
   * there is no starred zero in this row, go at once to Step 2. If there is a
   * starred zero Z in this row, cover this row and uncover the column of Z.
   * Repeat until all zeros are covered. Go to Step 3.
   *
   * time complexity: O(dim_^3)
   */

  while (true) { // O(dim_) uncovered zeros must be less than or equal to dim_
    auto [row, col] = uncovered_zero(); // O(dim_^2)
    if (row == -1 || col == -1) {
      next_step_ = 3;
      return;
    }

    mask_[row][col] = kPrime;
    auto starred_zero = starred_zero_in_row(row); // O(dim_)

    if (starred_zero == -1) {
      uncovered_primed_zero_ = {row, col};
      next_step_ = 2;
      return;
    }

    row_cover_.insert(row);
    col_cover_.erase(starred_zero);
  }
}

void HungarianAssigner::step2() {
  /**
   * @brief step 2: make sequence of alternating primed and starred zeros
   *
   * There is a sequence of alternating starred and primed zeros, constructed as
   * follows: Let Zo denote the uncovered 0'. [There is only one.] Let Z1 denote
   * the 0* in Z0's column (if any). Let Z2 denote the 0' in Z1's row (we must
   * prove that it exists). Let Z3 denote the 0* in Z2's column (if any).
   * Similarly continue until the sequence stops at a 0', Z2k, which has no 0*
   * in its column (this we must also prove).
   * Unstar each starred zero of the sequence and star each primed zero of the
   * sequence.
   * Erase all primes, uncover every row, and cover every column
   * containing a 0*. If all columns are covered, the starred zeros form the
   * desired independent set. Otherwise, return to Step 1.
   *
   */

  // 1. Z0: uncovered 0'
  auto [row, col] = uncovered_primed_zero_;

  assert(row != -1 && col != -1);

  std::vector<std::pair<int, int>> sequence = {{row, col}};

  while (true) { // O(dim_)

    // 2. Z1: 0* in Z0's column
    auto starred_row = starred_zero_in_col(sequence.back().second); // O(dim_)
    if (starred_row == -1) {
      break;
    }

    sequence.push_back({starred_row, sequence.back().second});

    // 3. Z2: 0' in Z1's row, (we must prove that it exists).
    auto primed_col = primed_zero_in_row(sequence.back().first); // O(dim_)
    sequence.push_back({sequence.back().first, primed_col});
  }

  row_cover_.clear();
  col_cover_.clear();
  for (const auto &[row, col] : sequence) {
    if (mask_[row][col] == kStar) {
      mask_[row][col] = 0;
    } else { // kPrime -> kStar
      mask_[row][col] = kStar;
    }
  }

  for (size_t row = 0; row < dim_; ++row) {
    for (size_t col = 0; col < dim_; ++col) {
      if (mask_[row][col] == kStar) {
        col_cover_.insert(col);
      }
    }
  }

  if (row_cover_.size() == dim_ || col_cover_.size() == dim_) {
    next_step_ = -1;
    return;
  }

  next_step_ = 3;
  return;
}

void HungarianAssigner::step3() {
  /**
   * @brief step 3: make new zeros uncovered
   *
   * Let h denote the smallest non-covered element of the matrix; it will be
   * positive. Add h to each covered row; then subtract h from each uncovered
   * column.
   *
   * Return to Step 1, without altering any asterisks, primes, or covered lines.
   *
   */

  auto h = min_uncovered(); // O(dim_^2)

  // O(dim_^2)
  for (size_t row = 0; row < dim_; row++) {
    for (size_t col = 0; col < dim_; col++) {
      if (row_cover_.count(row) > 0) {
        opt_matrix_[row][col] += h;
      }
      if (col_cover_.count(col) == 0) {
        opt_matrix_[row][col] -= h;
      }
    }
  }

  next_step_ = 1;
}

float HungarianAssigner::wrapUp(std::vector<int> *assignment_index) {
  assignment_index->clear();
  assignment_index->resize(n_, -1);

  float cost = 0.f;

  for (size_t row = 0; row < n_; ++row) {
    for (size_t col = 0; col < m_; ++col) {
      if (mask_[row][col] == kStar) {
        assignment_index->at(row) = col;
        cost += cost_matrix_[row][col];
      }
    }
  }

  return cost;
}

int HungarianAssigner::starred_zero_in_row(const int row) {
  /**
   * @brief return column of star in row
   * @param row
   * @return column of star in row
   * @return -1 if no star in row
   */

  for (size_t col = 0; col < dim_; ++col) {
    if (mask_[row][col] == kStar) {
      return col;
    }
  }
  return -1;
}

int HungarianAssigner::starred_zero_in_col(const int col) {
  for (size_t row = 0; row < dim_; ++row) {
    if (mask_[row][col] == kStar) {
      return row;
    }
  }
  return -1;
}

int HungarianAssigner::primed_zero_in_row(const int row) {
  for (size_t col = 0; col < dim_; ++col) {
    if (mask_[row][col] == kPrime) {
      return col;
    }
  }
  return -1;
}

std::pair<int, int> HungarianAssigner::uncovered_zero() {
  /**
   * @brief find uncovered zero
   * @return row of uncovered zero
   * @return -1 if no uncovered zero
   */
  for (size_t row = 0; row < dim_; ++row) {
    for (size_t col = 0; col < dim_; ++col) {
      if (std::abs(opt_matrix_[row][col]) <= EPSILON_FLOAT &&
          row_cover_.count(row) == 0 && col_cover_.count(col) == 0) {
        return {row, col};
      }
    }
  }
  return {-1, -1};
}

float HungarianAssigner::min_uncovered() {
  float min_val = std::numeric_limits<float>::max();

  for (size_t row = 0; row < dim_; ++row) {
    for (size_t col = 0; col < dim_; ++col) {
      if (row_cover_.count(row) == 0 && col_cover_.count(col) == 0) {
        if (opt_matrix_[row][col] < min_val) {
          min_val = opt_matrix_[row][col];
        }
      }
    }
  }
  return min_val;
}
