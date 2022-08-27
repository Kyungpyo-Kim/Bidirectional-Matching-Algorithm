/**
 * @file hungarian_algorithm.cc
 * @author kyungpyo-kim (kyungpyo94@gmail.com)
 * @date 2022-08-24
 * @brief Hungarian algorithm
 */

#include "hungarian_assigner.hpp"
#include <algorithm>
#include <cassert>
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
  is_working_ = true;
  next_step_ = 0;

  while (is_working_) {
    switch (next_step_) {
    case 0:
      preliminaries(); // O(dim_^2)
      show("preliminaries");
      break;

    case 1:
      step1(); // O(dim_^3)
      show("step1");
      break;

    case 2:
      step2(); // O(dim_^2)
      show("step2");
      break;

    case 3:
      step3(); // O(dim_^2)
      show("step3");
      break;

    default:
      break;
    }

    if (next_step_ == -1) {
      break;
    }
  }

  return wrapUp(assignment_index);
}

void HungarianAssigner::show(std::string name = "") {
  std::cout << "[ " << name << " ]" << std::endl;
  std::cout << "matrix: " << std::endl;
  for (size_t row = 0; row < dim_; ++row) {
    for (size_t col = 0; col < dim_; ++col) {
      std::cout << opt_matrix_[row][col] << " ";
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
  for (auto &row : opt_matrix_) {
    auto min_value = *std::min_element(row.begin(), row.end());
    for (auto &elem : row) {
      elem -= min_value;
    }
  }

  std::unordered_set<size_t> row_check;
  std::unordered_set<size_t> col_check;
  for (size_t row = 0; row < dim_; ++row) {
    for (size_t col = 0; col < dim_; ++col) {
      if (opt_matrix_[row][col] == 0 && row_check.count(row) == 0 &&
          col_check.count(col) == 0) {
        // mark starred zeros
        mask_[row][col] = kStar;
        row_check.insert(row);
        col_check.insert(col);

        // cover each column containing a starred zero
        col_cover_.insert(col);
      }
    }
  }

  next_step_ = 1;
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
      next_step_ = 2;
      return;
    }

    row_cover_.insert(row);
    col_cover_.erase(starred_zero);
  }

  if (row_cover_.size() == dim_) {
    next_step_ = -1;
    return;
  }
}

void HungarianAssigner::step2() {
  /**
   * @brief
   *
   * There is a sequence of alternating starred and primed zeros, constructed as
   * follows: Let Zo denote the uncovered 0'. [There is only one.] Let Z1 denote
   * the 0* in Z0's column (if any). Let Z2 denote the 0' in Z1's row (we must
   * prove that it exists). Let Z3 denote the 0* in Z2's column (if any).
   * Similarly continue until the sequence stops at a 0', Z2k, which has no 0*
   * in its column (this we must also prove). [Note that no column contains more
   * than one 0* and no row more than one 0', so that the sequence is uniquely
   * specified. The sequence may, however, contain only one element. Now the
   * column of Z1 is not covered, so its row must be covered; hence there is a
   * 0' in this row (see Step 1). This 0' serves as our Z2 . A similar argument
   * applies to show that, given Z2i-1, Z2j exists. Now let us index the primed
   * zeros 1, 2, 3, *.. in the order in which we primed them during Step 1. One
   * sees readily from the directions in Step 1 that the index of Z2i must be
   * smaller than the index of Z2i-2 It follows at once that the sequence does
   * stop and, furthermore, that all the elements of the sequence Z ... , Z2k
   * are distinct elements of the matrix.] Unstar each starred zero of the
   * sequence and star each primed zero of the sequence. [The resulting set of
   * starred zeros is easily seen to be independent. It is larger by one than
   * the previous set of independent starred zeros.] Erase all primes, uncover
   * every row, and cover every column containing a 0*. If all columns are
   * covered, the starred zeros form the desired independent set. Otherwise,
   * return to Step 1.
   *
   */

  // 1. Z0: uncovered 0'
  auto [row, col] = uncovered_primed_zero();
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
    auto primed_col = primed_zero_in_row(sequence.back().first);
    sequence.push_back({sequence.back().first, primed_col});
  }

  row_cover_.clear();
  col_cover_.clear();
  for (const auto &[row, col] : sequence) {
    if (mask_[row][col] == kStar) {
      mask_[row][col] = 0;
    } else { // kPrime -> kStar
      mask_[row][col] = kStar;
      col_cover_.insert(col);
    }
  }

  if (col_cover_.size() == dim_) {
    next_step_ = -1;
    return;
  }

  next_step_ = 3;
  return;
}

void HungarianAssigner::step3() {

  /**
   * @brief
   *
   * Let h denote the smallest non-covered element of the matrix; it will be
   * positive. Add h to each covered row; then subtract h from each uncovered
   * column.
   *
   * Return to Step 1, without altering any asterisks, primes, or covered lines.
   * [One might think one should "erase all primes, uncover every row, and cover
   * every column containing a 0*" before returning to Step 1, so that the input
   * of Step 1 is the standard one. That this is unnecessary may be seen from
   * the following argument: The effect of the transformation specified above is
   * to decrease each non-covered element of the matrix by h, increase each
   * twice-covered element by h, and leave each once-covered element unaltered.
   * Each 0* and 0' is once-covered, so each is still a zero This content
   * downloaded from 147.8.31.43 on Sun, 06 Dec 2015 11:06:56 UTC All use
   * subject to JSTOR Terms and Conditions ASSIGNMENT AND TRANSPORTATION
   * PROBLEMS 35 after the transformation. (Incidentally, this shows that nk+i _
   * nk , where ni denotes the maximal number of independent zeros of Ai. Ak
   * denotes the matrix before the transformation and Ak?1 denotes the
   * transformed matrix.) Let us index the primed zeros 1, 2, 3, ... in the
   * order in which they were primed previously. Imagine that we "erase all
   * primes, uncover every row, and cover every column containing a 0*" before
   * returning to Step 1. The directions in Step 1 tell us to find a non-covered
   * zero, but do not specify which one of possibly several non-covered ones we
   * should choose to prime. Hence we can consider the indexed zeros in the
   * order of their indices, priming each one in turn, without violating the
   * directions in Step 1. This will bring the configuration of asterisks,
   * primes, and covered lines back to precisely the same one we had at the
   * beginning of this paragraph. Hence there is no need to return to the
   * standard input when passing from Step 3 to Step 1; and indeed such a return
   * would be foolish.]
   *
   */

  auto h = min_uncovered();

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
      if (opt_matrix_[row][col] == 0 && row_cover_.count(row) == 0 &&
          col_cover_.count(col) == 0) {
        return {row, col};
      }
    }
  }
  return {-1, -1};
}

std::pair<int, int> HungarianAssigner::uncovered_primed_zero() {
  for (size_t row = 0; row < dim_; ++row) {
    for (size_t col = 0; col < dim_; ++col) {
      if (mask_[row][col] == kPrime && row_cover_.count(row) == 0 &&
          col_cover_.count(col) == 0) {
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
