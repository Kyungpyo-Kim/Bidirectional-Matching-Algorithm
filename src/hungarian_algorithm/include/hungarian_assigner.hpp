/**
 * @file hungarian_algorithm.hpp
 * @author kyungpyo-kim (kyungpyo94@gmail.com)
 * @date 2022-08-24
 * @brief Hungarian algorithm
 */

#include <limits>
#include <unordered_set>
#include <vector>

class HungarianAssigner {
  constexpr static float EPSILON_FLOAT = std::numeric_limits<float>::epsilon();
  constexpr static int kStar = 1;
  constexpr static int kPrime = 2;

public:
  HungarianAssigner() = default;
  ~HungarianAssigner() = default;

public:
  using CostType = std::vector<std::vector<float>>;
  using AssignmentType = std::vector<int>;
  using MaskType = std::vector<std::vector<int>>;

public:
  /**
   * @brief solve the assignment problem
   *
   * @param cost_matrix: cost matrix of the assignment problem (n x m, row
   * major)
   * @param n: number of rows of the cost matrix
   * @param m: number of columns of the cost matrix
   * @param mode: 0: minimize the total cost, 1: maximize the total cost
   * @param assignment_index: index of the assignment, -1 if unassigned
   * @return float: total cost of the assignment
   */
  float solve(const CostType &cost_matrix, const size_t n, const size_t m,
              const int mode, AssignmentType *assignment_index);

  /**
   * @brief
   *
   * @param name
   */
  void show(std::string name);

private:
  void init();
  void preliminaries();
  void step1();
  void step2();
  void step3();

private:
  float wrapUp(std::vector<int> *assignment_index);

private:
  int starred_zero_in_row(const int row);
  int starred_zero_in_col(const int col);
  int primed_zero_in_row(const int row);
  std::pair<int, int> uncovered_zero();
  float min_uncovered();

private:
  size_t n_;
  size_t m_;
  int mode_;
  size_t dim_;
  CostType cost_matrix_;

private:
  CostType opt_matrix_;
  MaskType mask_;
  std::unordered_set<size_t> row_cover_;
  std::unordered_set<size_t> col_cover_;
  std::pair<size_t, size_t> uncovered_primed_zero_;

private:
  int next_step_ = 0;
};
