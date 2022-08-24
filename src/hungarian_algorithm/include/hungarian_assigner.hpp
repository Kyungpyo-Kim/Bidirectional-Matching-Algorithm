/**
 * @file hungarian_algorithm.hpp
 * @author kyungpyo-kim (kyungpyo94@gmail.com)
 * @date 2022-08-24
 * @brief Hungarian algorithm
 */

#include <vector>
#include <limits>
#include <unordered_set>

class HungarianAssigner {
public:
  HungarianAssigner() = default;
  ~HungarianAssigner() = default;

public:
  using CostType = std::vector<std::vector<float>>;
  using MaskType = std::vector<std::vector<int>>;
  constexpr static int kStar = 1;
  constexpr static int kPrime = 2;

public:
  float solve(const CostType &cost_matrix, const size_t n, const size_t m,
              const int mode, std::vector<size_t> *assignment_index);
  void show(std::string name);

private:
  void step0();
  void step1();
  void step2();
  void step3();
  void step4();
  void step5();
  void step6();
  void step7();

private:
  // support step 4
  std::pair<int, int> find_zero_cost_element();
  bool star_in_row(const size_t row);
  size_t find_col_star_in_row(const size_t row);

private:
  // suport step 5
  int find_row_star_in_col(const size_t col);
  size_t find_col_prime_in_row(const size_t row);
  void augment_path();
  void clear_lines();
  void clear_primes();

private:
  float epsilon_ = std::numeric_limits<float>::epsilon();
  size_t n_;
  size_t m_;
  int mode_;
  size_t dim_;
  CostType cost_matrix_;
  CostType opt_matrix_;
  MaskType mask_;
  std::unordered_set<size_t> line_row_;
  std::unordered_set<size_t> line_col_;
  int path_row_0_; //temporary to hold the smallest uncovered value
  int path_col_0_;
  std::vector<std::vector<int>> path_;
  size_t cnt_path_;
  std::vector<size_t> results_;
private:
  bool is_working_ = true;
  int next_step_ = 0;
};
