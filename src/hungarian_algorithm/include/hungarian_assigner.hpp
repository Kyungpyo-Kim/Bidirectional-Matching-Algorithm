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
  using PrimeType = std::vector<std::vector<bool>>;
  using StarType = std::vector<std::vector<bool>>;

public:
  float solve(const CostType &cost_matrix, const size_t n, const size_t m,
              const int mode, std::vector<size_t> *assignment_index);
  void show();

private:
  void build_cost_matrix(const CostType &cost_matrix);

private:
  void step1();
  void step2();
  void step3();
  void step4();
  void step5();
  void step6();

private:
  float epsilon_ = std::numeric_limits<float>::epsilon();
  size_t n_;
  size_t m_;
  int mode_;
  size_t dim_;
  CostType cost_matrix_;
  CostType opt_matrix_;
  PrimeType prime_;
  StarType star_;
  std::unordered_set<size_t> line_row_;
  std::unordered_set<size_t> line_col_;
};
