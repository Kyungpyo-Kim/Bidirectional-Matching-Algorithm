/**
 * @file hungarian_algorithm.hpp
 * @author kyungpyo-kim (kyungpyo94@gmail.com)
 * @date 2022-08-24
 * @brief Hungarian algorithm
 */

#include <vector>
#include <limits>

class HungarianAssigner {
  using CostType = std::vector<std::vector<float>>;

public:
  float solve(const CostType &cost_matrix, const size_t n, const size_t m,
              const int mode, std::vector<size_t> *assignment_index);

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
  float epsilon = std::numeric_limits<float>::epsilon();
  size_t n_;
  size_t m_;
  int mode_;
  size_t dim_;
  CostType cost_matrix_;
};
