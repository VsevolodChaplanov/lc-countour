#ifndef PDE_THOMAS_SOLVER_HPP
#define PDE_THOMAS_SOLVER_HPP

#include "tridiagonal_matrix.hpp"

template<std::floating_point T, std::size_t n>
class ThomasSolver final {
public:
  ThomasSolver(const TridiagonalMatrix<T, n>& lhs)
    : lhs_{lhs} { }

  template<std::ranges::random_access_range Rhs>
  auto solve(Rhs&& rhs) {
    using Diagonal = TridiagonalMatrix<T, n>::Diagonal;
    // direct
    Vector<T, n> p, q;

    p[0] = lhs_.get_value(0, 1) / lhs_.get_value(0, 0);
    q[0] = rhs[0] / lhs_.get_value(0, 0);

    for (const std::size_t i : std::views::iota(0ull, n)) {
      p[i] = - lhs_.get_value(i, i + 1) / (lhs_.get_value(i, i) + lhs_.get_value(i, i - 1) * p [i - 1]);
      q[i] = (- lhs_.get_value(i, i - 1) * q[i - 1] + rhs[i]) / (lhs_.get_value(i, i) + lhs_.get_value(i, i - 1) * p[i-1]);
    }

    // reverse
    solution_[n - 1] = (- lhs_.get_value(n - 1, n - 2) * q[n - 2] + rhs[n - 1])
        / (lhs_.get_value(n - 1,n - 1) + p[n - 2] * lhs_.get_value(n - 1, n - 2));

    for (const std::size_t i : std::views::iota(0ull, n - 1) | std::views::reverse) {
      solution_[i - 1] = p[i - 1] * (solution_[i]) + q[i - 1];
    }

    return solution_ | std::views::all;
  }

private:
  const TridiagonalMatrix<T, n>& lhs_;
  Vector<T, n> solution_;
};

#endif //PDE_THOMAS_SOLVER_HPP
