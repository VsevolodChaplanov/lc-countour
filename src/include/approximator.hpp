#ifndef PDE_APPROXIMATOR_HPP
#define PDE_APPROXIMATOR_HPP

#include <concepts>
#include <vector>
#include <ranges>
#include <algorithm>
#include "boundary_conditions.hpp"

namespace detail {
  template<typename Matrix>
  concept GlobalMatrix = requires (Matrix mat) {
    mat.get_value(0ull, 0ull);
    mat.set_value(0., 0ull, 0ull);
    mat.dim();
    typename Matrix::value_type;
  };
}

// МКР
class Approximate final {
public:

  template<detail::GlobalMatrix Matrix>
  static void AddApproximationSecondDerivative(Matrix& lhs) {
    using T = Matrix::value_type;
    const auto size = lhs.dim();
    {
      std::vector<T> mass(size - 1, 1.);
      lhs.fill(mass | std::views::all, Matrix::Diagonal::Upper);
      lhs.fill(mass | std::views::all, Matrix::Diagonal::Lower);
    }
    {
      std::vector<T> mass(size, - 2.);
      lhs.fill(mass | std::views::all, Matrix::Diagonal::Main);
    }
  }

  template<detail::GlobalMatrix Matrix, std::ranges::viewable_range Range>
  static void AddApproximationStandalone(Matrix& lhs, Range&& values) {
    std::ranges::for_each(std::forward<Range>(values), [&lhs] (auto value) {
      static std::size_t index = 0;
      lhs.add_value(value, Matrix::Diagonal::Main, index++);
    });
  }

  template<detail::GlobalMatrix Matrix, typename Rhs, std::floating_point InitValue>
  static void SetTimeBoundaryCondition(Matrix& lhs, Rhs& rhs, InitValue value) {
    lhs.set_value(0, 0, 1);
    lhs.set_value(1, 0, 0);
    rhs[0] = value;
  }

private:
};

#endif //PDE_APPROXIMATOR_HPP
