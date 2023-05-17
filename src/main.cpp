#include <iostream>
#include <fstream>
#include "approximator.hpp"
#include "thomas_solver.hpp"
#include "tridiagonal_matrix.hpp"

namespace {
  constexpr double t_init = 0;
  constexpr double t_end = 100;
  constexpr std::size_t n = 100000;
  constexpr double tau = (t_end - t_init) / (n - 1);
}

int main() try {
  TridiagonalMatrix<double, n> lhs;
  Vector<double, n> rhs;

  // d2 q + h^2 / LC q = 0;
  Approximate::AddApproximationSecondDerivative(lhs);
  Approximate::AddApproximationStandalone(lhs,
    std::views::iota(0ull, n) | std::views::transform([](auto val) {
    return tau * tau / (inductance * capacity);
  }));
  Approximate::SetTimeBoundaryCondition(lhs, rhs, capacity * voltage);

  ThomasSolver solver{lhs};
  const auto solution = solver.solve(rhs);

  std::ofstream file{"solution"};
  std::ranges::for_each(solution, [&](auto val) {
    file << val << '\n';
  });

} catch (const std::exception& ex) {
  std::cerr << "catch: " << ex.what() << '\n';
}