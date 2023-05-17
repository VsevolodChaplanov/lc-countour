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

  Approximate::AddApproximationSecondDerivative(lhs);
  Approximate::AddApproximationStandalone(lhs,
                                          std::views::iota(0ull, n) | std::views::transform([](auto val) {
                                            return 1.;
                                          }));
  Approximate::SetTimeBoundaryCondition(lhs, rhs, 1.);

  ThomasSolver solver{lhs};
  const auto solution = solver.solve(rhs);

  std::ofstream file{"solution"};
  std::ranges::for_each(solution, [&](auto val) {
    file << val << '\n';
  });

  const auto exact = [] (double time) { return std::cos(time); };

  for (const std::size_t i : std::views::iota(0ull, n)) {
    auto time = t_init + i * tau;
    std::cout << solution[i] - exact(time) << std::endl;
  }

} catch (const std::exception& ex) {
  std::cerr << "catch: " << ex.what() << '\n';
}
