#ifndef PDE_TRIDIAGONAL_MATRIX_HPP
#define PDE_TRIDIAGONAL_MATRIX_HPP

#include <concepts>
#include <numeric>
#include <algorithm>
#include <ranges>
#include <array>
#include <cassert>

template<std::floating_point T, std::size_t n>
using Vector = std::array<T, n>;

template<std::floating_point T, std::size_t n>
class TridiagonalMatrix final {
  static_assert(n >= 1);
public:
  enum struct Diagonal {
    Upper,
    Main,
    Lower
  };

  using value_type = T;

  template<std::ranges::viewable_range Range>
  constexpr void fill(Range&& values, Diagonal diagonal) {
    switch (diagonal) {
      case Diagonal::Upper:
        assert(static_cast<std::size_t>(std::ranges::distance(values)) == n - 1);
        std::copy(std::ranges::begin(values), std::ranges::end(values), a_.begin());
        break;
      case Diagonal::Main:
        assert(static_cast<std::size_t>(std::ranges::distance(values)) == n);
        std::copy(std::ranges::begin(values), std::ranges::end(values), b_.begin());
        break;
      case Diagonal::Lower:
        assert(static_cast<std::size_t>(std::ranges::distance(values)) == n - 1);
        std::copy(std::ranges::begin(values), std::ranges::end(values), c_.begin());
        break;
    }
  }

  void set_value(value_type value, Diagonal diagonal, std::size_t place) {
    switch (diagonal) {
      case Diagonal::Upper:
        assert(place <= n - 1);
        a_[place] = value;
        break;
      case Diagonal::Main:
        assert(place <= n - 1);
        b_[place] = value;
        break;
      case Diagonal::Lower:
        assert(place <= n - 1);
        c_[place] = value;
        break;
    }
  }

  void set_value(value_type value, std::size_t i, std::size_t j) {
    if ( i == j) b_[i] = value;
    else if (i > j && i - j == 1) c_[i] = value;
    else if (i < j && j - i == 1) a_[i] = value;
  }

  void add_value(value_type value, Diagonal diagonal, std::size_t place) {
    const auto prev = get_value(diagonal, place);
    const auto fut = value + prev;
    set_value(fut, diagonal, place);
  }

  constexpr value_type get_value(Diagonal diagonal, std::size_t place) const {
    switch (diagonal) {
      case Diagonal::Upper:
        assert(place <= n - 1);
        return a_[place];
      case Diagonal::Main:
        assert(place <= n - 1);
        return a_[place];
      case Diagonal::Lower:
        assert(place <= n - 1);
        return a_[place];
    }
  }

  constexpr value_type get_value(std::size_t i, std::size_t j) const {
/*    assert(i <= n);
    assert(j <= n);*/

    if ( i == j) return b_[i];
    else if (i > j && i - j == 1) return c_[i];
    else if (i < j && j - i == 1) return a_[i];
    else return 0;
  }

  constexpr std::size_t dim() const { return n; }

protected:
  std::array<value_type, n - 1> a_;
  std::array<value_type, n> b_;
  std::array<value_type, n - 1> c_;
};

#endif //PDE_TRIDIAGONAL_MATRIX_HPP
