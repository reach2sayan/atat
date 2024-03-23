#include <catch2/catch_approx.hpp>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>

#include "base.h"
#include "testtools.h"
#include "zip.h"

namespace it = ATATIteratorTools;
namespace tt = ATATTestTools;

TEST_CASE("zip: Simple case, same length", "[zip]") {
  using Tu = std::tuple<int, char, double>;
  using ResVec = const std::vector<Tu>;
  std::vector<int> iv{10, 20, 30};
  std::string s{"hey"};
  double arr[] = {1.0, 2.0, 4.0};

  auto z = it::zip(iv, s, arr);
  ResVec v(std::begin(z), std::end(z));
  ResVec vc{Tu{10, 'h', 1.0}, Tu{20, 'e', 2.0}, Tu{30, 'y', 4.0}};
  REQUIRE(v == vc);
}  // namespace )

TEST_CASE("zip: const iteration", "[zip][const]") {
  using Tu = std::tuple<int, char, double>;
  using ResVec = const std::vector<Tu>;
  std::vector<int> iv{10, 20, 30};
  std::string s{"hey"};
  double arr[] = {1.0, 2.0, 4.0};

  const auto z = it::zip(iv, s, arr);

  ResVec v(std::begin(z), std::end(z));
  ResVec vc{Tu{10, 'h', 1.0}, Tu{20, 'e', 2.0}, Tu{30, 'y', 4.0}};
  REQUIRE(v == vc);
}

TEST_CASE("zip: const iterators can be compared to non-const iterators",
	  "[zip][const]") {
  std::vector<int> v;
  std::string s;
  auto z = it::zip(std::vector<int>{}, s);
  const auto& cz = z;
  (void)(std::begin(z) == std::end(cz));
}

TEST_CASE("zip: One empty, all empty", "[zip]") {
  std::vector<int> iv = {1, 2, 3};
  std::string s{};
  auto z = it::zip(iv, s);
  REQUIRE_FALSE(std::begin(z) != std::end(z));
  auto z2 = it::zip(s, iv);
  REQUIRE_FALSE(std::begin(z2) != std::end(z2));
}

TEST_CASE("zip: terminates on shortest sequence", "[zip]") {
  std::vector<int> iv{1, 2, 3, 4, 5};
  std::string s{"hi"};
  auto z = it::zip(iv, s);

  REQUIRE(std::distance(std::begin(z), std::end(z)) == 2);
}

TEST_CASE("zip: Modify sequence through zip", "[zip]") {
  std::vector<int> iv{1, 2, 3};
  std::vector<int> iv2{1, 2, 3, 4};
  for (auto&& [a, b] : it::zip(iv, iv2)) {
    a = -1;
    b = -1;
  }

  const std::vector<int> vc{-1, -1, -1};
  const std::vector<int> vc2{-1, -1, -1, 4};
  REQUIRE(iv == vc);
  REQUIRE(iv2 == vc2);
}

TEST_CASE("zip: binds reference when it should", "[zip]") {
  tt::IterableType<char> bi{'x', 'y', 'z'};
  it::zip(bi);
  REQUIRE_FALSE(bi.was_moved_from());
}

TEST_CASE("zip: moves rvalues", "[zip]") {
  tt::IterableType<char> bi{'x', 'y', 'z'};
  it::zip(std::move(bi));
  REQUIRE(bi.was_moved_from());
}

TEST_CASE("zip: Can bind ref and move in single zip", "[zip]") {
  tt::IterableType<char> b1{'x', 'y', 'z'};
  tt::IterableType<char> b2{'a', 'b'};
  it::zip(b1, std::move(b2));
  REQUIRE_FALSE(b1.was_moved_from());
  REQUIRE(b2.was_moved_from());
}

TEST_CASE("zip: doesn't move or copy elements of iterable", "[zip]") {
  constexpr tt::MonolithObject<int> arr[] = {{6}, {7}, {8}};
  for (auto&& t : it::zip(arr)) {
    (void)std::get<0>(t);
  }
}

TEST_CASE("zip: postfix ++", "[zip]") {
  const std::vector<int> v = {1};
  auto z = it::zip(v);
  auto it = std::begin(z);
  it++;
  REQUIRE(it == std::end(z));
}

TEST_CASE("zip: iterator meets requirements", "[zip]") {
  std::string s{};
  auto c = it::zip(s);
  REQUIRE(it::is_iterator_v<decltype(std::begin(c))>);
  REQUIRE(tt::reference_t_matches_deref_t<decltype(std::begin(c))>::value);
  auto c2 = it::zip(s, s);
  REQUIRE(it::is_iterator_v<decltype(std::begin(c2))>);
  REQUIRE(tt::reference_t_matches_deref_t<decltype(std::begin(c2))>::value);
}

