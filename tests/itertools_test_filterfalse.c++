#include <catch2/catch_approx.hpp>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>

#include "base.h"
#include "filter.h"
#include "testtools.h"

namespace it = ATATIteratorTools;
namespace tt = ATATTestTools;

using Vec = const std::vector<int>;
bool filter_false_less_than_five(int i) { return i < 5; }
class FilterFalseLessThanValue {
 private:
  int compare_val;

 public:
  FilterFalseLessThanValue(int v) : compare_val(v) {}
  bool operator()(int i) { return i < this->compare_val; }
};

TEST_CASE("filterfalse: handles different callable types", "[filterfalse]") {
  Vec ns = {1, 2, 5, 6, 3, 1, 7, -1, 5};
  Vec vc = {5, 6, 7, 5};
  SECTION("with function pointer") {
    auto f = it::filter(filter_false_less_than_five, ns, true);
    Vec v(std::begin(f), std::end(f));
    REQUIRE(v == vc);
  }

  SECTION("with callable object") {
    std::vector<int> v;
    SECTION("Normal call") {
      auto f = it::filter(FilterFalseLessThanValue{5}, ns, true);
      v.assign(std::begin(f), std::end(f));
    }
    REQUIRE(v == vc);
  }

  SECTION("with lambda") {
    auto ltf = [](int i) { return i < 5; };
    auto f = it::filter(ltf, ns, true);
    Vec v(std::begin(f), std::end(f));
    REQUIRE(v == vc);
  }
}

/*
TEST_CASE("filterfalse: handles pointer to member", "[filterfalse]") {
  using itertest::Point;
  const std::vector<Point> ps = {{0, 3}, {4, 0}, {0, 1}, {-1, -1}};
  std::vector<Point> v;
  SECTION("with pointer to data member") {
    auto f = filterfalse(&Point::x, ps);
    v.assign(std::begin(f), std::end(f));
  }

  SECTION("with pointer to member function") {
    auto f = filterfalse(&Point::get_x, ps);
    v.assign(std::begin(f), std::end(f));
  }

  const std::vector<Point> vc = {{0, 3}, {0, 1}};
  REQUIRE(v == vc);
}
*/

TEST_CASE("filterfalse: const iteration", "[filterfalse][const]") {
  Vec ns = {1, 2, 5, 6, 3, 1, 7, -1, 5};
  const auto f = it::filter(FilterFalseLessThanValue{5}, ns, true);
  Vec v(std::begin(f), std::end(f));
  Vec vc = {5, 6, 7, 5};
  REQUIRE(v == vc);
}

TEST_CASE("filterfalse: const iterator and non-const iterator can be compared",
	  "[filterfalse][const]") {
  auto f = it::filter(FilterFalseLessThanValue{5}, Vec{});
  const auto& cf = f;
  (void)(std::begin(f) == std::end(cf));
}

TEST_CASE("filterfalse: using identity", "[filterfalse]") {
  Vec ns{0, 1, 2, 0, 3, 0, 0, 0, 4, 5, 0};
  std::vector<int> v;
  SECTION("Normal call") {
    auto f = it::filter(ns, true);
    v.assign(std::begin(f), std::end(f));
  }

  Vec vc = {0, 0, 0, 0, 0, 0};
  REQUIRE(v == vc);
}

TEST_CASE("filterfalse: binds to lvalues, moves rvales", "[filterfalse]") {
  tt::IterableType<int> bi{1, 2, 3, 4};

  SECTION("one-arg binds to lvalues") {
    it::filter(bi, true);
    REQUIRE_FALSE(bi.was_moved_from());
  }

  SECTION("two-arg binds to lvalues") {
    it::filter(filter_false_less_than_five, bi, true);
    REQUIRE_FALSE(bi.was_moved_from());
  }

  SECTION("one-arg moves rvalues") {
    it::filter(std::move(bi), true);
    REQUIRE(bi.was_moved_from());
  }

  SECTION("two-arg moves rvalues") {
    it::filter(filter_false_less_than_five, std::move(bi), true);
    REQUIRE(bi.was_moved_from());
  }
}

TEST_CASE("filterfalse: all elements pass predicate", "[filterfalse]") {
  Vec ns{0, 1, 2, 3, 4};
  auto f = it::filter(filter_false_less_than_five, ns, true);

  REQUIRE(std::begin(f) == std::end(f));
}

TEST_CASE("filterfalse: iterator meets requirements", "[filterfalse]") {
  std::string s{};
  auto c = it::filter([] { return true; }, s, true);
  REQUIRE(it::is_iterator_v<decltype(std::begin(c))>);
}

template <typename T, typename U>
using ImpT = decltype(it::filter(std::declval<T>(), std::declval<U>(), true));
TEST_CASE("filterfalse: has correct ctor and assign ops", "[filterfalse]") {
  using T1 = ImpT<bool (*)(char c), std::string&>;
  auto lam = [](char) { return false; };
  using T2 = ImpT<decltype(lam), std::string>;
  REQUIRE(it::is_move_constructible_only<T1>::value);
  REQUIRE(it::is_move_constructible_only<T2>::value);
}
