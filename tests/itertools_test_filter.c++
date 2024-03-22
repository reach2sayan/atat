#include <catch2/catch_approx.hpp>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>

#include "base.h"
#include "filter.h"
#include "testtools.h"

namespace it = ATATIteratorTools;
namespace tt = ATATTestTools;

using Vec = const std::vector<int>;

bool less_than_five(int i) { return i < 5; }
class LessThanValue {
 private:
  int compare_val;

 public:
  LessThanValue(int v) : compare_val(v) {}
  bool operator()(int i) { return i < this->compare_val; }
};

TEST_CASE("filter: handles different callable types", "[filter]") {
  Vec ns = {1, 2, 5, 6, 3, 1, 7, -1, 5};
  Vec vc = {1, 2, 3, 1, -1};
  SECTION("with function pointer") {
    auto f = it::filter(less_than_five, ns);
    Vec v(std::begin(f), std::end(f));
    REQUIRE(v == vc);
  }

  SECTION("with callable object") {
    auto f = it::filter(LessThanValue{5}, ns);
    Vec v(std::begin(f), std::end(f));
    REQUIRE(v == vc);
  }

  SECTION("with lambda") {
    auto ltf = [](int i) { return i < 5; };
    auto f = it::filter(ltf, ns);
    Vec v(std::begin(f), std::end(f));
    REQUIRE(v == vc);
  }
}

/*
TEST_CASE("filter: handles pointer to member", "[filter]") {
  using itertest::Point;
  const std::vector<Point> ps = {{0, 3}, {4, 0}, {0, 1}, {-1, -1}};
  std::vector<Point> v;
  SECTION("with pointer to data member") {
    auto f = filter(&Point::x, ps);
    v.assign(std::begin(f), std::end(f));
  }

  SECTION("with pointer to member function") {
    auto f = filter(&Point::get_x, ps);
    v.assign(std::begin(f), std::end(f));
  }

  const std::vector<Point> vc = {{4, 0}, {-1, -1}};
  REQUIRE(v == vc);
}
*/

TEST_CASE("filter: const iteration", "[filter][const]") {
  Vec ns = {1, 2, 5, 6, 3, 1, 7, -1, 5};
  const auto f = it::filter(LessThanValue{5}, ns);
  Vec v(std::begin(f), std::end(f));
  Vec vc = {1, 2, 3, 1, -1};
  REQUIRE(v == vc);
}

TEST_CASE("filter: const iterator can be compared to non-const iterator",
	  "[filter][const]") {
  auto f = it::filter(LessThanValue{5}, Vec{});
  const auto& cf = f;
  (void)(std::begin(f) == std::end(cf));
}

TEST_CASE("filter: iterator with lambda can be assigned", "[filter]") {
  Vec ns{};
  auto ltf = [](int i) { return i < 5; };
  auto f = it::filter(ltf, ns);
  auto it = std::begin(f);
  it = std::begin(f);
}

TEST_CASE("filter: using identity", "[filter]") {
  Vec ns{0, 1, 2, 0, 3, 0, 0, 0, 4, 5, 0};
  auto f = it::filter(ns);
  Vec v(std::begin(f), std::end(f));

  Vec vc = {1, 2, 3, 4, 5};
  REQUIRE(v == vc);
}

TEST_CASE("filter: skips null pointers", "[filter]") {
  int a = 1;
  int b = 2;
  const std::vector<int*> ns = {0, &a, nullptr, nullptr, &b, nullptr};

  auto f = it::filter(ns);
  const std::vector<int*> v(std::begin(f), std::end(f));
  const std::vector<int*> vc = {&a, &b};
  REQUIRE(v == vc);
}

TEST_CASE("filter: binds to lvalues, moves rvales", "[filter]") {
  tt::IterableType<int> bi{1, 2, 3, 4};

  SECTION("one-arg binds to lvalues") {
    it::filter(bi);
    REQUIRE_FALSE(bi.was_moved_from());
  }

  SECTION("two-arg binds to lvalues") {
    it::filter(less_than_five, bi);
    REQUIRE_FALSE(bi.was_moved_from());
  }

  SECTION("one-arg moves rvalues") {
    it::filter(std::move(bi));
    REQUIRE(bi.was_moved_from());
  }

  SECTION("two-arg moves rvalues") {
    it::filter(less_than_five, std::move(bi));
    REQUIRE(bi.was_moved_from());
  }
}

TEST_CASE("filter: operator->", "[filter]") {
  std::vector<std::string> vs = {"ab", "abc", "abcdef"};
  auto f =
      it::filter([](const std::string& str) { return str.size() > 4; }, vs);
  auto it = std::begin(f);
  REQUIRE(it->size() == 6);
}

TEST_CASE("filter: all elements fail predicate", "[filter]") {
  Vec ns{10, 20, 30, 40, 50};
  auto f = it::filter(less_than_five, ns);

  SECTION("normal compare") { REQUIRE(std::begin(f) == std::end(f)); }
  SECTION("reversed compare") { REQUIRE(std::end(f) == std::begin(f)); }
}

TEST_CASE("filter: doesn't move or copy elements of iterable", "[filter]") {
  constexpr tt::MonolithObject<int> arr[] = {{1}, {0}, {2}};
  for (auto&& i : it::filter(
	   [](const tt::MonolithObject<int>& si) { return si.get(); }, arr)) {
    (void)i;
  }
}

TEST_CASE("filter: iterator meets requirements", "[filter]") {
  std::string s{};
  auto c = it::filter([] { return true; }, s);
  REQUIRE(it::is_iterator_v<decltype(std::begin(c))>);
}

template <typename T, typename U>
using Imp_t = decltype(it::filter(std::declval<T>(), std::declval<U>()));
TEST_CASE("filter: has correct ctor and assign ops", "[filter]") {
  using T1 = Imp_t<bool (*)(char c), std::string&>;
  auto lam = [](char) { return false; };
  using T2 = Imp_t<decltype(lam), std::string>;
  REQUIRE(it::is_move_constructible_only<T1>::value);
  REQUIRE(it::is_move_constructible_only<T2>::value);
}

