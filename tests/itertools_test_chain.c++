#include <catch2/catch_approx.hpp>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <list>

#include "base.h"
#include "chain.h"
#include "testtools.h"

namespace it = ATATIteratorTools;
namespace tt = ATATTestTools;
using Vec = const std::vector<int>;

TEST_CASE("chain: three strings", "[chain]") {
  std::string s1{"abc"};
  std::string s2{"mno"};
  std::string s3{"xyz"};
  auto ch = it::chain(s1, s2, s3);

  Vec v(std::begin(ch), std::end(ch));
  Vec vc{'a', 'b', 'c', 'm', 'n', 'o', 'x', 'y', 'z'};

  REQUIRE(v == vc);
}

TEST_CASE("chain: const iteration", "[chain][const]") {
  std::string s1{"abc"};
  const std::string s2{"mno"};
  const auto ch = it::chain(s1, s2, std::string{"xyz"});

  Vec v(std::begin(ch), std::end(ch));
  Vec vc{'a', 'b', 'c', 'm', 'n', 'o', 'x', 'y', 'z'};

  REQUIRE(v == vc);
}

// TODO make this work
/*
TEST_CASE("chain: const iterators can be compared to non-const itertors",
	  "[chain][const]") {
  auto ch = it::chain(std::string{}, std::string{});
  const auto& cch = ch;
  (void)(std::begin(ch) == std::end(cch));
}*/

TEST_CASE("chain: with different container types", "[chain]") {
  std::string s1{"abc"};
  std::list<char> li{'m', 'n', 'o'};
  std::vector<char> vec{'x', 'y', 'z'};
  auto ch = it::chain(s1, li, vec);

  Vec v(std::begin(ch), std::end(ch));
  Vec vc{'a', 'b', 'c', 'm', 'n', 'o', 'x', 'y', 'z'};

  REQUIRE(v == vc);
}

TEST_CASE("chain: handles empty containers", "[chain]") {
  std::string emp;
  std::string a{"a"};
  std::string b{"b"};
  std::string c{"c"};
  Vec vc{'a', 'b', 'c'};

  SECTION("Empty container at front") {
    auto ch = it::chain(emp, a, b, c);
    Vec v(std::begin(ch), std::end(ch));

    REQUIRE(v == vc);
  }

  SECTION("Empty container at back") {
    auto ch = it::chain(a, b, c, emp);
    Vec v(std::begin(ch), std::end(ch));

    REQUIRE(v == vc);
  }

  SECTION("Empty container in middle") {
    auto ch = it::chain(a, emp, b, emp, c);
    Vec v(std::begin(ch), std::end(ch));

    REQUIRE(v == vc);
  }

  SECTION("Consecutive empty containers at front") {
    auto ch = it::chain(emp, emp, a, b, c);
    Vec v(std::begin(ch), std::end(ch));

    REQUIRE(v == vc);
  }

  SECTION("Consecutive empty containers at back") {
    auto ch = it::chain(a, b, c, emp, emp);
    Vec v(std::begin(ch), std::end(ch));

    REQUIRE(v == vc);
  }

  SECTION("Consecutive empty containers in middle") {
    auto ch = it::chain(a, emp, emp, b, emp, emp, c);
    Vec v(std::begin(ch), std::end(ch));

    REQUIRE(v == vc);
  }
}

TEST_CASE("chain: with only empty containers", "[chain]") {
  std::string emp{};
  SECTION("one empty container") {
    auto ch = it::chain(emp);
    REQUIRE_FALSE(std::begin(ch) != std::end(ch));
  }

  SECTION("two empty containers") {
    auto ch = it::chain(emp, emp);
    REQUIRE_FALSE(std::begin(ch) != std::end(ch));
  }

  SECTION("three empty containers") {
    auto ch = it::chain(emp, emp, emp);
    REQUIRE_FALSE(std::begin(ch) != std::end(ch));
  }
}

TEST_CASE("chain: doesn't move or copy elements of iterable", "[chain]") {
  constexpr tt::MonolithObject<int> arr[] = {{6}, {7}, {8}};
  for (auto&& i : it::chain(arr, arr)) {
    (void)i;
  }
}

TEST_CASE("chain: binds reference to lvalue and moves rvalue", "[chain]") {
  tt::IterableType<char> bi{'x', 'y', 'z'};
  tt::IterableType<char> bi2{'a', 'j', 'm'};
  SECTION("First moved, second ref'd") {
    it::chain(std::move(bi), bi2);
    REQUIRE(bi.was_moved_from());
    REQUIRE_FALSE(bi2.was_moved_from());
  }
  SECTION("First ref'd, second moved") {
    it::chain(bi, std::move(bi2));
    REQUIRE_FALSE(bi.was_moved_from());
    REQUIRE(bi2.was_moved_from());
  }
}

TEST_CASE("chain: operator==", "[chain]") {
  std::string emp{};
  auto ch = it::chain(emp);
  REQUIRE(std::begin(ch) == std::end(ch));
}

TEST_CASE("chain: postfix ++", "[chain]") {
  std::string s1{"a"}, s2{"b"};
  auto ch = it::chain(s1, s2);
  auto it = std::begin(ch);
  it++;
  REQUIRE(*it == 'b');
}

TEST_CASE("chain: iterator meets requirements", "[chain]") {
  Vec ns{};
  auto c = it::chain(ns, ns);
  REQUIRE(it::is_iterator_v<decltype(std::begin(c))>);
}

template <typename... Ts>
using Imp_T = decltype(it::chain(std::declval<Ts>()...));

TEST_CASE("chain: has correct ctor and assign ops", "[chain]") {
  using T = Imp_T<std::string&, std::vector<char>, char(&)[10]>;
  REQUIRE(it::is_move_constructible_only<T>::value);
}
