#include <catch2/catch_approx.hpp>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>

#include "base.h"
#include "imap.h"
#include "testtools.h"

namespace it = ATATIteratorTools;
namespace tt = ATATTestTools;

using Vec = const std::vector<int>;

namespace {
int plusone(int i) { return i + 1; }

class PlusOner {
 public:
  int operator()(int i) { return i + 1; }
};

int power(int b, int e) {
  int acc = 1;
  for (int i = 0; i < e; ++i) {
    acc *= b;
  }
  return acc;
}
}  // namespace

TEST_CASE("imap: works with lambda, callable, and function", "[imap]") {
  Vec ns = {10, 20, 30};
  std::vector<int> v;
  SECTION("with lambda") {
    auto im = it::imap([](int i) { return i + 1; }, ns);
    v.assign(std::begin(im), std::end(im));
  }

  SECTION("with callable") {
    SECTION("Normal call") {
      auto im = it::imap(PlusOner{}, ns);
      v.assign(std::begin(im), std::end(im));
    }
  }

  SECTION("with function") {
    auto im = it::imap(PlusOner{}, ns);
    v.assign(std::begin(im), std::end(im));
  }

  Vec vc = {11, 21, 31};
  REQUIRE(v == vc);
}

TEST_CASE("imap: supports const iteration", "[imap][const]") {
  Vec ns = {10, 20, 30};
  const auto m = it::imap(PlusOner{}, ns);
  Vec v(std::begin(m), std::end(m));
  Vec vc = {11, 21, 31};
  REQUIRE(v == vc);
}

TEST_CASE("imap: const iterators can be compared to non-const iterators",
	  "[imap][const]") {
  auto m = it::imap(PlusOner{}, Vec{});
  const auto& cm = m;
  (void)(std::begin(m) == std::end(cm));
}

TEST_CASE("imap: operator->", "[imap]") {
  std::vector<std::string> vs = {"ab", "abcd", "abcdefg"};
  {
    auto m = it::imap([](std::string& s) { return s; }, vs);
    auto it = std::begin(m);
    REQUIRE(it->size() == 2);
  }

  {
    auto m = it::imap([](std::string& s) -> std::string& { return s; }, vs);
    auto it = std::begin(m);
    REQUIRE(it->size() == 2);
  }
}

TEST_CASE("imap: empty sequence gives nothing", "[imap]") {
  Vec v{};
  auto im = it::imap(plusone, v);
  REQUIRE(std::begin(im) == std::end(im));
}

TEST_CASE("imap: binds to lvalues, moves rvalues", "[imap]") {
  tt::IterableType<int> bi{1, 2};
  SECTION("binds to lvalues") {
    it::imap(plusone, bi);
    REQUIRE_FALSE(bi.was_moved_from());
  }

  SECTION("moves rvalues") {
    it::imap(plusone, std::move(bi));
    REQUIRE(bi.was_moved_from());
  }
}

TEST_CASE("imap: doesn't move or copy elements of iterable", "[imap]") {
  constexpr tt::MonolithObject<int> arr[] = {{1}, {0}, {2}};
  for (auto&& i : it::imap(
	   [](const tt::MonolithObject<int>& si) { return si.get(); }, arr)) {
    (void)i;
  }
}

TEST_CASE("imap: postfix ++", "[imap]") {
  Vec ns = {10, 20};
  auto im = it::imap(plusone, ns);
  auto it = std::begin(im);
  it++;
  REQUIRE((*it) == 21);
  it++;
  REQUIRE(it == std::end(im));
}

TEST_CASE("imap: iterator meets requirements", "[imap]") {
  std::string s{};
  auto c = it::imap([](char) { return 1; }, s);
  REQUIRE(it::is_iterator_v<decltype(std::begin(c))>);
}

template <typename T, typename U>
using Imp_t = decltype(it::imap(std::declval<T>(), std::declval<U>()));
TEST_CASE("imap: has correct ctor and assign ops", "[imap]") {
  using T1 = Imp_t<bool (*)(char c), std::string&>;
  auto lam = [](char) { return false; };
  using T2 = Imp_t<decltype(lam), std::string>;
  REQUIRE(it::is_move_constructible_only<T1>::value);
  REQUIRE(it::is_move_constructible_only<T2>::value);
}
