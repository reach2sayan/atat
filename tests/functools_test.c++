/*
 * atft C++ higher order convenience functions
 *
 * Copyright © Björn Fahller 2017,2018
 *
 *  Use, modification and distribution is subject to the
 *  Boost Software License, Version 1.0. (See accompanying
 *  file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 *
 * Project home: https://github.com/rollbear/atft
 */

#include "functools/functools.h"

#include <catch2/catch_test_macros.hpp>
#include <ostream>
namespace atft = ATATFunctionalTools;

// constexpr tests

template <auto N>
constexpr auto eq = atft::equal(N);

template <auto N>
constexpr auto gt = atft::greater_than(N);

template <auto N>
constexpr auto ne = atft::not_equal(N);

template <typename T, typename = std::enable_if_t<std::is_same<T, bool>{}>>
bool func(T);

TEST_CASE("constexpr tests") {
  static_assert(eq<3>(3), "equal is constexpr");
  static_assert(ne<3>(1), "not_equal is constexpr");
  static_assert(atft::less_than(3)(2), "less_than is constexpr");
  static_assert(atft::less_than_equal(3)(2), "less_equal is constexpr");
  static_assert(gt<3>(4), "greater_than is constexpr");
  static_assert(atft::greater_than_equal(3)(4), "greater_equal is constexpr");
  static_assert(atft::negate(eq<3>)(2), "negate is constexpr");
  static_assert(atft::when_all(eq<3>, ne<4>)(3), "when_all is constexpr");
  static_assert(atft::when_any(eq<3>, eq<4>)(4), "when_any is constexpr");
  static_assert(atft::when_none(eq<3>, eq<4>)(5), "when_none is constexpr");
  static_assert(atft::if_then_else(gt<3>, eq<4>, eq<5>)(4),
		"if_then_else is constexpr");
  static_assert(atft::compose(gt<2>, std::plus<>{})(1, 2),
		"compose is constexpr");

  auto constexpr sumgt2 = atft::compose(gt<2>, std::negate<>{}, std::plus<>{});
  static_assert(std::is_invocable_r_v<bool, decltype(sumgt2), int, int>);
  static_assert(std::is_same<bool, decltype(func(sumgt2(1, 2)))>{});
}

TEST_CASE("compose") {
  auto to_string = [](auto t) { return std::to_string(t); };
  WHEN("called with several unary functions") {
    auto twice = [](int i) { return i + i; };
    auto add_one = [](int i) { return i + 1; };
    auto string_plus_one_twice = atft::compose(to_string, twice, add_one);
    THEN("they are chained last to first") {
      REQUIRE(string_plus_one_twice(2) == "6");
    }
  }
  AND_WHEN("called with a unary and a binary function") {
    auto add = [](int x, int y) { return x + y; };
    auto string_add = atft::compose(to_string, add);
    THEN(
	"the composition is binary, calling the last with two values and the "
	"first with the result") {
      REQUIRE(string_add(3, 2) == "5");
    }
  }
  AND_WHEN("called with a binary and a unary function") {
    auto pick_first = [](const auto& x) { return x.first; };
    auto cmp = atft::compose(std::less<>{}, pick_first);
    THEN(
	"the composition is binary, calling the last with each value, and the "
	"first with the two results") {
      REQUIRE(cmp(std::pair(1, 3), std::pair(3, 1)));
    }
  }
  AND_WHEN("functions are non-copyable") {
    auto f1 = [x = std::make_unique<int>(3)](int p) { return *x + p; };
    auto f2 = [y = std::make_unique<std::string>("foo")](auto p) {
      return *y + std::to_string(p);
    };
    THEN("they are moved") {
      auto func = atft::compose(std::move(f2), std::move(f1));
      REQUIRE(func(5) == "foo8");
    }
  }
}
TEST_CASE("not_equal") {
  int* p = new int{3};
  REQUIRE(atft::not_equal(std::unique_ptr<int>(p))(std::unique_ptr<int>()));
  REQUIRE_FALSE(atft::not_equal(std::unique_ptr<int>())(nullptr));
}
TEST_CASE("less_than") {
  REQUIRE_FALSE(atft::less_than(3)(3));
  REQUIRE(atft::less_than(3)(2));
  REQUIRE_FALSE(atft::less_than(3)(4));
}

TEST_CASE("less_equal") {
  REQUIRE(atft::less_than_equal(3)(3));
  REQUIRE(atft::less_than_equal(3)(2));
  REQUIRE_FALSE(atft::less_than_equal(3)(4));
}

TEST_CASE("greater_than") {
  REQUIRE_FALSE(atft::greater_than(3)(3));
  REQUIRE_FALSE(atft::greater_than(3)(2));
  REQUIRE(atft::greater_than(3)(4));
}

TEST_CASE("greater_equal") {
  REQUIRE(atft::greater_than_equal(3)(3));
  REQUIRE_FALSE(atft::greater_than_equal(3)(2));
  REQUIRE(atft::greater_than_equal(3)(4));
}

TEST_CASE("when_all") {
  WHEN("all predicates are true") {
    THEN("they are called in order") {
      int num = 0;

      REQUIRE(atft::when_all(
	  [&](int i) {
	    num += i;
	    return num == 1;
	  },
	  [&](int i) {
	    num += i;
	    return num == 2;
	  },
	  [&](int i) {
	    num += i;
	    return num == 3;
	  })(1));
      REQUIRE(num == 3);
    }
  }
  AND_WHEN("some predicates are false") {
    THEN("if first is false it's the only one tested") {
      int num = 0;
      REQUIRE_FALSE(atft::when_all([&](int i) { return ++num == i; },
				   [&](int i) { return ++num == i; },
				   [&](int i) { return ++num == i; })(0));
      REQUIRE(num == 1);
    }
    AND_THEN("if all but last are true, all are tested") {
      int num = 0;
      REQUIRE_FALSE(atft::when_all([&](int i) { return ++num < i; },
				   [&](int i) { return ++num < i; },
				   [&](int i) { return ++num < i; })(3));
      REQUIRE(num == 3);
    }
    AND_THEN("if first is the only true one, the two first are tested") {
      int num = 0;
      REQUIRE_FALSE(atft::when_all([&](int i) { return ++num == i; },
				   [&](int i) { return ++num == i; },
				   [&](int i) { return ++num == i; })(1));
      REQUIRE(num == 2);
    }
  }
  AND_WHEN("a predicate is not copyable") {
    auto pred = atft::when_all(
	[x = std::make_unique<int>(3)](auto& p) { return p == x; });
    THEN("it is moved") { REQUIRE_FALSE(pred(nullptr)); }
  }
}

TEST_CASE("when_any") {
  WHEN("all predicates are false") {
    THEN("they are called in order") {
      int num = 0;

      REQUIRE_FALSE(atft::when_any(
	  [&](int i) {
	    num += i;
	    return num == 0;
	  },
	  [&](int i) {
	    num += i;
	    return num == 0;
	  },
	  [&](int i) {
	    num += i;
	    return num == 0;
	  })(1));
      REQUIRE(num == 3);
    }
  }
  AND_WHEN("some predicates are true") {
    THEN("if first is true it's the only one tested") {
      int num = 0;
      REQUIRE(atft::when_any([&](int i) { return ++num == i; },
			     [&](int i) { return ++num == i; },
			     [&](int i) { return ++num == i; })(1));
      REQUIRE(num == 1);
    }
    AND_THEN("if all but last are false, all are tested") {
      int num = 0;
      REQUIRE(atft::when_any([&](int i) { return ++num == i; },
			     [&](int i) { return ++num == i; },
			     [&](int i) { return ++num == i; })(3));
      REQUIRE(num == 3);
    }
    AND_THEN("if first is the only false one, the two first are tested") {
      int num = 0;
      REQUIRE(atft::when_any([&](int i) { return ++num == i; },
			     [&](int i) { return ++num == i; },
			     [&](int i) { return ++num == i; })(2));
      REQUIRE(num == 2);
    }
  }
  AND_WHEN("a predicate is not copyable") {
    auto pred = atft::when_any(
	[x = std::make_unique<int>(3)](auto& p) { return p == x; });
    THEN("it is moved") { REQUIRE_FALSE(pred(nullptr)); }
  }
}

TEST_CASE("when_none") {
  WHEN("all predicates are false") {
    THEN("they are called in order") {
      int num = 0;

      REQUIRE(atft::when_none(
	  [&](int i) {
	    num++;
	    return num == i;
	  },
	  [&](int i) {
	    num++;
	    return num == i;
	  },
	  [&](int i) {
	    num++;
	    return num == i;
	  })(0));
      REQUIRE(num == 3);
    }
  }
  AND_WHEN("some predicates are true") {
    THEN("if first is true it's the only one tested") {
      int num = 0;
      REQUIRE_FALSE(atft::when_none([&](int i) { return ++num == i; },
				    [&](int i) { return ++num == i; },
				    [&](int i) { return ++num == i; })(1));
      REQUIRE(num == 1);
    }
    AND_THEN("if all but last are false, all are tested") {
      int num = 0;
      REQUIRE_FALSE(atft::when_none([&](int i) { return ++num == i; },
				    [&](int i) { return ++num == i; },
				    [&](int i) { return ++num == i; })(3));
      REQUIRE(num == 3);
    }
    AND_THEN("if first is the only false one, the two first are tested") {
      int num = 0;
      REQUIRE_FALSE(atft::when_none([&](int i) { return ++num > i; },
				    [&](int i) { return ++num > i; },
				    [&](int i) { return ++num > i; })(1));
      REQUIRE(num == 2);
    }
  }
  AND_WHEN("a predicate is not copyable") {
    THEN("it is moved") {
      auto pred = atft::when_none(
	  [x = std::make_unique<int>(3)](auto& p) { return p == x; });
      REQUIRE(pred(nullptr));
    }
  }
}

TEST_CASE("if_then") {
  WHEN("predicate is true") {
    THEN("action is called with value") {
      int num = 0;
      atft::if_then(atft::equal(3), [&](int i) { num = i; })(3);
      REQUIRE(num == 3);
    }
    AND_THEN("multi parameter action is called with all values") {
      int num = 0;
      atft::if_then(
	  atft::compose(atft::equal(3), [](int x, int y) { return x + y; }),
	  [&](int x, int y) { num = x - y; })(4, -1);
      REQUIRE(num == 5);
    }
    AND_THEN("action can mutate its capture") {
      std::unique_ptr<int> m;
      atft::if_then(atft::not_equal(3),
		    [x = std::make_unique<int>(0), &m](auto n) mutable {
		      *x = n;
		      m = std::move(x);
		    })(4);
      REQUIRE(m);
      REQUIRE(*m == 4);
    }
  }
  AND_WHEN("predicate is false") {
    THEN("action is not called") {
      int num = 0;
      atft::if_then(atft::equal(3), [&](int i) { num = i; })(4);
      REQUIRE(num == 0);
    }
  }
}

TEST_CASE("if_then_else") {
  int tnum = 0;
  int fnum = 0;
  auto condition_if_3 = atft::if_then_else(
      atft::equal(3), [&](int i) { tnum = i; }, [&](int i) { fnum = i; });
  WHEN("predicate is true") {
    THEN("true_action is called with value") {
      condition_if_3(3);
      REQUIRE(tnum == 3);
      REQUIRE(fnum == 0);
    }
  }
  AND_WHEN("predicate is false") {
    THEN("action is not called") {
      condition_if_3(4);
      REQUIRE(tnum == 0);
      REQUIRE(fnum == 4);
    }
  }
  AND_WHEN("actions have return values") {
    THEN("the returned type is the common type") {
      auto op = atft::if_then_else(
	  atft::less_than(0), [](auto) { return 0U; },
	  [](auto x) { return x; });
      auto p1 = op(-1);
      static_assert(
	  std::is_same<decltype(p1), std::common_type_t<int, unsigned>>{});
      REQUIRE(p1 == 0U);
      auto p2 = op(1);
      static_assert(
	  std::is_same<decltype(p2), std::common_type_t<int, unsigned>>{});
      REQUIRE(p2 == 1U);
    }
  }
}

TEST_CASE("do_all") {
  WHEN("there are several functions") {
    THEN("they are called in sequence") {
      int num = 0;
      atft::do_all(
	  [&](int i) {
	    num += i;
	    REQUIRE(num == 1);
	  },
	  [&](int i) {
	    num += i;
	    REQUIRE(num == 2);
	  },
	  [&](int i) {
	    num += i;
	    REQUIRE(num == 3);
	  })(1);
      REQUIRE(num == 3);
    }
  }
  AND_WHEN("functions are non copyable") {
    THEN("they are moved") {
      int n = 0;
      atft::do_all([x = std::make_unique<int>(3), &n](auto p) { n = p + *x; })(
	  5);
      REQUIRE(n == 8);
    }
  }
}
