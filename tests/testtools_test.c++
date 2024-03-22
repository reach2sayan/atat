#include "testtools.h"

#include <catch2/catch_test_macros.hpp>
#include <utility>

#include "itertools/base.h"

using ATATIteratorTools::is_iterator_v;
using ATATIteratorTools::is_move_constructible_only;
using MonolistInt = ATATTestTools::MonolithObject<int>;

namespace {

class ValidIter {
 private:
  int i;

 public:
  ValidIter& operator++();    // prefix
  ValidIter operator++(int);  // postfix
  bool operator==(const ValidIter&) const;
  bool operator!=(const ValidIter&) const;
  int operator*();
  void* operator->();
};
}  // namespace

TEST_CASE("IsIterator fails when missing prefix ++", "[helpers]") {
  struct InvalidIter : ValidIter {
    InvalidIter& operator++() = delete;
  };

  REQUIRE(!is_iterator_v<InvalidIter>);
}

TEST_CASE("IsIterator fails when missing postfix ++", "[helpers]") {
  struct InvalidIter : ValidIter {
    InvalidIter operator++(int) = delete;
  };

  REQUIRE(!is_iterator_v<InvalidIter>);
}

TEST_CASE("IsIterator fails when missing ==", "[helpers]") {
  struct InvalidIter : ValidIter {
    bool operator==(const InvalidIter&) const = delete;
  };

  REQUIRE(!is_iterator_v<InvalidIter>);
}

TEST_CASE("IsIterator fails when missing !=", "[helpers]") {
  struct InvalidIter : ValidIter {
    bool operator!=(const InvalidIter&) const = delete;
  };

  REQUIRE(!is_iterator_v<InvalidIter>);
}

TEST_CASE("IsIterator fails when missing *", "[helpers]") {
  struct InvalidIter : ValidIter {
    int operator*() = delete;
  };

  REQUIRE(!is_iterator_v<InvalidIter>);
}

TEST_CASE("IsIterator fails when missing copy-ctor", "[helpers]") {
  struct InvalidIter : ValidIter {
    InvalidIter(const InvalidIter&) = delete;
  };

  REQUIRE(!is_iterator_v<InvalidIter>);
}

TEST_CASE("IsIterator fails when missing copy assignment", "[helpers]") {
  struct InvalidIter : ValidIter {
    InvalidIter& operator=(const InvalidIter&) = delete;
  };

  REQUIRE(!is_iterator_v<InvalidIter>);
}

TEST_CASE("IsIterator passes a valid iterator", "[helpers]") {
  REQUIRE(is_iterator_v<ValidIter>);
}

struct HasNothing {
  HasNothing(const HasNothing&) = delete;
  HasNothing& operator=(const HasNothing&) = delete;
};

struct HasMoveAndCopyCtor {
  HasMoveAndCopyCtor(const HasMoveAndCopyCtor&);
  HasMoveAndCopyCtor(HasMoveAndCopyCtor&&);
};

struct HasMoveCtorAndAssign {
  HasMoveCtorAndAssign(HasMoveCtorAndAssign&&);
  HasMoveCtorAndAssign& operator=(HasMoveCtorAndAssign&&);
};

struct HasMoveCtorAndCopyAssign {
  HasMoveCtorAndCopyAssign(HasMoveCtorAndCopyAssign&&);
  HasMoveCtorAndCopyAssign& operator=(const HasMoveCtorAndCopyAssign&);
};

struct HasMoveCtorOnly {
  HasMoveCtorOnly(HasMoveCtorOnly&&);
};

TEST_CASE("IsMoveConstructibleOnly false without move ctor", "[helpers]") {
  REQUIRE_FALSE(is_move_constructible_only<HasNothing>::value);
}

TEST_CASE("IsMoveConstructibleOnly false with copy ctor", "[helpers]") {
  REQUIRE_FALSE(is_move_constructible_only<HasMoveAndCopyCtor>::value);
}

TEST_CASE("IsMoveConstructibleOnly false with move assign", "[helpers]") {
  REQUIRE_FALSE(is_move_constructible_only<HasMoveCtorAndAssign>::value);
}
TEST_CASE("IsMoveConstructibleOnly false with copy assign", "[helpers]") {
  REQUIRE_FALSE(is_move_constructible_only<HasMoveCtorAndCopyAssign>::value);
}

TEST_CASE("IsMoveConstructibleOnly true when met", "[helpers]") {
  REQUIRE(is_move_constructible_only<HasMoveCtorOnly>::value);
}
