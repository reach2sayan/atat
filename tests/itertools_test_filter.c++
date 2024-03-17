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
