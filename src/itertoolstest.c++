#include <iostream>

#include "enumerate.h"
#include "range.h"

int main() {
  // print [0, 10)
  std::cout << "range(10): { ";
  for (auto i : itertools::range(10)) {
    std::cout << i << ' ';
  }
  std::cout << "}\n";

  std::cout << "enumerating the characters of a string \"hello\":\n";
  const std::string const_string("hello");
  for (auto&& [i, c] : itertools::enumerate(const_string)) {
    std::cout << '(' << i << ", " << c << ") ";
  }
  std::cout << '\n';

  std::vector<int> vec = {20, 30, 50};
  std::cout << "enumerating a vector of {20, 30, 50}:\n";
  for (auto&& [i, n] : itertools::enumerate(vec)) {
    std::cout << '(' << i << ", " << n << ") ";
    n = 0;
  }
}
