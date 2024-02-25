#include <iostream>

#include "range.h"

int main() {
  // print [0, 10)
  std::cout << "range(10): { ";
  for (auto i : itertools::range(10)) {
    std::cout << i << ' ';
  }
  std::cout << "}\n";
}
