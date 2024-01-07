#include "thermofunctions.h"

double sroCorrectionFunction(double x, double a1, double b1, double a2) {
  return a1 - a1 * exp(b1 / x) - a2 / x + (a2 / x) * exp(b1 / x);
}

double gaussian(double x, double a, double b, double c) {
  const double z = (x - b) / c;
  return a * std::exp(-0.5 * z * z);
}

Array<double> linspace(double a, double b, size_t n) {
  Array<double> res(n);
  const auto step = (b - a) / (n - 1);
  auto val = a;
  for (int i = 0; i < n; i++) {
    res[i] = val;
    val += step;
  }
  return res;
}
