#include <fstream>

#include "getvalue.h"
#include "linalg.h"
#include "parse.h"
#include "version.h"

Real pot(const Array<Real> &x) {
  return pow(abs(1. + cos(x(0)) * cos(1.2 * x(1)) + 0.4 * x(1)), 1.7);
}

void pot_der(Array<Real> *pf, const Array<Real> &x) {
  Real dx = 1e-4;
  pf->resize(x.get_size());
  zero_array(pf);
  for (int i = 0; i < x.get_size(); i++) {
    Array<Real> y;
    Array<Real> dy(x.get_size());
    zero_array(&dy);
    dy(i) = dx;
    sum(&y, x, dy);
    (*pf)(i) = pot(y);
    diff(&y, x, dy);
    (*pf)(i) -= pot(y);
    (*pf)(i) /= (2. * dx);
  }
}

int main(int argc, char *argv[]) {
  int dim = 2;
  int n = 20;
  Real lfact = 50.;
  Real dlfact = 4e-1;
  Real dfact = 1e-2;
  Real totall = M_PI * 0.5;
  Real told = 1e-4;

  //  for (totall=M_PI/8.; totall<M_PI*2.0; totall+=M_PI/32.) {
  Array<Array<Real> > band(n + 1);
  for (int i = 0; i <= n; i++) {
    band(i).resize(dim);
  }
  band(0)(0) = M_PI / 4;
  band(0)(1) = M_PI / 10.;
  band(n)(0) = M_PI;
  band(n)(1) = 0.0;
  for (int i = 1; i < n; i++) {
    Real r = (Real)i / (Real)(n);
    Array<Real> acc1, acc2;
    product(&acc1, band(0), 1. - r);
    product(&acc2, band(n), r);
    sum(&band(i), acc1, acc2);
  }
  Real sumfm;
  while (1) {
    sumfm = 0.;
    Real onel = totall / (Real)(n + 1);
    Array<Array<Real> > dband(n + 1);
    zero_array(&dband);
    Array<Real> l(n);
    Array<Array<Real> > dx(n);
    for (int i = 0; i < n; i++) {
      diff(&dx(i), band(i + 1), band(i));
      l(i) = norm(dx(i));
    }
    Real dtotall = 0.;
    for (int i = 1; i < n; i++) {
      Array<Real> fl, fl2;
      product(&fl, dx(i - 1), (onel - l(i - 1)) / l(i - 1));
      product(&fl2, dx(i), (l(i) - onel) / l(i));
      sum(&fl, fl, fl2);
      product(&fl, fl, lfact);
      Array<Real> u(dim);
      diff(&u, band(i + 1), band(i - 1));
      normalize(&u);
      Array<Real> f, fu;
      pot_der(&f, band(i));
      product(&fu, u, inner_product(u, f));
      diff(&f, f, fu);
      dtotall =
	  (inner_product(f, dx(i - 1)) - inner_product(f, dx(i))) * dlfact;
      sumfm += norm(f);
      diff(&dband(i), fl, f);
      product(&dband(i), dband(i), dfact);
    }
    Real tol = 0.;
    for (int i = 1; i < n; i++) {
      sum(&band(i), band(i), dband(i));
      tol += norm(dband(i));
    }
    totall += dtotall;
    cerr << sumfm << " " << totall << endl;
    for (int i = 0; i < n + 1; i++) {
      for (int j = 0; j < dim; j++) {
	cout << band(i)(j) << " ";
      }
      cout << pot(band(i)) << endl;
    }
    cout << endl;
    if (tol < told) {
      break;
    }
  }
}

