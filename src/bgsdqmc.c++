#include <fstream>

#include "array.h"
#include "normal.h"

class FuncofVector {
 public:
  virtual Real operator()(const Array<Real> &) const { return 0; };
};

Real posherm01(void) { return (sqrt(-log(uniform01()))); }

int find_in_sorted_array(Real w, const Array<Real> &cumw) {
  int l = 0;
  int h = cumw.get_size() - 1;
  while (l < h) {
    int t = (l + h) / 2;
    if (cumw(t) > w) {
      h = t;
    } else {
      l = t + 1;
    }
  }
  return l;
}

Real iterate_QMC(Array<Array<Real> > *pnew_walker,
		 const Array<Array<Real> > &walker, const FuncofVector &pot,
		 Real spread, Real mult) {
  int n = walker.get_size();
  int d = walker(0).get_size();
  Array<Real> cumw(n);
  Real neww = 0.;
  Real E = 0.;
  //  cerr << n << endl;
  for (int i = 0; i < n; i++) {
    Real V = pot(walker(i));
    neww += max(1. - mult * V, 0.);
    E += V;
    cumw(i) = neww;
    //    cerr << neww << endl;
  }
  E /= (Real)n;
  pnew_walker->resize(n);
  for (int i = 0; i < n; i++) {
    Real w = neww * uniform01();
    (*pnew_walker)(i) = walker(find_in_sorted_array(w, cumw));
  }
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < d; j++) {
      (*pnew_walker)(i)(j) = (*pnew_walker)(i)(j) + spread * normal01();
    }
  }
  return E;
}

class MyPotential : public FuncofVector {
 public:
  Real operator()(const Array<Real> &x) const {
    Real V = ipow(x(0), 2) + ipow(x(1), 2) + ipow(x(2), 2) +
	     0.1 * (ipow(x(0), 4) + ipow(x(1), 4) + ipow(x(2), 4));
    //(2.)*(0.5*ipow(x(0),3)
    //-1.5*x(0)*ipow(x(1),2))+(6./16.)*ipow(ipow(x(0),2)+ipow(x(1),2),2);
    Real r = 0.;
    for (int i = 3; i < x.get_size(); i++) {
      r += ipow(x(i), 2);
    }
    return V + r;
    /*
    return
    (ipow(x(0),2)-4.*ipow(x(1),2)+ipow(x(1),4)+ipow(x(0),2)*ipow(x(2),2));
    return (ipow(x(0),2)-4.*ipow(x(1),2)+ipow(x(1),4));
    Real V=0;
    for (int i=0; i<x.get_size(); i++) {V+=0.5*ipow(x(i),2);}
    return V;
    return (0.5*ipow(x(0),2)+0.5*4*ipow(x(1),2)+0.5*ipow(x(2),2));
    */
  };
};

int main(int argc, char *argv[]) {
  rndseed(0);
  int dim = 24;
  Real hbar = 1.;
  Real mass = 1.;
  Real dt = 0.01;
  Real spread = sqrt(2. * hbar / (2. * mass) * dt);
  Real mult = dt / hbar;
  int nwalker = 5000;
  int tstep = 50;
  Real wid0 = 1.0;
  Array<Array<Real> > *pwalker = new Array<Array<Real> >(nwalker);
  Array<Array<Real> > *pnew_walker = new Array<Array<Real> >(nwalker);
  for (int i = 0; i < pwalker->get_size(); i++) {
    (*pwalker)(i).resize(dim);
    (*pnew_walker)(i).resize(dim);
    for (int j = 0; j < dim; j++) {
      if (j < 3) {
	if (j == 0) {
	  (*pwalker)(i)(j) = wid0 * posherm01();
	} else {
	  (*pwalker)(i)(j) = wid0 * normal01();
	}
      } else {
	(*pwalker)(i)(j) = wid0 * normal01();
      }
    }
  }
  MyPotential mypotential;
  for (int t = 0; t < tstep; t++) {
    Real energy = iterate_QMC(pnew_walker, *pwalker, mypotential, spread, mult);
    cout << energy << endl;
    swap(&pwalker, &pnew_walker);
  }
  {
    ofstream file("qmc.out");
    for (int i = 0; i < nwalker; i++) {
      for (int j = 0; j < dim; j++) {
	file << (*pwalker)(i)(j) << " ";
      }
      file << endl;
    }
  }
}
