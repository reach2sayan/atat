#include <complex.h>

#include <fstream>

#include "getvalue.h"
#include "parse.h"
#include "version.h"

Real dPhi(Real p, Real k, Real kT) {
  Real acc = 0.;
  Real k0val = 1. / (1. + exp(-p / kT));
  Real dk0val = k * kT * log(1 + exp(p / kT));
  if (dk0val / k0val < 0.1) {
    acc = k0val - dk0val;
  } else {
    kT *= M_PI;
    for (Real s = -1; s < 1.5; s += 2.) {
      Complex r = sqrt(k) / sqrt(Complex(p, -s * kT));
      Complex ftfermi(1., 0.);
      if (abs(kT * r) > 0) {
	ftfermi = kT * r / sinh(kT * r);
      }
      acc += real(Complex(0, s) * sqrt(Complex(0, -s)) / (2. * sqrt(M_PI)) *
		  pow(k * Complex(p, -s * kT), -0.25) *
		  exp(-Complex(kT, 2 * s * p) * r) * ftfermi);
    }
  }
  return acc;
}

Complex dPhi_complex(Real p, Real k, Real kT) {
  Complex acc = 0.;
  Real k0val = 1. / (1. + exp(-p / kT));
  Real dk0val = k * kT * log(1 + exp(p / kT));
  if (dk0val / k0val < 0.1) {
    acc = k0val - dk0val;
  } else {
    kT *= M_PI;
    Complex r = sqrt(k) / sqrt(Complex(p, -kT));
    Complex ftfermi(1., 0.);
    if (abs(kT * r) > 0) {
      ftfermi = kT * r / sinh(kT * r);
    }
    acc += Complex(0, 1.) * sqrt(Complex(0, -1.)) / (2. * sqrt(M_PI)) *
	   pow(k * Complex(p, -kT), -0.25) * exp(-Complex(kT, 2 * p) * r) *
	   ftfermi;
  }
  return acc;
}

Real dPhi_pos(Real p, Real k, Real kT) {
  Real acc = 0.;
  Real k0val = 1. / (1 + exp(-p / kT));
  Real dk0val = k * kT * log(1 + exp(p / kT));
  if (dk0val / k0val < 0.1) {
    acc = k0val - dk0val;
  } else {
    kT *= M_PI;
    for (Real s = -1; s < 1.5; s += 2.) {
      Complex r = sqrt(k) / sqrt(Complex(p, -s * kT));
      Real ftfermi = 1.;
      // if (real(kT*r)>0.) {ftfermi=abs(kT*r)/sinh(real(kT*r));}
      if (real(kT * r) > 0.) {
	ftfermi = abs(kT * r / sinh(kT * r));
      }
      acc += (1. / (2. * sqrt(M_PI))) *
	     pow(abs(k * Complex(p, -s * kT)), -0.25) *
	     exp(real(-Complex(kT, 2 * s * p) * r)) * ftfermi;
    }
  }
  return acc;
}

Real dPhi_raw(Real p, Real k, Real kT, Real prange, int pmesh) {
  Real tiny = 1e-4;
  Real dp = 2. * prange / (Real)(2 * pmesh + 1);
  int plmesh = pmesh;
  Real minp = p - (Real)plmesh * dp;
  if (minp < 0) {
    plmesh = (int)floor(p / dp);
    minp = p - (Real)plmesh * dp;
  }
  int p2mesh = plmesh + pmesh + 1;
  Real sk = 2. * sqrt(k);
  Real acc = 0;
  int i;
  Real u;
  for (i = p2mesh, u = minp + dp / 2.; i > 0; i--, u += dp) {
    acc += j0(sk * sqrt(u)) * pow(cosh((u - p) / kT / 2.), -2.);
  }
  acc *= dp / kT / 4.;
  return acc;
}

Real dPhi(Real p, Real k, Real kT, Real prec) {
  Real prange = kT;
  int pmesh = 1;
  int done = 0;
  while (!done) {
    done = 1;
    if (fabs(dPhi_raw(p, k, kT, prange, 2 * pmesh) -
	     dPhi_raw(p, k, kT, prange, pmesh)) > prec) {
      pmesh *= 2;
      done = 0;
    }
    if (fabs(dPhi_raw(p, k, kT, 2. * prange, 2 * pmesh) -
	     dPhi_raw(p, k, kT, prange, pmesh)) > prec) {
      pmesh *= 2;
      prange *= 2.;
      done = 0;
    }
    //    cerr << dPhi_raw(p,k,kT,prange,pmesh) << " " << prange << " " << pmesh
    //    << endl;
  }
  return dPhi_raw(p, k, kT, prange, pmesh);
}

/*
int binomial(int n, int k) {
  return factorial(n)/(factorial(n-k)*factorial(k));
}

Real dPhi(Real p, Real k, int dpmdk, Real a, int nn) {
  Real acc=0.;
  for (int n=1; n<=nn; n++) {
    Real acc2=0.;
    for (int m=1; m<=n; m++) {
      Real rm=(Real)m;
      acc2+=(Real)binomial(n,m)*ipow(rm,dpmdk)/ipow(-a,m)*exp(-k/rm-rm*p);
    }
    Real rn=(Real)n;
    acc+=acc2/(rn*pow(1.+1./a,rn));
  }
  return ipow(-1.,dpmdk)*acc;
}
*/

template <class R, class T>
class Function {
 public:
  virtual R operator()(const T &x) const {}
};

class ActualPot : public Function<Real, rVector3d> {
  Real actualstiffo2;

 public:
  ActualPot(Real _actualstiff)
      : actualstiffo2(_actualstiff / 2.), Function<Real, rVector3d>() {}
  Real operator()(const rVector3d &x) const { return actualstiffo2 * norm2(x); }
};

Real path_pot(Array<rVector3d> xs, const Function<Real, rVector3d> &pot) {
  int maxn = xs.get_size();
  Real dt = 1. / (Real)maxn;
  Real acc = 0.;
  for (int n = 0; n < maxn; n++) {
    acc += pot(xs(n));
  }
  return acc * dt;
}

Real path_mean_pot(Array<rVector3d> xs, const Function<Real, rVector3d> &pot) {
  int maxn = xs.get_size();
  Real dt = 1. / (Real)maxn;
  rVector3d acc(0., 0., 0.);
  for (int n = 0; n < maxn; n++) {
    acc += xs(n);
  }
  return pot(acc * dt);
}

Real path_kin(Array<rVector3d> xs) {
  int maxn = xs.get_size();
  Real dt = 1. / (Real)maxn;
  Real acc = norm2(xs(0) - xs(maxn - 1));
  for (int n = 1; n < maxn; n++) {
    acc += norm2(xs(n) - xs(n - 1));
  }
  return acc / dt;
}

class SamplingFunc : public Function<Real, rVector2d> {
  Real kT;
  Real wscale;

 public:
  SamplingFunc(Real _kT, Real _wscale)
      : kT(_kT), wscale(_wscale), Function<Real, rVector2d>() {}
  Real operator()(const rVector2d &x) const { return dPhi_pos(x(0), x(1), kT); }
};

class MCState {
 public:
  Real pot;
  Real meanpot;
  Real kin;
  Real weight;
  Array<rVector3d> xs;
  MCState(void) : pot(0.), meanpot(0.), kin(0.), weight(0.), xs(0) {}
  MCState(const MCState &old)
      : pot(old.pot),
	meanpot(old.meanpot),
	kin(old.kin),
	weight(old.weight),
	xs(old.xs) {}
};

Complex sgn(Complex z) { return z / abs(z); }

Real determinant(const Array2d<Real> &_a) {
  Array<int> tmp;
  Array2d<Real> a(_a);
  lu_decomposition(&a, &tmp);
  Real acc = 1;
  for (int i = 0; i < a.get_size()(0); i++) {
    acc *= a(i, i);
  }
  return acc;
}

Complex determinant_2x2(const Array2d<Complex> &a) {
  return a(0, 0) * a(1, 1) - a(0, 1) * a(1, 0);
}

void inverse_2x2(Array2d<Complex> *pinv_a, const Array2d<Complex> &a) {
  pinv_a->resize(2, 2);
  Complex det = determinant_2x2(a);
  Complex tmp = a(0, 0) / det;
  (*pinv_a)(0, 0) = a(1, 1) / det;
  (*pinv_a)(1, 1) = tmp;
  (*pinv_a)(0, 1) = -a(0, 1) / det;
  (*pinv_a)(1, 0) = -a(1, 0) / det;
}

Complex c_quadratic_form(const Array2d<Complex> &a, const Array<Complex> &x) {
#ifdef DEBUG
  if (a.get_size()(0) != a.get_size()(1)) ERRORQUITDUMP("Matrix is not square");
  if (a.get_size()(0) != x.get_size())
    ERRORQUITDUMP("Matrices not conformable");
#endif
  Complex result(0., 0.);
  for (int i = 0; i < x.get_size(); i++) {
    for (int j = 0; j < x.get_size(); j++) {
      result += (a(i, j) * conj(x(i)) * x(j));
    }
  }
  return result;
}

void product(Array2d<Complex> *result, const Array2d<Complex> &a,
	     const Array2d<Complex> &b) {
#ifdef DEBUG
  if (a.get_size()(1) != b.get_size()(0))
    ERRORQUITDUMP("Matrices not conformable");
#endif
  result->resize(iVector2d(a.get_size()(0), b.get_size()(1)));
  zero_array(result);
  for (int i = 0; i < a.get_size()(0); i++) {
    for (int j = 0; j < b.get_size()(1); j++) {
      for (int k = 0; k < a.get_size()(1); k++) {
	(*result)(i, j) += a(i, k) * b(k, j);
      }
    }
  }
}

int main(int argc, char *argv[]) {
  Real dx = 0.1;
  Real kT = 1.0;
  Real shift = 0.5;
  Real mu0 = 1.0;
  Real mu1 = 1.0;
  Real dmu = 0.1;
  Real actualstiff = 1.0;
  Real wscale = 1.0;
  Real strstiff = 1.0;
  Real lowcut = exp(8.);
  int nlink = 12;
  int nstepequil = 50;
  int nstepmax = 100;
  int verbstep = 10;
  int printstr = 0;
  int dummy = 0;
  Real maxk = 4.0;
  Real maxp = 4.0;
  Real minp = 0.0;
  int kstep = 20;
  int pstep = 20;
  Real prec = 1e-4;
  Real demean = 0.;
  int usemeanpot = 0;
  int dogrid = 0;
  Real dp = 1e-4;
  Real dk = 1e-4;
  int readmv = 0;
  AskStruct options[] = {
      {"",
       "Path integral electronic structure code " MAPS_VERSION
       " (EXPERIMENTAL), by Axel van de Walle",
       TITLEVAL, NULL},
      CMDLINEPARAM(dx, REALVAL),
      CMDLINEPARAM(kT, REALVAL),
      CMDLINEPARAM(mu0, REALVAL),
      CMDLINEPARAM(mu1, REALVAL),
      CMDLINEPARAM(dmu, REALVAL),
      CMDLINEPARAM(actualstiff, REALVAL),
      CMDLINEPARAM(wscale, REALVAL),
      CMDLINEPARAM(strstiff, REALVAL),
      CMDLINEPARAM(nlink, INTVAL),
      CMDLINEPARAM(nstepmax, INTVAL),
      CMDLINEPARAM(nstepequil, INTVAL),
      CMDLINEPARAM(verbstep, INTVAL),
      CMDLINEPARAM(maxk, REALVAL),
      CMDLINEPARAM(maxp, REALVAL),
      CMDLINEPARAM(minp, REALVAL),
      CMDLINEPARAM(kstep, INTVAL),
      CMDLINEPARAM(pstep, INTVAL),
      CMDLINEPARAM(prec, REALVAL),
      CMDLINEPARAM(dogrid, BOOLVAL),
      CMDLINEPARAM(usemeanpot, BOOLVAL),
      CMDLINEPARAM(printstr, BOOLVAL),
      CMDLINEPARAM(demean, REALVAL),
      CMDLINEPARAM(dp, REALVAL),
      CMDLINEPARAM(dk, REALVAL),
      CMDLINEPARAM(readmv, BOOLVAL),
      {"-d", "Use all default values", BOOLVAL, &dummy}};
  if (!get_values(argc, argv, countof(options), options)) {
    display_help(countof(options), options);
    return 1;
  }
  if (dogrid) {
    for (Real k = 0.; k <= maxk + zero_tolerance; k += maxk / (Real)kstep) {
      for (Real p = minp; p <= maxp + zero_tolerance;
	   p += (maxp - minp) / (Real)pstep) {
	cout << p << " " << k << " " << dPhi(p, k, kT, prec) << " "
	     << dPhi(p, k, kT) << " " << dPhi_pos(p, k, kT) << " "
	     << j0(2 * sqrt(k * p)) << endl;
      }
      cout << endl;
    }
    exit(1);
  }

  ActualPot actualpot(actualstiff);
  SamplingFunc samplingfunc(kT, wscale);
  MCState curstate;
  MCState newstate;
  ifstream dumpin;
  ofstream dumpout;

  if (readmv) {
    dumpin.open("dump.in");
  } else {
    dumpout.open("dump.out");
  }
  curstate.xs.resize(nlink);
  fill_array(&(curstate.xs), rVector3d(0., 0., 0.));
  newstate = curstate;
  for (Real mu = mu0; mu <= mu1 + dmu / 2.; mu += dmu) {
    curstate.weight = 1e-30;
    Real acc = 0.;
    Real nstepacc = 0;
    int nstep = 0;
    Array<Real> mean_pk(2);
    Array2d<Real> var_pk(iVector2d(2, 2));
    if (readmv) {
      Real tmp;
      dumpin >> tmp;
      dumpin >> mean_pk(0) >> mean_pk(1);
      dumpin >> var_pk(0, 0) >> var_pk(0, 1) >> var_pk(1, 0) >> var_pk(1, 1);
    } else {
      zero_array(&mean_pk);
      zero_array(&var_pk);
      while (1) {
	for (int i = 0; i < curstate.xs.get_size(); i++) {
	  newstate.xs(i) = curstate.xs(i) + rVector3d(dx * (uniform01() - 0.5),
						      dx * (uniform01() - 0.5),
						      dx * (uniform01() - 0.5));
	}
	newstate.pot = path_pot(newstate.xs, actualpot) - mu;
	newstate.meanpot = path_mean_pot(newstate.xs, actualpot) - mu;
	newstate.kin = strstiff * path_kin(newstate.xs);
	newstate.weight = samplingfunc(rVector2d(-newstate.pot, newstate.kin));
	Real r = newstate.weight / curstate.weight;
	if (uniform01() <= r) {
	  curstate = newstate;
	}
	Real cur_dPhi = dPhi(-curstate.pot, curstate.kin, kT);
	Real cur_dPhi_mean = dPhi(-curstate.meanpot, curstate.kin, kT);
	Real curweightmean =
	    samplingfunc(rVector2d(-curstate.meanpot, curstate.kin));
	Real cur_dPhiow = cur_dPhi / curstate.weight;
	nstep++;
	if (nstep > nstepequil) {
	  acc += cur_dPhiow;
	  nstepacc += 1.;
	  mean_pk(0) += curstate.pot;
	  mean_pk(1) += curstate.kin;
	  var_pk(0, 0) += curstate.pot * curstate.pot;
	  var_pk(0, 1) += curstate.pot * curstate.kin;
	  var_pk(1, 0) += curstate.pot * curstate.kin;
	  var_pk(1, 1) += curstate.kin * curstate.kin;
	  if (nstep % verbstep == 0) {
	    if (printstr) {
	      for (int i = 0; i < curstate.xs.get_size(); i++) {
		cout << curstate.xs(i) << endl;
	      }
	      cout << curstate.xs(0) << endl;
	      cout << endl;
	    } else {
	      cout << mu << " " << nstep << " " << curstate.pot << " "
		   << curstate.meanpot << " " << curstate.kin << " " << cur_dPhi
		   << " " << curstate.weight << " " << curweightmean << " "
		   << cur_dPhiow << " " << acc / nstepacc << endl;
	    }
	  }
	}
	if (nstep >= nstepmax) break;
      }
      product(&mean_pk, mean_pk, 1. / nstepacc);
      product(&var_pk, var_pk, 1. / nstepacc);
      Array2d<Real> tmp;
      outer_product(&tmp, mean_pk, mean_pk);
      diff(&var_pk, var_pk, tmp);
      dumpout << mu << " ";
      dumpout << mean_pk(0) << " " << mean_pk(1) << " ";
      dumpout << var_pk(0, 0) << " " << var_pk(0, 1) << " " << var_pk(1, 0)
	      << " " << var_pk(1, 1) << " " << endl;
    }
    Array<Real> freq(2);
    Real p = mean_pk(0);
    Real k = mean_pk(1);
    Real sdk = dk * sqrt(mean_pk(1));
    freq(0) =
	arg(dPhi_complex(-(p + dp), k, kT) / dPhi_complex(-(p - dp), k, kT)) /
	(2. * dp);
    freq(1) =
	arg(dPhi_complex(-p, k + sdk, kT) / dPhi_complex(-p, k - sdk, kT)) /
	(2. * sdk);
    Array2d<Real> quadph(2, 2);
    quadph(0, 0) =
	arg(dPhi_complex(-(p + dp), k, kT) * dPhi_complex(-(p - dp), k, kT) /
	    sqr(dPhi_complex(-(p), k, kT))) /
	sqr(dp);
    quadph(1, 1) =
	arg(dPhi_complex(-(p), k + sdk, kT) * dPhi_complex(-(p), k - sdk, kT) /
	    sqr(dPhi_complex(-(p), k, kT))) /
	sqr(sdk);
    quadph(0, 1) = arg((dPhi_complex(-(p + dp), k + sdk, kT) /
			dPhi_complex(-(p + dp), k - sdk, kT)) /
		       (dPhi_complex(-(p - dp), k + sdk, kT) /
			dPhi_complex(-(p - dp), k - sdk, kT))) /
		   (4. * dp * sdk);
    quadph(1, 0) = quadph(0, 1);
    Real ph = arg(sgn(dPhi_complex(-p, k, kT)));
    cout << "ph(x,y)=(" << ph << ")+(" << freq(0) << ")*x+(" << freq(1)
	 << ")*y+(" << quadph(0, 0) / 2. << ")*x*x+(" << quadph(1, 1) / 2.
	 << ")*y*y+2.*(" << quadph(0, 1) / 2. << ")*x*y" << endl;
    cout << "phs(x,y)=ph(x-(" << p << "),y-(" << k << "))" << endl;
    Array2d<Complex> W(2, 2);
    for (MultiDimIterator<iVector2d> ii(iVector2d(2, 2)); ii; ii++) {
      W(ii) = var_pk(ii);
    }
    cout << W << endl;
    inverse_2x2(&W, W);
    cout << W << endl;
    for (MultiDimIterator<iVector2d> ii(iVector2d(2, 2)); ii; ii++) {
      W(ii) += Complex(0., quadph(ii));
    }
    cout << W << endl;
    Array2d<Complex> tmp1, tmp2;
    tmp1 = W;
    inverse_2x2(&W, W);
    product(&tmp2, W, tmp1);
    cout << "I=" << tmp2 << endl;
    cout << W << endl;
    Array<Complex> cfreq(2);
    for (int i = 0; i < 2; i++) {
      cfreq(i) = freq(i);
    }

    cout << cfreq << endl;
    //    cout << "ne= " << mu << " " <<
    //    real(2.*M_PI*exp(Complex(0.,ph))*exp(-0.5*c_quadratic_form(W,cfreq))*sqrt(determinant_2x2(W)))
    //    << endl;
    cout << "ne= " << mu << " "
	 << real(2. * M_PI * exp(-0.5 * c_quadratic_form(W, cfreq)) *
		 sqrt(determinant_2x2(W)))
	 << endl;
  }
}
