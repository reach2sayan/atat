#include <fstream>

#include "linalg.h"
#include "opti.h"
#include "parse.h"
#include "stringo.h"
#include "strinterf.h"
#include "version.h"
#include "xtalutil.h"

class MotionEpicycle : public FunctionWithGrad {
  FunctionWithGrad *pfunc;
  OptimizerParam sparam;
  OptimizerParam sparam_save;
  Real slen;
  Real constif;
  Real nearcutoff;
  Array<Real> s0;
  Array<Real> oldc;
  Real energy;
  Array<Real> gradenergy;

 public:
  MotionEpicycle(FunctionWithGrad *_pfunc)
      : FunctionWithGrad(), sparam(), s0(), oldc(), gradenergy() {
    pfunc = _pfunc;
  }
  void init(const OptimizerParam &_sparam, Real _slen, Real _constif,
	    Real _nearcutoff) {
    sparam = _sparam;
    sparam_save = _sparam;
    slen = _slen;
    constif = _constif;
    nearcutoff = _nearcutoff;
  }
  void set_arg(const Array<Real> &c) {
    cout << "Ion Motion" << endl;
    SmallEpicycle sepi(pfunc);
    sepi.init(c);
    if (oldc.get_size() == 0) {
      oldc = c;
    }
    energy = pfunc->get_val();	// assumes set_arg was called in sepi;
    Array<Real> f0;
    pfunc->get_grad(&f0);

    system("cp str.out str_current.out");

    cout << "grad(c)= ";
    tracev(f0);
    if (s0.get_size() == 0) {
      ifstream file("epidir.out");
      if (file) {
	file >> s0;
      } else {
	product(&s0, f0, -slen / norm(f0));
	//	for (int i=0; i<s0.get_size(); i++) {
	//	  s0(i)+=(2*uniform01()-1.)*sparam.trialstep*slen;
	//	}
      }
      product(&s0, s0, slen / norm(s0));
    }
    Array<Real> ds0(s0);
    conjugate_gradient_on_sphere(&s0, sparam, &sepi);
    diff(&ds0, ds0, s0);
    if (norm(ds0) != 0.0) {
      Array<Real> dc;
      diff(&dc, oldc, c);
      if (norm(dc) <= nearcutoff) {
	sparam.trialstep = max(min(1.5 * norm(ds0), sparam_save.trialstep),
			       sparam.trialstep / 2.);
      } else {
	sparam.trialstep = sparam_save.trialstep;
      }
    }

    cout << "epi_trial_step= " << sparam.trialstep << endl;
    cout << "epidir= " << s0.get_size() << " ";
    tracev(s0);
    {
      ofstream file("epidir.out");
      file << s0;
    }

    {
      ofstream file("epipos.out");
      file << c;
    }
    oldc = c;
    cout << "epipos= " << c.get_size() << " ";
    tracev(c);

    Array<Real> f1, f2;
    Real e1, e2;
    Array<Real> cs;

    sum(&cs, c, s0);
    // pfunc->set_arg(cs);
    pfunc->get_grad(&f1);
    e1 = pfunc->get_val();
    cout << "plotend1= ";
    tracevl(cs);
    cout << e1 << endl;

    diff(&cs, c, s0);
    pfunc->set_arg(cs);
    pfunc->get_grad(&f2);
    e2 = pfunc->get_val();
    cout << "plotend2= ";
    tracevl(cs);
    cout << e2 << endl << endl << endl;

    cout << "grad(c+s)= ";
    tracev(f1);
    cout << "grad(c-s)= ";
    tracev(f2);
    Real s2 = inner_product(s0, s0);
    Real mcurv = -(e1 + e2 - 2. * energy) / s2;
    Array<Real> gradmcurv;
    product(&gradmcurv, f0, -2.);
    sum(&gradmcurv, gradmcurv, f1);
    sum(&gradmcurv, gradmcurv, f2);
    product(&gradmcurv, gradmcurv, -1. / s2);

    {
      Array<Real> tmp;
      cout << "plotend3= ";
      tracevl(c);
      cout << energy << endl;
      product(&tmp, gradmcurv, -slen / norm(gradmcurv));
      sum(&tmp, c, tmp);
      cout << "plotend4= ";
      tracevl(tmp);
      cout << energy << endl << endl << endl;
    }

    Array<Real> fperp;
    ortho_project(&fperp, gradmcurv, f0);

    // Real zerocurv=-mcurv/norm(gradmcurv);
    if (constif == 0.) {
      //      constif=norm(fperp)/max(fabs(zerocurv),1e-2);
      constif = norm(fperp) / fabs(mcurv);
      cout << "constraint_stiffness= " << constif << endl;
    }
    Array<Real> fpar;

    product(&fpar, gradmcurv, mcurv * constif / norm(gradmcurv));

    sum(&gradenergy, fpar, fperp);
    // cout << "constraint_closeness= " << zerocurv << endl;
    cout << "grad_par= ";
    tracev(fpar);
    cout << "grad_perp= ";
    tracev(fperp);
    cout << "grad_tot= ";
    tracev(gradenergy);
    cout << "grad_norm par= " << norm(fpar) << " perp= " << norm(fperp) << endl;
    cout << "mincurv= " << -mcurv << "  energy= " << energy
	 << "  grad_norm= " << norm(gradenergy) << endl;
  }
  Real get_val(void) const { return energy; }
  void get_grad(Array<Real> *df) const { *df = gradenergy; }
};

Real find_min_in_unstable(Array<Real> *pc, const OptimizerParam &sparam,
			  Real slen, const OptimizerParam &steepparam,
			  const OptimizerParam &ionparam, Real constif,
			  FunctionWithGrad *pfunc, int *psuccess = NULL) {
  MotionEpicycle mepi(pfunc);
  mepi.init(sparam, slen, constif, ionparam.trialstep);
  if (steepest_descent(pc, steepparam, &mepi)) {
    return mepi.get_val();
  }
  cout << "STARTCG" << endl;
  int success = conjugate_gradient(pc, ionparam, &mepi);
  if (psuccess) {
    *psuccess = success;
  }
  return mepi.get_val();
}

// these are dummy functions for debugging;

Real myfr(Real r) { return ipow(r, 2) / 2. - ipow(r, 3) / 6.; }

Real myf(Real x, Real y) {
  //  return pow(1+cos(x),1.5)+pow(1+cos(y),1.5)+sqr(z);
  //  return
  //  myfr(sqrt(sqr(x)+sqr(y)))+1.0*sqr(y)+x/8.+3.*exp(-sqr(x-2.3)-sqr(y-0.5))+sqr(x)+0.2*y;
  return myfr(sqrt(sqr(x) + sqr(y))) + 1.0 * sqr(y) + x / 8. +
	 3. * exp(-sqr(x - 2.3) - sqr(0.75 * (y - 0.5))) + sqr(x) + 0.2 * y;
}

class MyFunc : public FunctionWithGrad {
  Array<Real> x;
  Array<Real> df;

 public:
  MyFunc(void) : FunctionWithGrad(), x() {}
  void set_arg(const Array<Real> &_x) {
    x = _x;
    df.resize(x.get_size());
    Real dx = 0.01;
    df(0) = (myf(x(0) + dx, x(1)) - myf(x(0) - dx, x(1))) / (2. * dx);
    df(1) = (myf(x(0), x(1) + dx) - myf(x(0), x(1) - dx)) / (2. * dx);
    // cout << "gnorm= " << norm(df) << endl;
  }
  Real get_val(void) const { return myf(x(0), x(1)); }
  void get_grad(Array<Real> *pdf) const { *pdf = df; }
};

class MyHiDFunc : public FunctionWithGrad {
  Array<Real> lambda;
  Array<Real> slope;
  Real f;
  Array<Real> df;

 public:
  int dim;
  MyHiDFunc(void) : FunctionWithGrad(), f(0), df(), lambda(), slope() {
    dim = 100;
    lambda.resize(dim);
    slope.resize(dim);
    for (int i = 0; i < lambda.get_size(); i++) {
      lambda(i) = uniform01();
      slope(i) = 0 * uniform01();
    }
    cout << "lambda= ";
    tracev(lambda);
    cout << endl;
  }
  void set_arg(const Array<Real> &x) {
    product_diag(&df, lambda, x);
    f = inner_product(df, x) / 2.;
    sum(&df, df, slope);
    for (int i = 0; i < lambda.get_size(); i++) {
      df(i) += 0 * (uniform01() - 0.5);
    }
    f += inner_product(slope, x);
  }
  Real get_val(void) const { return f; }
  void get_grad(Array<Real> *pdf) const { *pdf = df; }
};

class MyFunk : public Function1DWithGrad {
  Real x;
  Real f(Real y) const {
    return -8. / (1. + ipow(y - 1, 8)) - y * exp(-0.5 * y);
  }

 public:
  void set_arg(Real _x) {
    x = _x;
    cout << x << " " << f(x) << endl;
  }
  Real get_val(void) const { return f(x); }
  Real get_grad(void) const {
    Real dx = 1e-3;
    return (f(x + dx) - f(x - dx)) / (2. * dx);
  }
};

// end dummy functions;

extern const char *helpstring;

int main(int argc, char *argv[]) {
  // parse command line arguments or display help (see getvalue.h);
  int dohelp = 0;

  int sigdig = 5;
  OptimizerParam sparam, steepparam, ionparam;
  sparam.gradtol = 0.04;
  sparam.trialstep = 0.1;
  sparam.maxstep = 5.;
  sparam.multline = 2.0;
  sparam.maxbrakline = 5;
  sparam.maxitline = 3;
  sparam.maxitcg = 20;
  sparam.badstrike = 2;
  sparam.itreset = 5;
  Real slen = 0.1;

  steepparam.trialstep = 3e-3;
  steepparam.maxitcg = 0;

  ionparam.gradtol = 0.1;
  ionparam.trialstep = 0.005;
  ionparam.maxstep = 5.;
  ionparam.multline = 2.0;
  ionparam.maxbrakline = 5;
  ionparam.maxitline = 3;
  ionparam.maxitcg = 40;
  ionparam.badstrike = 2;
  ionparam.itreset = 5;
  Real constif = 30.0;
  Real forcefact = 3.;
  int dimstrain = 3;
  int noscale = 0;
  int sleeptime = 5;

  int outputstr = 0;
  int relax = 0;
  int debug = 0;
  int analex = 0;
  int defaults = 0;
  Array<Real> tmpparam(2);
  tmpparam(0) = 1.9;
  tmpparam(1) = -2.0;

  AskStruct options[] = {
      {"", "INFlection DETection method " MAPS_VERSION ", by Axel van de Walle",
       TITLEVAL, NULL},
      {"-h", "Display more help", BOOLVAL, &dohelp},
      {"-egt", "Epicycle Gradient Tolerance", REALVAL, &sparam.gradtol},
      {"-ets", "Epicycle Trial Step", REALVAL, &sparam.trialstep},
      {"-ems", "Epicycle Maximum Step multiplier", REALVAL, &sparam.maxstep},
      {"-eml", "Epicycle Multiplier in Line minimization", REALVAL,
       &sparam.multline},
      {"-ebl",
       "Epicycle maximum number of Bracketing steps in Line minimizations",
       INTVAL, &sparam.maxbrakline},
      {"-eil", "Epicycle maximum number of Iterations in Line minimizations",
       INTVAL, &sparam.maxitline},
      {"-eic", "Epicycle maximum number of Iterations in Conjugate gradient",
       INTVAL, &sparam.maxitcg},
      {"-ebf", "Epicycle motion Bad step Forgiveness", INTVAL,
       &sparam.badstrike},
      {"-eir", "Epicycle iterations between conjugate gradient reset", INTVAL,
       &sparam.itreset},
      {"-el", "Epicycle Length", REALVAL, &slen},

      {"-ism", "Ion motion Steepest descent Multiplier", REALVAL,
       &steepparam.trialstep},
      {"-isi", "Ion motion Steepest maximum number of Iterations", INTVAL,
       &steepparam.maxitcg},
      {"-igt", "Ion motion Gradient Tolerance", REALVAL, &ionparam.gradtol},
      {"-its", "Ion motion Trial Step", REALVAL, &ionparam.trialstep},
      {"-ims", "Ion motion Maximum Step multiplier", REALVAL,
       &ionparam.maxstep},
      {"-iml", "Ion motion Multiplier in Line minimization", REALVAL,
       &ionparam.multline},
      {"-ibl",
       "Ion motion maximum number of Bracketing steps in Line minimizations",
       INTVAL, &ionparam.maxbrakline},
      {"-iil", "Ion motion maximum number of Iterations in Line minimizations",
       INTVAL, &ionparam.maxitline},
      {"-iic", "Ion motion maximum number of Iterations in Conjugate gradient",
       INTVAL, &ionparam.maxitcg},
      {"-ibf", "Ion motion Bad step Forgiveness", INTVAL, &ionparam.badstrike},
      {"-iir", "Ion motion iterations between conjugate gradient reset", INTVAL,
       &ionparam.itreset},
      {"-ics", "Ion motion Curvature constraint Stiffness (0=automatic)",
       REALVAL, &constif},

      {"-ff", "Force scale Factor", REALVAL, &forcefact},
      {"-ds",
       "Dimension of the strain relaxation allowed (0: none, 1: along x, 2: "
       "along y,z, 3: along x,y,z (default)) ",
       INTVAL, &dimstrain},
      {"-ns", "Do Not Scale tolerance with number of atoms", BOOLVAL, &noscale},
      {"-st", "Sleep Time between read access", INTVAL, &sleeptime},
      {"-d", "Use all default values", BOOLVAL, &defaults},
      {"-os", "Output Structure defined by epipos.out into str.out", BOOLVAL,
       &outputstr},
      {"-r", "Relax to minimum", BOOLVAL, &relax},
      {"", "For debugging only:", TITLEVAL, NULL},
      {"-db", "debug", BOOLVAL, &debug},
      {"-aex", "Analytical Example (for debugging)", BOOLVAL, &analex},
      {"-t1", "Temp param 1", REALVAL, &tmpparam(0)},
      {"-t2", "Temp param 2", REALVAL, &tmpparam(1)}};
  if (!get_values(argc, argv, countof(options), options)) {
    display_help(countof(options), options);
    return 1;
  }
  if (dohelp) {
    cout << helpstring;
    return 1;
  }

  steepparam.gradtol = ionparam.gradtol;

  if (analex) {
    {
      MyHiDFunc f;
      Array<Real> s0(f.dim);
      for (int i = 0; i < s0.get_size(); i++) {
	s0(i) = uniform01();
      }
      product(&s0, s0, slen / norm(s0));
      //      conjugate_gradient_on_sphere(&s0,sparam,&f);
      conjugate_gradient(&s0, sparam, &f);
      cout << "optimized dir:";
      tracev(s0);
      cout << endl;
      exit(0);
    }

    {
      Real dx = 0.01;
      ofstream file("ftmp.out");
      Array<Real> x(2);
      for (x(0) = -2.; x(0) <= 2.; x(0) += 0.05) {
	for (x(1) = -2.; x(1) <= 2.; x(1) += 0.05) {
	  Real gxx = (myf(x(0) + dx, x(1)) + myf(x(0) - dx, x(1)) -
		      2. * myf(x(0), x(1))) /
		     sqr(dx);
	  Real gyy = (myf(x(0), x(1) + dx) + myf(x(0), x(1) - dx) -
		      2. * myf(x(0), x(1))) /
		     sqr(dx);
	  Real gxy = (myf(x(0) + dx, x(1) + dx) - myf(x(0) + dx, x(1) - dx) -
		      (myf(x(0) - dx, x(1) + dx) - myf(x(0) - dx, x(1) - dx))) /
		     sqr(dx) / 4.0;
	  Real t = (gxx + gyy) / 2.;
	  Real d = gxx * gyy - gxy * gxy;  // 1*x^2-2*t*x+d;
	  Real mincurv = t - sqrt(sqr(t) - d);
	  file << x(0) << " " << x(1) << " " << myf(x(0), x(1)) << " "
	       << mincurv << endl;
	}
	file << endl;
      }
    }

    MyFunc f;
    Array<Real> x(2);
    //  x(0)=1.2; x(1)=0.5;
    //    x(0)=1.9; x(1)=-0.2;
    x = tmpparam;

    if (relax) {
      conjugate_gradient(&x, ionparam, &f);
    } else {
      Real e = find_min_in_unstable(&x, sparam, slen, steepparam, ionparam,
				    constif, &f);
    }
    exit(0);
  }

  rMatrix3d axes;
  Structure str;
  Array<AutoString> atom_label;
  {
    Array<Array<int> > site_type_list;
    ifstream file("str.in");
    if (!file) {
      ERRORQUIT("Unable to open str.in");
    }
    parse_lattice_file(&str.cell, &str.atom_pos, &str.atom_type,
		       &site_type_list, &atom_label, file, &axes);
    fix_atom_type(&str, site_type_list);
  }

  if (!noscale) {
    Real scaletol = (Real)str.atom_pos.get_size();
    sparam.gradtol *= scaletol;
    ionparam.gradtol *= scaletol;
    steepparam.gradtol *= scaletol;
  }

  OptimizedStructure optstr;
  Array<Real> x;
  optstr.init(&x, axes, str, atom_label, forcefact, dimstrain);
  optstr.sleeptime = sleeptime;

  Real e;
  if (outputstr) {
    ifstream file("epipos.out");
    if (file) {
      file >> x;
    }
    Structure curstr;
    optstr.vect_to_str(&curstr, x);
    {
      ofstream file("str.out");
      file.setf(ios::fixed);
      file.precision(5);
      write_structure(curstr, atom_label, axes, file, 0);
    }
  } else if (relax) {
    conjugate_gradient(&x, ionparam, &optstr);
    e = optstr.get_val();
  } else if (debug) {
    for (int j = 0; j < 6; j++) {
      x(j) = (uniform01() - 0.5) * 0.1;
    }
    for (int i = 0; i < ionparam.maxitline; i++) {
      optstr.set_arg(x);
      cout << "energy: " << optstr.get_val() << endl;
      Array<Real> g;
      optstr.get_grad(&g);
      cout << "pred: "
	   << optstr.get_val() + g(ionparam.maxitcg) * ionparam.trialstep
	   << endl;
      x(ionparam.maxitcg) += ionparam.trialstep;
    }
    /*
    MotionEpicycle epi(&optstr);
    epi.init(sparam, slen, constif);
    for (int i=0; i<ionparam.maxitline; i++) {
      epi.set_arg(x);
      cout << "energy: " << epi.get_val() << endl;
      Array<Real> g;
      epi.get_grad(&g);
      cout << "pred: " << epi.get_val()+g(ionparam.maxitcg)*ionparam.trialstep
    << endl; x(ionparam.maxitcg)+=ionparam.trialstep;
    }
    */
  } else {
    ifstream file("epipos.out");
    if (file) {
      file >> x;
    }
    int success = 0;
    e = find_min_in_unstable(&x, sparam, slen, steepparam, ionparam, constif,
			     &optstr, &success);
    if (!success) {
      ofstream file("error");
      file << "Error during conjugate gradient optimization." << endl;
    }
  }

  {
    ofstream ofile("cenergy.out");
    ofile.setf(ios::fixed);
    ofile.precision(sigdig);
    ofile << e << endl;
  }

  Structure relstr;
  optstr.vect_to_str(&relstr, x);
  {
    ofstream ofile("cstr_relax.out");
    ofile.setf(ios::fixed);
    ofile.precision(sigdig);
    write_structure(relstr, atom_label, axes, ofile, 0);
  }
  cout << "infdet terminated normally" << endl;

  /*
  MyFunc f;

  OptimizerParam ionparam;
  ionparam.gradtol=1e-4;
  ionparam.trialstep=0.02;
  ionparam.maxstep=1.0;
  ionparam.multline=2.0;
  ionparam.maxitline=20;
  ionparam.maxitcg=40;

  Array<Real> x(2);
  x(0)=1.8; x(1)=-0.3;

  conjugate_gradient(&x,ionparam,&f);
  */

  /*
    OptimizerParam param;
    param.gradtol=1e-3;
    param.trialstep=0.02;
    param.maxstep=0.3;
    param.multline=2.0;
    param.maxitline=20;
    MyFunk func;
    Real x=0;
    func.set_arg(x);
    min_func1d(&x,param,&func);
    cout << "Min= " << x << endl;
  */

  /*
  {
    rMatrix3d s;
    Array<Real> x(7);
    for (int i=0; i<7; i++) {x(i)=i;}
    for (int d=3; d>=0; d--) {
      Array6_to_Matrix3d(&s,x,d);
      cout << d << endl;
      cout << s << endl;
    }
    for (int i=0; i<3; i++) {
      for (int j=0; j<3; j++) {
	s(i,j)=i*3+j;
      }
    }
    for (int d=3; d>=0; d--) {
      cout << num_strain_dim(d) << endl;
      for (int i=0; i<x.get_size(); i++) {x(i)=1.1;}
      Matrix3d_to_Array6(&x,s,d);
      cout << d << endl;
      cout << x << endl;
    }

    exit(1);
  }
  */
}
