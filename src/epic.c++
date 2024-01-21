#include <fstream>

#include "linalg.h"
#include "opti.h"
#include "parse.h"
#include "stringo.h"
#include "strinterf.h"
#include "version.h"
#include "xtalutil.h"

class SaddleEpicycle : public FunctionWithGrad {
  FunctionWithGrad *pfunc;
  OptimizerParam sparam;
  OptimizerParam sparam_save;
  Real slen;
  Array<Real> s0;
  Real energy;
  Array<Real> gradenergy;

 public:
  SaddleEpicycle(FunctionWithGrad *_pfunc)
      : FunctionWithGrad(), sparam(), s0(), gradenergy() {
    pfunc = _pfunc;
  }
  void init(const OptimizerParam &_sparam, Real _slen) {
    sparam = _sparam;
    sparam_save = _sparam;
    slen = _slen;
  }
  void set_arg(const Array<Real> &c) {
    cout << "Ion Motion" << endl;
    SmallEpicycle sepi(pfunc);
    sepi.init(c);
    energy = pfunc->get_val();	// assumes set_arg was called in sepi;
    cout << "energy= " << energy << endl;
    Array<Real> f0, f1;
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
      sparam.trialstep = max(min(1.5 * norm(ds0), sparam_save.trialstep),
			     sparam.trialstep / 10.);
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
    cout << "epipos= " << c.get_size() << " ";
    tracev(c);

    Real b = inner_product(s0, f0) / inner_product(s0, s0);
    Array<Real> fpar, fpar2, fperp;
    product(&fpar, s0, b);
    product(&fpar2, s0, 2. * b);
    diff(&gradenergy, f0, fpar2);

    diff(&fperp, f0, fpar);
    pfunc->get_grad(&f1);  // assume set_arg was called in sepi;
    Real curv =
	(inner_product(f1, s0) - inner_product(f0, s0)) / inner_product(s0, s0);

    cout << "grad_par= ";
    tracev(fpar);
    cout << "grad_perp= ";
    tracev(fperp);
    cout << "grad_tot= ";
    tracev(gradenergy);

    cout << "mincurv= " << curv << "  energy= " << energy
	 << "  grad_norm total= " << norm(gradenergy) << " par= " << norm(fpar)
	 << " perp= " << norm(fperp) << endl;
  }
  Real get_val(void) const { return energy; }
  void get_grad(Array<Real> *df) const { *df = gradenergy; }
};

Real find_saddle(Array<Real> *pc, const OptimizerParam &sparam, Real slen,
		 const OptimizerParam &steepparam,
		 const OptimizerParam &ionparam, FunctionWithGrad *pfunc) {
  SaddleEpicycle mepi(pfunc);
  mepi.init(sparam, slen);
  if (steepest_descent(pc, steepparam, &mepi)) {
    return mepi.get_val();
  }
  cout << "STARTCG" << endl;
  conjugate_gradient(pc, ionparam, &mepi);
  return mepi.get_val();
}

// these are dummy functions for debugging;

Real myfr(Real r) { return -1. / (1. + ipow(r, 2)); }

Real mylen(Real x, Real y) { return sqrt(sqr(x) + sqr(y)); }

Real myf(Real x, Real y) { return myfr(mylen(x, y)) + myfr(mylen(x - 3, y)); }

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
  }
  Real get_val(void) const { return myf(x(0), x(1)); }
  void get_grad(Array<Real> *pdf) const { *pdf = df; }
};

// end dummy functions;

extern char *helpstring;

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
  Real forcefact = 3.;
  int dimstrain = 0;
  int noscale = 0;
  int sleeptime = 5;

  int outputstr = 0;
  int relax = 0;
  int debug = 0;
  int analex = 0;
  int defaults = 0;
  Array<Real> tmpparam(2);
  tmpparam(0) = 0.3;
  tmpparam(1) = 0.2;

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

      {"-ff", "Force scale Factor", REALVAL, &forcefact},
      {"-ds",
       "Dimension of the strain relaxation allowed (0: none (default), 1: "
       "along x, 2: along y,z 3: along x,y,z) ",
       INTVAL, &dimstrain},
      {"-ns", "Do Not Scale tolerance with number of atoms", BOOLVAL, &noscale},
      {"-st", "Sleep Time between read access", INTVAL, &sleeptime},
      {"-d", "Use all default values", BOOLVAL, &defaults},
      {"-os", "Output Structure defined by epipos.out into str.out", BOOLVAL,
       &outputstr},
      {"-r", "Relax to minimum (instead of finding saddle point)", BOOLVAL,
       &relax},
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
      Real dx = 0.01;
      ofstream file("ftmp.out");
      Array<Real> x(2);
      for (x(0) = -2.; x(0) <= 6.; x(0) += 0.05) {
	for (x(1) = -2.; x(1) <= 2.; x(1) += 0.05) {
	  file << x(0) << " " << x(1) << " " << myf(x(0), x(1)) << endl;
	}
	file << endl;
      }
    }

    MyFunc f;
    Array<Real> x(2);
    x = tmpparam;

    if (relax) {
      Real e = conjugate_gradient(&x, ionparam, &f);
    } else {
      Real e = find_saddle(&x, sparam, slen, steepparam, ionparam, &f);
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
    e = conjugate_gradient(&x, ionparam, &optstr);
  } else {
    ifstream file("epipos.out");
    if (file) {
      file >> x;
    }
    e = find_saddle(&x, sparam, slen, steepparam, ionparam, &optstr);
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
  cout << "epic terminated normally" << endl;
}
