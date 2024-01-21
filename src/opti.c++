#include "opti.h"

OptimizerParam::OptimizerParam(void) {
  gradtol = 1e-3;
  trialstep = 0.1;
  maxstep = 5.0;
  multline = 2.0;
  maxbrakline = 5;
  maxitline = 5;
  maxitcg = 50;
  badstrike = 1;
  itreset = 0;
}

OptimizerParam::OptimizerParam(const OptimizerParam &a) {
  gradtol = a.gradtol;
  trialstep = a.trialstep;
  maxstep = a.maxstep;
  multline = a.multline;
  maxbrakline = a.maxbrakline;
  maxitline = a.maxitline;
  maxitcg = a.maxitcg;
  badstrike = a.badstrike;
  itreset = a.itreset;
}

void OptimizerParam::operator=(const OptimizerParam &a) {
  gradtol = a.gradtol;
  trialstep = a.trialstep;
  maxstep = a.maxstep;
  multline = a.multline;
  maxbrakline = a.maxbrakline;
  maxitline = a.maxitline;
  maxitcg = a.maxitcg;
  badstrike = a.badstrike;
  itreset = a.itreset;
}

FunctionAlongLine::FunctionAlongLine(const Array<Real> &_org,
				     const Array<Real> &_dir,
				     FunctionWithGrad *_pfunc)
    : Function1DWithGrad(), org(_org), dir(_dir) {
  pfunc = _pfunc;
  val = pfunc->get_val();
  Array<Real> vgrad;
  pfunc->get_grad(&vgrad);
  grad = inner_product(vgrad, dir);
}

void FunctionAlongLine::scalar_to_vector(Array<Real> *px, Real u) {
  product(px, dir, u);
  sum(px, *px, org);
}

void FunctionAlongLine::set_arg(Real u) {
  Array<Real> x;
  scalar_to_vector(&x, u);
  pfunc->set_arg(x);
  val = pfunc->get_val();
  Array<Real> vgrad;
  pfunc->get_grad(&vgrad);
  grad = inner_product(vgrad, dir);
}

FunctionAlongGreatCircle::FunctionAlongGreatCircle(const Array<Real> &_org,
						   const Array<Real> &_dir,
						   FunctionWithGrad *_pfunc)
    : Function1DWithGrad(), org(_org), dir(_dir) {
  pfunc = _pfunc;
  val = pfunc->get_val();
  Array<Real> vgrad;
  pfunc->get_grad(&vgrad);
  grad = inner_product(vgrad, dir);
  cdir = norm(org) / norm(dir);
}

void FunctionAlongGreatCircle::scalar_to_vector(Array<Real> *px, Real u) {
  Array<Real> y;
  product(px, org, cos(u));
  product(&y, dir, sin(u) * cdir);
  sum(px, *px, y);
}

void FunctionAlongGreatCircle::set_arg(Real u) {
  Array<Real> x;
  scalar_to_vector(&x, u);
  pfunc->set_arg(x);
  val = pfunc->get_val();
  Array<Real> vgrad;
  pfunc->get_grad(&vgrad);
  Array<Real> dx;
  scalar_to_vector(&dx, u + M_PI / 2.);
  product(&dx, dx, 1. / cdir);
  grad = inner_product(vgrad, dx);
}

int min_func1d(Real *px, const OptimizerParam &param,
	       Function1DWithGrad *pfunc) {
  Real a, b;
  b = *px;
  a = b;
  Real dfb = pfunc->get_grad();	 // at *px;
  if (fabs(dfb) < param.gradtol) {
    cout << "1d min found immediately" << endl;
    return 1;
  }
  Real sdf0 = sgn(dfb);
  Real dfa = 0.;
  Real curstep = param.trialstep;
  int it;
  for (it = 0; it < param.maxbrakline; it++) {
    cout << "trialbracket\n";
    cout << "a= " << a << " dfa= " << dfa << endl;
    cout << "b= " << b << " dfb= " << dfb << endl;
    Real dx;
    if (fabs(dfb) < fabs(dfa)) {
      dx = min(1.5 * fabs(dfb) * fabs(a - b) / (fabs(dfa) - fabs(dfb)),
	       param.trialstep * param.maxstep);
      cout << "smart_bracket_guess\n";
    } else {
      dx = curstep;
      curstep *= param.multline;
      curstep = min(curstep, param.trialstep * param.maxstep);
    }
    a = b;
    dfa = dfb;
    b += -dx * sdf0;
    pfunc->set_arg(b);
    dfb = pfunc->get_grad();
    if (fabs(dfb) < param.gradtol) {
      cout << "found line min at " << b << " dfx= " << dfb << endl;
      *px = b;
      return 1;
    }
    if (dfb * sdf0 < -param.gradtol) {
      cout << "found bracket:\n";
      cout << "a= " << a << " dfa= " << dfa << endl;
      cout << "b= " << b << " dfb= " << dfb << endl;
      break;
    }
  }
  if (it == param.maxbrakline) {
    *px = b;
    cout << "Too many iterations in bracketing." << endl;
    return 0;
  }
  int badstrike = param.badstrike;
  Real dfx;
  for (it = 0; it < param.maxitline; it++) {
    Real z = dfb / (dfb - dfa);
    Real x = a * z + b * (1 - z);
    if (z < 0.25 || z > 0.75) {
      if (badstrike > 0) {
	badstrike--;
	cout << "badstrike forgiven\n";
      } else {
	cout << "midpoint" << endl;
	x = (a + b) / 2.;
      }
    }
    pfunc->set_arg(x);
    dfx = pfunc->get_grad();
    cout << "a= " << a << " dfa= " << dfa << endl;
    cout << "x= " << x << " dfx= " << dfx << endl;
    cout << "b= " << b << " dfb= " << dfb << endl;
    *px = x;
    if (fabs(dfx) < param.gradtol) {
      cout << "found line min at " << x << " dfx= " << dfx << endl;
      return 1;
    }
    if (dfx * dfa < 0) {
      b = x;
      dfb = dfx;
    } else {
      a = x;
      dfa = dfx;
    }
  }
  if (fabs(dfa) < fabs(dfx)) {
    *px = a;
    dfx = dfa;
  }
  if (fabs(dfb) < fabs(dfx)) {
    *px = b;
    dfx = dfb;
  }
  cout << "Too many iterations in bisection, taking x= " << *px
       << " dfx= " << dfx << endl;
  return 0;
}

void line_search(Array<Real> *px, const Array<Real> &dir,
		 const OptimizerParam &param, FunctionWithGrad *pfunc) {
  Array<Real> ndir(dir);
  normalize(&ndir);
  Array<Real> org(*px);
  FunctionAlongLine func1d(org, ndir, pfunc);
  Real u = 0.;
  min_func1d(&u, param, &func1d);
  func1d.scalar_to_vector(px, u);
}

void great_circle_search(Array<Real> *px, const Array<Real> &dir,
			 const OptimizerParam &param, FunctionWithGrad *pfunc) {
  Array<Real> ndir(dir);
  normalize(&ndir);
  Array<Real> org(*px);
  FunctionAlongGreatCircle func1d(org, ndir, pfunc);
  Real u = 0.;
  min_func1d(&u, param, &func1d);
  func1d.scalar_to_vector(px, u);
}

int conjugate_gradient(Array<Real> *px, const OptimizerParam &_param,
		       FunctionWithGrad *pfunc) {
  OptimizerParam param(_param);
  if (param.itreset == 0) {
    param.itreset = px->get_size();
  }
  pfunc->set_arg(*px);
  Array<Real> dir, olddir;
  pfunc->get_grad(&dir);
  product(&dir, dir, -1.);
  olddir = dir;
  int itreset = 0;
  for (int it = 0; it < param.maxitcg; it++) {
    Array<Real> dx(*px);
    // cout << "LS\n";
    // cout << "dir= "; tracev(dir);
    param.gradtol /= sqrt((Real)px->get_size());
    line_search(px, dir, param, pfunc);
    param.gradtol = _param.gradtol;
    cout << "CG_step " << it << endl;
    diff(&dx, dx, *px);
    Real ndx = norm(dx);
    if (ndx != 0) {
      param.trialstep = min(ndx, _param.trialstep);
    }
    //    if (fabs(pfunc->get_val()-curf)<tol) {return 1;}
    //    curf=pfunc->get_val();
    pfunc->get_grad(&dir);
    product(&dir, dir, -1.);
    // cout << "l_dir= "; tracev(dir);
    cout << "l_gnorm= " << norm(dir) << endl;
    if (norm(dir) < param.gradtol) {
      return 1;
    }
    Array<Real> diffdir, betaolddir;
    diff(&diffdir, dir, olddir);
    Real beta = inner_product(dir, diffdir) / inner_product(olddir, olddir);
    // Real beta=-inner_product(dir,diffdir)/inner_product(olddir,olddir);
    itreset++;
    if (beta < 0 || itreset > param.itreset) {
      cout << "CG Reset: beta=0" << endl;
      beta = 0.;
      itreset = 0;
    }
    product(&betaolddir, olddir, beta);
    sum(&dir, dir, betaolddir);
    olddir = dir;
  }
  cout << "Too many iteration in conjugate gradient" << endl;
  return 0;
}

void ortho_project(Array<Real> *px, const Array<Real> &p,
		   const Array<Real> &x) {
  Real a = inner_product(p, x) / inner_product(p, p);
  Array<Real> ap;
  product(&ap, p, a);
  diff(px, x, ap);  //*px=x-p <p|x>/<p|p)>;
}

int conjugate_gradient_on_sphere(Array<Real> *px, const OptimizerParam &_param,
				 FunctionWithGrad *pfunc) {
  cout << "begin on_sphere" << endl;
  OptimizerParam param(_param);
  if (param.itreset == 0) {
    param.itreset = px->get_size();
  }
  pfunc->set_arg(*px);
  //  Real curf=pfunc->get_val();
  Array<Real> dir, olddir;
  pfunc->get_grad(&dir);
  product(&dir, dir, -1.);
  ortho_project(&dir, *px, dir);
  olddir = dir;
  int itreset = 0;
  for (int it = 0; it < param.maxitcg; it++) {
    Array<Real> dx(*px);
    param.gradtol /= sqrt((Real)px->get_size());
    great_circle_search(px, dir, param, pfunc);
    param.gradtol = _param.gradtol;
    diff(&dx, dx, *px);
    Real ndx = norm(dx);
    if (ndx != 0) {
      param.trialstep = min(ndx / norm(*px), _param.trialstep);
    }
    //    if (fabs(pfunc->get_val()-curf)<tol) {return 1;}
    //    curf=pfunc->get_val();
    pfunc->get_grad(&dir);
    product(&dir, dir, -1.);
    ortho_project(&dir, *px, dir);
    // cout << "s_dir= "; tracev(dir);
    cout << "s_gnorm= " << norm(dir) << endl;
    if (norm(dir) < param.gradtol) {
      cout << "end on_sphere" << endl;
      return 1;
    }
    Array<Real> diffdir, betaolddir;
    diff(&diffdir, dir, olddir);
    Real beta = inner_product(dir, diffdir) / inner_product(olddir, olddir);
    itreset++;
    if (beta < 0 || itreset > param.itreset || norm(dx) > norm(*px)) {
      cout << "sCG Reset: beta=0" << endl;
      beta = 0.;
      itreset = 0;
    }
    product(&betaolddir, olddir, beta);
    sum(&dir, dir, betaolddir);
    ortho_project(&dir, *px, dir);
    olddir = dir;

    cout << "CG_on_sphere step " << it << endl;
  }
  cout << "Too many iteration in conjugate gradient on sphere" << endl;
  return 0;
}

int steepest_descent(Array<Real> *px, const OptimizerParam &param,
		     FunctionWithGrad *pfunc) {
  for (int it = 0; it < param.maxitcg; it++) {
    cout << "Steepest_descent " << it << endl;
    Array<Real> g;
    pfunc->set_arg(*px);
    pfunc->get_grad(&g);
    if (norm(g) < param.gradtol) {
      return 1;
    }
    product(&g, g, -param.trialstep);
    sum(px, *px, g);
  }
  return 0;
}
