#include "linalg.h"
#ifndef __OPTI_H__
#define __OPTI_H__

class OptimizerParam {
public: 
  Real gradtol;
  Real trialstep;
  Real maxstep;
  Real multline;
  int maxbrakline;
  int maxitline;
  int maxitcg;
  int badstrike;
  int itreset;
  OptimizerParam(void);
  OptimizerParam(const  OptimizerParam &a);
  void operator = (const  OptimizerParam &a);
};

class FunctionWithGrad {
public:
  virtual void set_arg(const Array<Real> &x) {}
  virtual Real get_val(void) const {return 0.;}
  virtual void get_grad(Array<Real> *df) const {}
};

class Function1DWithGrad {
public: 
  virtual void set_arg(Real x) {}
  virtual Real get_val(void) const {return 0.;}
  virtual Real get_grad(void) const {return 0.;}
};

class FunctionAlongLine: public  Function1DWithGrad {
  Array<Real> org;
  Array<Real> dir;
  FunctionWithGrad *pfunc;
  Real val;
  Real grad;
public:
  FunctionAlongLine(const Array<Real> &_org, const Array<Real> &_dir,  FunctionWithGrad *_pfunc);
  virtual void scalar_to_vector(Array<Real> *px, Real u);
  virtual void set_arg(Real u);
  virtual Real get_val(void) const {return val;}
  virtual Real get_grad(void) const {return grad;}
};

class FunctionAlongGreatCircle: public  Function1DWithGrad {
  Array<Real> org;
  Array<Real> dir;
  Real cdir;
  FunctionWithGrad *pfunc;
  Real val;
  Real grad;
public:
  FunctionAlongGreatCircle(const Array<Real> &_org, const Array<Real> &_dir,  FunctionWithGrad *_pfunc);
  virtual void scalar_to_vector(Array<Real> *px, Real u);
  virtual void set_arg(Real u);
  virtual Real get_val(void) const {return val;}
  virtual Real get_grad(void) const {return grad;}
};

int min_func1d(Real *px, const OptimizerParam &param, Function1DWithGrad *pfunc);

void line_search(Array<Real> *px, const Array<Real> &dir, const OptimizerParam &param, FunctionWithGrad *pfunc);

void great_circle_search(Array<Real> *px, const Array<Real> &dir, const OptimizerParam &param, FunctionWithGrad *pfunc);

int conjugate_gradient(Array<Real> *px, const OptimizerParam &_param, FunctionWithGrad *pfunc);

void ortho_project(Array<Real> *px, const Array<Real> &p, const Array<Real> &x);

int conjugate_gradient_on_sphere(Array<Real> *px, const OptimizerParam &_param, FunctionWithGrad *pfunc);

int steepest_descent(Array<Real> *px, const OptimizerParam &param, FunctionWithGrad *pfunc);

#endif
