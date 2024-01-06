#ifndef __STRINTERF_H__
#define __STRINTERF_H__

#include "opti.h"
#include "stringo.h"
#include "xtalutil.h"
#include <fstream>

int num_strain_dim(int n);

void Array6_to_Matrix3d(rMatrix3d *pm, const Array<Real> &x, int dimstrain);

void Matrix3d_to_Array6(Array<Real> *px, const rMatrix3d &m, int dimstrain);

class OptimizedStructure : public FunctionWithGrad {
  Real val;
  Array<Real> grad;
  rVector3d sumpos;
  Real stressscale;
  Real strainscale;
  int dimstrain;

public:
  rMatrix3d axes;
  Structure str;
  Array<AutoString> atom_label;
  int sleeptime;

public:
  OptimizedStructure(void);
  ~OptimizedStructure(void);
  void init(Array<Real> *px, const rMatrix3d &_axes, const Structure &_str,
            const Array<AutoString> &_atom_label, Real forcefact = 1.,
            int _dimstrain = 3);
  void vect_to_str(Structure *pstr, const Array<Real> &x);
  void set_arg(const Array<Real> &x);
  Real get_val(void) const { return val; }
  void get_grad(Array<Real> *pgrad) const { *pgrad = grad; }
};

class SmallEpicycle : public FunctionWithGrad {
  FunctionWithGrad *pfunc;
  Array<Real> c;
  Array<Real> s;
  Array<Real> g0;
  Array<Real> df;
  Real e;

public:
  SmallEpicycle(FunctionWithGrad *_pfunc)
      : FunctionWithGrad(), c(), s(), g0(), df() {
    pfunc = _pfunc;
  }
  void init(const Array<Real> _c);
  void get_grad_org(Array<Real> *pg) { *pg = g0; }
  void set_arg(const Array<Real> &_s);
  Real get_val(void) const { return e; }
  void get_grad(Array<Real> *pdf) const { *pdf = df; }
};

#endif
