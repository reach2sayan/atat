#ifndef __KMECI_H__
#define __KMECI_H__

#include "getvalue.h"
#include "linklist.h"
#include "xtalutil.h"
#include <complex>

class KSpaceECI {
public:
  virtual void init(const Structure &_lattice,
                    const Array<Array<int>> &_site_type_list,
                    const Array<AutoString> &atom_label,
                    const iVector3d &_supercell,
                    const Array<Array<Array<Real>>> &_corrfunc) {}
  virtual void
  get_k_space_eci(Array<Array<Array<Array<Array<Complex>>>>> *p_ft_eci,
                  const Array<Real> &x) {}
};

class MultiKSpaceECI : public KSpaceECI, public LinkedList<KSpaceECI> {
public:
  MultiKSpaceECI(void) : KSpaceECI(), LinkedList<KSpaceECI>() {}
  void get_k_space_eci(Array<Array<Array<Array<Array<Complex>>>>> *p_ft_eci,
                       const Array<Real> &x);
};

#endif
