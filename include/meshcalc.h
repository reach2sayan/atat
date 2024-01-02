#ifndef _MESHCALC_H_
#define _MESHCALC_H_

#include <iostream>
#include "version.h"
#include "integer.h"
#include "linalg.h"
#include "tensor.h"

void find_perp_to_vector_array(Array<Real> *pperp, const Array<Array<Real> > & v);

int calc_circumsphere(Array<Real> *pcenter, Real *pradius2, const Array<Array<Real> > &pts);

class SideFlag {
public:
  int done;
  Array<int> ipts;
  Array<Real> forward;
  Array<Real> normal;
  SideFlag(void): done(0),ipts(0),forward(0),normal(0) {}
  SideFlag(int _done, const Array<int>& _ipts): done(_done),ipts(_ipts),forward(0),normal(0) {}
};

void debug_mesh(const LinkedList<Array<int> > & polylist, const Array<Array<Real> > &pts);
void debug_side(const LinkedList<SideFlag > & edgelist, const Array<Array<Real> > &pts);

class RealIntPair {
public:
  Real r;
  int i;
  RealIntPair(void): r(0.), i(0) {}
  RealIntPair(Real _r, int _i): r(_r), i(_i) {}
};

template<class T>
class OrderBy_r {
 public:
  int operator () (const T& a, const T& b) const {
    return (a.r<b.r);
  }
};

class PointGridIndex {
public:
  Array<Real> low;
  Real boxsize;
  Tensor<LinkedList<int> > indices;
public:
  PointGridIndex(const Array<Array<Real> > &pts);
};

class PointFromCenterIteratorGridIndex {
  const Array<Array<Real> > &pts;
  Array<Real> center;
  Array<int> box_center;
  Array<RealIntPair> distindex;
  int cur_index;
  const PointGridIndex &grid_index;
  Real cur_rad;
  int cur_boxrad;
  
  int fill_buffer(int new_boxrad);
  void update_if_needed(void);
public:
  PointFromCenterIteratorGridIndex(const Array<Array<Real> > &_pts, const PointGridIndex &_grid_index, const Array<Real> &_center);
  void init(void);
  operator void * ();
  int operator () (void);
  void operator++(int);
};

class PointFromCenterIterator {
  const Array<Array<Real> > &pts;
  Array<Real> center;
  Array<RealIntPair> distindex;
  int curindex;
public:
  PointFromCenterIterator(const Array<Array<Real> > &_pts, const Array<Real> &_center);
  void init(void);
  operator void * ();
  int operator () (void);
  void operator++(int);
};

void calc_normal_vector(Array<Real> *pnormal, const Array<int> &sideipts, int ptipts, const Array<Array<Real> > &pts);
void create_mesh(LinkedList<Array<int> > *ppolylist, const Array<Array<Real> > &pts, int fulldim, Real maxr);
int remove_nearby(Array<Array<Real> > *ppts, const Array<Array<Real> > &orgpts, Real cutoff);
void cross_section_simplexes(Array<Array<Real> > *pcutpts,  LinkedList<Array<int> > *pcutpoly, const Array<Array<Real> > &axes, const Array<Real> &org, const Array<Array<Real> > &pts, const LinkedList<Array<int> > &poly);
void cross_section_simplexes(Array<Array<Real> > *pcutpts,  LinkedList<Array<int> > *pcutpoly, Array<Array<Real> > *pcurcolors, const Array<Array<Real> > &axes, const Array<Real> &org, const Array<Array<Real> > &pts, const LinkedList<Array<int> > &poly, const  Array<Array<Real> > &colors);
void simple_cross_section_simplexes(Array<Array<Real> > *pcutpts,  LinkedList<Array<int> > *pcutpoly, Array<Array<Real> > *pcutcolors, const Array<int> &iaxes, const Array<Real> &org, const Array<Array<Real> > &pts, const LinkedList<Array<int> > &poly, const  Array<Array<Real> > &colors);

#endif
