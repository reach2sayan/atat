#include <fstream>
#include "linalg.h"

Real func_pol(Real r, Real theta) {
  Real A=0.5;
  Real sigma=0.5;
  return (cos(r)-A*exp(-(sqr(r-M_PI)+sqr(theta))/(2.*sqr(sigma))));
}

Array2d<Real> z;

Real func(const Array<Real> &x) {
  iVector2d ix((int)round(x(0)),(int)round(x(1)));
  return z(ix);
  /*
  Real r=norm(x);
  Real theta=atan2(x(1),x(0));
  return func_pol(r,theta);
  */
  //  return ipow(cos(x(0)),3)*cos(x(1));
  //return sqr(x(0))+sqr(x(1))-sqr(x(0))*sqr(x(1));
  //return ((x(1)-sqr(x(0)))-pow(x(1)-sqr(x(0)),3.)/2.);
}

//2D only
Real curv(const Array<Real> &x, int perp) {
  Real dx=1;
  Array<Real> grad(x.get_size()),pos;
  Array2d<Real> hess(x.get_size(),x.get_size());
  for (int i=0; i<x.get_size(); i++) {
    pos=x;
    hess(i,i)=-2.*func(pos);
    pos(i)=x(i)+dx;
    grad(i)=func(pos);
    hess(i,i)+=func(pos);
    pos(i)=x(i)-dx;
    grad(i)-=func(pos);
    hess(i,i)+=func(pos);
  }
  product(&grad,grad,dx/norm(grad));
  if (perp) {
    swap(&grad(0),&grad(1));
    grad(0)*=-1.;
  }
  pos=x; pos(0)+=dx; pos(1)+=dx; hess(0,1) =func(pos);
  pos=x; pos(0)-=dx; pos(1)+=dx; hess(0,1)-=func(pos);
  pos=x; pos(0)+=dx; pos(1)-=dx; hess(0,1)-=func(pos);
  pos=x; pos(0)-=dx; pos(1)-=dx; hess(0,1)+=func(pos);
  hess(1,0)=hess(0,1);
  Real c=quadratic_form(hess,grad);
  return c;
}

/*
Real curv(const Array<Real> &x) {
  Real dx=2;
  //Real dx=1e-1;
  Array<Real> grad(x.get_size()),pos;
  for (int i=0; i<x.get_size(); i++) {
    pos=x;
    pos(i)=x(i)+dx;
    grad(i)=func(pos);
    pos(i)=x(i)-dx;
    grad(i)-=func(pos);
  }
  product(&grad,grad,dx/norm(grad));
  Real c=-2.*func(x);
  sum(&pos,x,grad);
  c+=func(pos);
  diff(&pos,x,grad);
  c+=func(pos);
  return c;
}

Real curvp(const Array<Real> &x) {
  Real dx=1;
  Array<Real> grad(x.get_size()),pos;
  for (int i=0; i<x.get_size(); i++) {
    pos=x;
    pos(i)=x(i)+dx;
    grad(i)=func(pos);
    pos(i)=x(i)-dx;
    grad(i)-=func(pos);
  }
  swap(&grad(0),&grad(1));
  grad(0)*=-1.;
  product(&grad,grad,dx/norm(grad));
  Real c=-2.*func(x);
  sum(&pos,x,grad);
  c+=func(pos);
  diff(&pos,x,grad);
  c+=func(pos);
  return c;
}
*/

Real smoothstep(Real x) {
  return atan(x)/(M_PI/2.);
}

int main(int argc, char *argv[]) {
  iVector2d sz;
  rVector2d xmin,xmax;
  char tmpc;
  cin >> tmpc;
  cin >> sz;
  cin >> xmin;
  cin >> xmax;
  //  cout << tmpc << " " << sz << " " << xmin << " " << xmax << endl;
  z.resize(sz);
  MultiDimIterator<iVector2d> i(sz);
  for (; i; i++) {
    Real tmp;
    cin >> tmp >> tmp >> z(i);
    //    cout << tmp << " " << tmp << " " << z(i) << endl;
  }

  Array<Real> x(2);
  for (x(0)=1; x(0)<sz(0)-1.; x(0)+=1.) {
    for (x(1)=1.; x(1)<sz(1)-1.; x(1)+=1.) {
      cout << xmin(0)+(xmax(0)-xmin(0))*(Real)x(0)/(Real)(sz(0)-1) << " " << xmin(1)+(xmax(1)-xmin(1))*(Real)x(1)/(Real)(sz(1)-1) << " " << func(x) << " " << curv(x,0) << " " << curv(x,1) << endl;
    }
    cout << endl;
  }
}
