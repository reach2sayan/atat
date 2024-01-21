#include <fstream>

#include "getvalue.h"
#include "linalg.h"
#include "parse.h"
#include "version.h"

template <class T>
void project(Array<T> *pv, const Array2d<T> &proj, const Array<T> &v,
	     int paral) {
  Array2d<T> A, iA;
  Array<T> b, c;
  if (proj.get_size()(0) == 0) {
    if (paral) {
      pv->resize(v.get_size());
      zero_array(pv);
    } else {
      *pv = v;
    }
  }
  inner_product(&A, proj, proj);
  invert_matrix(&iA, A);
  inner_product(&b, proj, v);
  product(&c, iA, b);
  if (paral) {
    product(pv, proj, c);
  } else {
    Array<Real> dv;
    product(&dv, proj, c);
    diff(pv, v, dv);
  }
}

class FuncofVector {
 public:
  virtual Real operator()(const Array<Real> &) const { return 0; };
};

class PolyForm {
 public:
  Array<Real> center;
  Real coef0;
  Array<rTensor> coef;
  Real operator()(const Array<Real> &x) const {
    Array<Real> dx;
    diff(&dx, x, center);
    Real a = coef0;
    for (int n = 0; n < coef.get_size(); n++) {
      Array<int> size = coef(n).get_size();
      MultiDimIterator<Array<int> > i(size);
      for (; i; i++) {
	Real p = coef(n)(i);
	for (int j = 0; j < size.get_size(); j++) {
	  p *= dx(((Array<int> &)i)(j));
	}
	a += p;
      }
    }
    return a;
  }
};

class InterpolPolyForm : public FuncofVector {
  Array<PolyForm> polyform;
  Array<Array<int> > connect;
  Array<Array<int> > border;
  Array<Array<Array<Real> > > plane_n;
  Array<Array<Real> > plane_c;
  Array<Array<Real> > maxdist;

 public:
  InterpolPolyForm(const Array<PolyForm> &_polyform,
		   const Array<Array<int> > &_connect) {
    polyform = _polyform;
    connect = _connect;
    border.resize(connect.get_size());
    Array<Array<Array<Real> > > plane_o;
    plane_o.resize(connect.get_size());
    plane_n.resize(connect.get_size());
    for (int s = 0; s < connect.get_size(); s++) {
      int nv = connect(s).get_size();
      border(s).resize(nv);
      plane_o(s).resize(nv);
      plane_n(s).resize(nv);
      for (int x = 0; x < nv; x++) {
	border(s)(x) = 1;
	plane_o(s)(x) = polyform(connect(s)((x + 1) % nv)).center;
	Array<Real> u;
	diff(&u, polyform(connect(s)(x)).center,
	     polyform(connect(s)((x + 1) % nv)).center);
	Array2d<Real> proj(iVector2d(u.get_size(), nv - 2));
	for (int j = 0; j < proj.get_size()(1); j++) {
	  Array<Real> spanv;
	  diff(&spanv, polyform(connect(s)((x + 2 + j) % nv)).center,
	       polyform(connect(s)((x + 1) % nv)).center);
	  set_column(&proj, spanv, j);
	}
	project(&plane_n(s)(x), proj, u, 0);
	normalize(&plane_n(s)(x));
      }
    }

    iVector2d s(-1, -1);
    for (s(0) = 0; s(0) < connect.get_size(); s(0)++) {
      iVector2d x(-1, -1);
      for (x(0) = 0; x(0) < connect(s(0)).get_size(); x(0)++) {
	for (s(1) = 0; s(1) < s(0); s(1)++) {
	  int allin = 1;
	  for (int v = 0; v < connect(s(0)).get_size(); v++) {
	    if (v != x(0)) {
	      if (!is_in_array(connect(s(1)), connect(s(0))(v))) {
		allin = 0;
	      }
	    }
	  }
	  if (allin) {
	    for (x(1) = 0; x(1) < connect(s(1)).get_size(); x(1)++) {
	      if (!is_in_array(connect(s(0)), x(1))) break;
	    }
	    border(s(0))(x(0)) = 0;
	    border(s(1))(x(1)) = 0;
	    Array<Real> n;
	    diff(&n, plane_n(s(0))(x(0)), plane_n(s(1))(x(1)));
	    normalize(&n);
	    plane_n(s(0))(x(0)) = n;
	    product(&n, n, -1.);
	    plane_n(s(1))(x(1)) = n;
	    break;
	  }
	}
      }
    }
    plane_c.resize(connect.get_size());
    maxdist.resize(connect.get_size());
    for (int s = 0; s < connect.get_size(); s++) {
      plane_c(s).resize(connect(s).get_size());
      maxdist(s).resize(connect(s).get_size());
      for (int x = 0; x < connect(s).get_size(); x++) {
	plane_c(s)(x) = inner_product(plane_n(s)(x), plane_o(s)(x));
	maxdist(s)(x) =
	    inner_product(plane_n(s)(x), polyform(connect(s)(x)).center) -
	    plane_c(s)(x);
      }
    }
  }
  Real smooth_func(Real t) const { return pow(t, 3.); }
  Real operator()(const Array<Real> &x) const {
    for (int s = 0; s < plane_n.get_size(); s++) {
      Array<Real> p(plane_n(s).get_size());
      int i;
      for (i = 0; i < p.get_size(); i++) {
	p(i) = inner_product(plane_n(s)(i), x) - plane_c(s)(i);
	if (!border(s)(i) && p(i) < 0) {
	  break;
	}
      }
      if (i == p.get_size()) {
	Real num = 0.;
	Real den = 0.;
	for (int j = 0; j < p.get_size(); j++) {
	  if (p(j) >= 0) {
	    Real w = smooth_func(p(j) / maxdist(s)(j));
	    num += w * polyform(connect(s)(j))(x);
	    den += w;
	  }
	}
	return num / den;
      }
    }
  }
};

int read_vector(Array<Real> *pa, istream &file) {
  LinkedList<Real> l;
  while (1) {
    Real r;
    file >> r;
    if (file.eof()) break;
    l << new Real(r);
  }
  LinkedList_to_Array(pa, l);
}

int read_vector(Array<Real> *pa, const char *filename) {
  ifstream file(filename);
  return read_vector(pa, file);
}

int read_square_matrix(Array2d<Real> *pa, istream &file) {
  LinkedList<Real> l;
  while (1) {
    Real r;
    file >> r;
    if (file.eof()) break;
    l << new Real(r);
  }
  int n = (int)sqrt((Real)l.get_size());
  iVector2d size(n, n);
  pa->resize(size);
  MultiDimIterator<iVector2d> j(size);
  LinkedListIterator<Real> i(l);
  for (; i; i++, j++) {
    (*pa)(j) = *i;
  }
}

int read_square_matrix(Array2d<Real> *pa, const char *filename) {
  ifstream file(filename);
  return read_square_matrix(pa, file);
}

int read_equi_tensor(Tensor<Real> *pa, int dim, istream &file) {
  LinkedList<Real> l;
  while (1) {
    Real r;
    file >> r;
    if (file.eof()) break;
    l << new Real(r);
  }
  int n = (int)pow((Real)l.get_size(), 1. / (Real)dim);
  Array<int> size(dim);
  for (int d = 0; d < dim; d++) {
    size(d) = n;
  }
  pa->resize(size);
  MultiDimIterator<Array<int> > j(size);
  LinkedListIterator<Real> i(l);
  for (; i; i++, j++) {
    (*pa)(j) = *i;
  }
}

int read_equi_tensor(Tensor<Real> *pa, int dim, const char *filename) {
  ifstream file(filename);
  return read_equi_tensor(pa, dim, file);
}

char *helpstring = "Insert more help here";

int main(int argc, char *argv[]) {
  char *hesfilename = "hessian.out";
  int n = 0;
  int sigdig = 5;
  int dohelp = 0;
  int dummy = 0;
  AskStruct options[] = {
      {"", "Skeleton for atat-like codes" MAPS_VERSION ", by Axel van de Walle",
       TITLEVAL, NULL},
      {"-h", "Input file defining the lattice (defaults: lat.in)", STRINGVAL,
       &hesfilename},
      {"-n", "An integer parameter", INTVAL, &n},
      {"-sig", "Number of significant digits printed (default: 5)", INTVAL,
       &sigdig},
      {"-h", "Display more help", BOOLVAL, &dohelp},
      {"-d", "Use all default values", BOOLVAL, &dummy}};
  if (!get_values(argc, argv, countof(options), options)) {
    display_help(countof(options), options);
    return 1;
  }
  if (dohelp) {
    cout << helpstring;
    return 1;
  }
  Array<PolyForm> node(3);
  read_vector(&node(0).center, "pos0.in");
  {
    ifstream file("e0.in");
    file >> node(0).coef0;
  }
  node(0).coef.resize(2);
  read_equi_tensor(&node(0).coef(0), 1, "der0.in");
  read_equi_tensor(&node(0).coef(1), 2, "hes0.in");

  read_vector(&node(1).center, "pos1.in");
  {
    ifstream file("e1.in");
    file >> node(1).coef0;
  }
  node(1).coef.resize(2);
  read_equi_tensor(&node(1).coef(0), 1, "der1.in");
  read_equi_tensor(&node(1).coef(1), 2, "hes1.in");

  read_vector(&node(2).center, "pos2.in");
  {
    ifstream file("e2.in");
    file >> node(2).coef0;
  }
  node(2).coef.resize(2);
  read_equi_tensor(&node(2).coef(0), 1, "der2.in");
  read_equi_tensor(&node(2).coef(1), 2, "hes2.in");
  /*
  read_vector(&node(3).center,"pos3.in");
  {ifstream file("e3.in"); file >> node(3).coef0;}
  node(3).coef.resize(2);
  read_equi_tensor(&node(3).coef(0),1,"der3.in");
  read_equi_tensor(&node(3).coef(1),2,"hes3.in");
  */
  Array<Array<int> > connect;
  {
    ifstream cfile("con.in");
    cfile >> connect;
  }
  InterpolPolyForm pot(node, connect);
  Array<Real> x(3);
  for (x(0) = -1; x(0) <= 1; x(0) += 0.05) {
    for (x(1) = -1; x(1) <= 1; x(1) += 0.05) {
      for (x(2) = 0; x(2) <= 0; x(2) += 0.05) {
	cout << x(0) << " " << x(1) << " " << pot(x) << endl;
      }
    }
    cout << endl;
  }
}

/*
  ifstream hesfile(hesfilename);
  Array2d<Real> hes(iVector2d(n,n));
  MultiDimIterator<iVector2d> i(iVector2d(n,n));
  for (; i; i++) {
    hesfile >> hes(i);
  }
  //  cout << hes << endl;
  Array<Real> lambda;
  Array2d<Real> vect;
  diagonalize_symmetric_matrix(&lambda,&vect,hes);
  for (int j=0; j<lambda.get_size(); j++) {
    cout << lambda(j) << endl;
  }

*/
