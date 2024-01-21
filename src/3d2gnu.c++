#include <fstream>

#include "getvalue.h"
#include "linalg.h"
#include "meshutil.h"
#include "version.h"
//#include "integer.h"

template <class T, class C>
class OrderByLookup {
  const Array<T> &lookup;
  C c;

 public:
  OrderByLookup(const Array<Real> &_lookup) : lookup(_lookup) {}
  int operator()(int i, int j) const { return (c(lookup(i), lookup(j))); }
};

/*
class CompareDistance {
  int operator() (const rVector3d &a,const rVector3d &b) {return (a(2)<b(2));}
}
*/

class IsGreater {
 public:
  int operator()(Real a, Real b) const { return (a > b); }
};

void make_2dview(Array<Array<rVector2d> > *ptri, Array<Array<Real> > *prgb,
		 const rVector3d &eye, const rVector3d &lookat,
		 const rVector3d &up, Real fov, const rVector3d &light,
		 const Real ambiant, const Array<rVector3d> &pts,
		 const Array<Array<int> > &polys,
		 const Array<Array<Real> > &colors) {
  rVector3d front = lookat - eye;
  front.normalize();
  rVector3d side = front ^ up;
  side.normalize();
  rVector3d myup = side ^ front;
  myup.normalize();
  rMatrix3d view;
  view.set_row(0, side);
  view.set_row(1, myup);
  view.set_row(2, front);
  Real dtoplane = 1. / tan((M_PI / 180.) * fov / 2.);
  Array<rVector2d> trans_pts(pts.get_size());
  Array<Real> dpts(pts.get_size());
  for (int i = 0; i < pts.get_size(); i++) {
    rVector3d tpts = view * (pts(i) - eye);
    dpts(i) = tpts(2);
    trans_pts(i)(0) = tpts(0) * dtoplane / tpts(2);
    trans_pts(i)(1) = tpts(1) * dtoplane / tpts(2);
  }
  Array<int> order(polys.get_size());
  Array<Real> dpolys(polys.get_size());
  int nb_back = 0;
  for (int i = 0; i < polys.get_size(); i++) {
    order(i) = i;
    dpolys(i) = 0.;
    for (int j = 0; j < polys(i).get_size(); j++) {
      dpolys(i) += dpts(polys(i)(j));
      if (dpts(polys(i)(j)) < 0) {
	dpolys(i) = -1.;
	nb_back++;
	break;
      }
    }
    dpolys(i) /= (Real)(polys(i).get_size());
  }

  OrderByLookup<Real, IsGreater> greater(dpolys);
  sort_array(&order, greater);

  Array<Real> shades(polys.get_size());
  for (int i = 0; i < shades.get_size(); i++) {
    if (polys(i).get_size() < 3) {
      shades(i) = ambiant;
    } else {
      rVector3d a = pts(polys(i)(1)) - pts(polys(i)(0));
      rVector3d b = pts(polys(i)(2)) - pts(polys(i)(0));
      rVector3d normal = a ^ b;
      normal.normalize();
      if (normal * front > 0.) {
	normal = -normal;
      }
      shades(i) = min(max(0., ambiant + (normal * light) / 2.), 1.);
    }
  }

  ptri->resize(polys.get_size() - nb_back);
  prgb->resize(polys.get_size() - nb_back);
  for (int i = 0; i < ptri->get_size(); i++) {
    (*ptri)(i).resize(polys(order(i)).get_size());
    for (int j = 0; j < (*ptri)(i).get_size(); j++) {
      (*ptri)(i)(j) = trans_pts(polys(order(i))(j));
    }
    (*prgb)(i) = colors(order(i));
    product(&((*prgb)(i)), (*prgb)(i), shades(order(i)));
  }
}

void write_gnuplot(ostream &file, const char *cmd,
		   const Array<Array<rVector2d> > &tri,
		   const Array<Array<Real> > &rgb) {
  file << cmd << endl;
  file << "set size square" << endl;
  file << "unset key" << endl;
  file << "unset tics" << endl;
  file << "unset border" << endl;
  file << "set style fill solid 1 noborder" << endl;
  file << "plot [-1:1] [-1:1] '-' w filledcurve lc rgb var" << endl;
  for (int i = 0; i < tri.get_size(); i++) {
    int col = 0;
    for (int c = 0; c < 3; c++) {
      col = col * 256 + (int)(rgb(i)(c));
    }
    for (int j = 0; j < tri(i).get_size(); j++) {
      file << tri(i)(j)(0) << " " << tri(i)(j)(1) << " "
	   << "0x" << setw(6) << setfill('0') << hex << col << endl;
    }
    file << endl;
  }
  file << "e" << endl;
  file << "unset out " << endl;
}

char *helpstring = "";

int main(int argc, char *argv[]) {
  // parsing command line. See getvalue.hh for details;
  char *infilename = "tri.bi";
  rVector3d eye(0., -5, 0.);
  rVector3d lookat(0., 0., 0.);
  rVector3d up(0., 0., 1.);
  rVector3d light(0., sqrt(2.) / 2., -sqrt(2.) / 2.);
  Real ambiant = 0.5;
  Real fov = 30.;
  char *termcmd = "set term postscript eps color ; set out 'test.eps'";
  int dohelp = 0;
  int dummy = 0;
  Array<Real> eyev, lookatv, upv, lightv;
  FixedVector_to_Array(&eyev, eye);
  FixedVector_to_Array(&lookatv, lookat);
  FixedVector_to_Array(&upv, up);
  FixedVector_to_Array(&lightv, light);
  AskStruct options[] = {
      {"", "3D to GNUplot " MAPS_VERSION ", by Axel van de Walle", TITLEVAL,
       NULL},
      {"-if", "Input table file", STRINGVAL, &infilename},
      {"-tc", "set terminal command", STRINGVAL, &termcmd},
      {"-eye", "Eye position", ARRAYRVAL, &eyev},
      {"-at", "Lookat position", ARRAYRVAL, &lookatv},
      {"-up", "Up direction", ARRAYRVAL, &upv},
      {"-fov", "Field of View (in degrees)", REALVAL, &fov},
      {"-lit", "Light direction", ARRAYRVAL, &lightv},
      {"-amb", "Ambiant light level", REALVAL, &ambiant},
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
  Array_to_FixedVector(&eye, eyev);
  Array_to_FixedVector(&lookat, lookatv);
  Array_to_FixedVector(&up, upv);
  Array_to_FixedVector(&light, lightv);

  LinkedList<Array<Real> > ptslist;
  LinkedList<Array<int> > polyslist;
  LinkedList<Array<Real> > colorslist;
  {
    ifstream file(infilename);
    if (!file) ERRORQUIT("Unable to open input file.");
    while (1) {
      Array<Array<Real> > one_pts;
      LinkedList<Array<int> > one_polys;
      Array<Array<Real> > one_colors;
      file >> one_pts;
      // cerr << one_pts.get_size() << endl;
      if (one_pts.get_size() == 0) break;
      file >> one_polys;
      file >> one_colors;
      combine_mesh(&ptslist, &polyslist, one_pts, one_polys);
      LinkedListIterator<Array<int> > ip(one_polys);
      for (int i = 0; ip; ip++, i++) {
	colorslist << new Array<Real>(one_colors(i % one_colors.get_size()));
      }
    }
  }
  Array<Array<Real> > raw_pts;
  Array<rVector3d> pts;
  Array<Array<int> > polys;
  Array<Array<Real> > colors;
  LinkedList_to_Array(&raw_pts, ptslist);
  ArrayArray_to_ArrayFixedVector(&pts, raw_pts);
  LinkedList_to_Array(&polys, polyslist);
  LinkedList_to_Array(&colors, colorslist);

  Array<Array<rVector2d> > tri;
  Array<Array<Real> > rgb;
  make_2dview(&tri, &rgb, eye, lookat, up, fov, light, ambiant, pts, polys,
	      colors);
  write_gnuplot(cout, termcmd, tri, rgb);
}
