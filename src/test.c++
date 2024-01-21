//#include <strstream.h>
#include <fstream>

#include "chull.h"
#include "integer.h"
#include "linalg.h"

int main(void) {
  int sd;
  cin >> sd;
  rndseed(sd);
  zero_tolerance = 1e-5;
  int d = 3;
  int n = 55;
  Array<Array<Real> > x(n);
  for (int i = 0; i < n; i++) {
    x(i).resize(d);
    for (int j = 0; j < d; j++) {
      x(i)(j) = uniform01();
    }
  }
  Array<Real> ground(d);
  zero_array(&ground);
  //  ground(d-1)=-1;
  LinkedList<PolytopeFace> hull;
  calc_convex_hull(&hull, x, ground);
  ofstream pfile("pts.out");
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < d; j++) {
      pfile << x(i)(j) << " ";
    }
    pfile << endl;
  }
  ofstream hfile("hull.out");
  LinkedListIterator<PolytopeFace> h(hull);
  for (; h; h++) {
    for (int j = 0; j <= h->pts.get_size(); j++) {
      int i = h->pts(j % (h->pts.get_size()));
      for (int j = 0; j < d; j++) {
	hfile << x(i)(j) << " ";
      }
      hfile << endl;
    }
    hfile << endl << endl;
  }
}
