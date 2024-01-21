#include <fstream>

#include "arraylist.h"
#include "integer.h"
#include "linalg.h"
#include "vectmac.h"

class AtomCluster {
 public:
  AtomCluster(void) : atom_pos(0), atom_type(0) {}
  Array<rVector3d> atom_pos;
  Array<int> atom_type;
};

int which_atom(const Array<rVector3d> &atom_pos, const rVector3d &pos,
	       Real cutoff = MAXFLOAT, Real *pdist = NULL) {
  int besti = -1;
  Real bestd = MAXFLOAT;
  for (int i = 0; i < atom_pos.get_size(); i++) {
    Real d = norm(atom_pos(i) - pos);
    if (d < bestd) {
      besti = i;
      bestd = d;
    }
  }
  if (bestd > cutoff) {
    besti = -1;
  }
  if (pdist) {
    *pdist = bestd;
  }
  return besti;
}

void enum_triplets(LinkedList<Array<rVector3d> > *trip_list,
		   LinkedList<Array<Real> > *dist_list,
		   const Array<rVector3d> &atom_pos, Real cutoff = MAXFLOAT) {
  iVector3d i;
  for (i(0) = 0; i(0) < atom_pos.get_size(); i(0)++) {
    for (i(1) = 0; i(1) < i(0); i(1)++) {
      for (i(2) = 0; i(2) < i(1); i(2)++) {
	Array<Real> dist(3);
	int j;
	for (j = 0; j < 3; j++) {
	  dist(j) = norm(atom_pos(i((j + 1) % 3)) - atom_pos(i(j)));
	  if (dist(j) > cutoff) break;
	}
	rVector3d v1 = atom_pos(i(1)) - atom_pos(i(0));
	v1.normalize();
	rVector3d v2 = atom_pos(i(2)) - atom_pos(i(0));
	v2.normalize();
	if (j == 3 && norm(v1 ^ v2) > zero_tolerance) {
	  int rot = index_max(dist);
	  int flip = (dist((rot + 1) % 3) > dist((rot + 2) % 3) ? 1 : -1);
	  Array<rVector3d> *ptrip = new Array<rVector3d>(3);
	  Array<Real> *pdist = new Array<Real>(3);
	  for (int k = 0; k < 3; k++) {
	    (*ptrip)(k) = atom_pos(i((3 + rot + flip * k) % 3));
	    (*pdist)(k) = dist((3 + rot + flip * k) % 3);
	  }
	  (*dist_list) << pdist;
	  (*trip_list) << ptrip;
	}
      }
    }
  }
}

void find_overlapping_clusters(LinkedList<Array<rVector3d> > *cluslist,
			       const Array<rVector3d> &seed,
			       const Array<rVector3d> &clus, Real cutoff) {
  LinkedList<Array<int> > sitelist;
  LinkedList<Array<rVector3d> > tripclus;
  LinkedList<Array<Real> > distclus;
  enum_triplets(&tripclus, &distclus, clus);
  Real maxdistclus = 0.;
  LinkedListIterator<Array<Real> > i(distclus);
  for (; i; i++) {
    Real m = max(*i);
    maxdistclus = MAX(maxdistclus, m);
  }
  LinkedList<Array<rVector3d> > tripseed;
  LinkedList<Array<Real> > distseed;
  enum_triplets(&tripseed, &distseed, seed, maxdistclus);

  LinkedListIterator<Array<rVector3d> > iseed(tripseed);
  LinkedListIterator<Array<Real> > dseed(distseed);
  for (; iseed; iseed++, dseed++) {
    LinkedListIterator<Array<rVector3d> > iclus(tripclus);
    LinkedListIterator<Array<Real> > dclus(distclus);
    for (; iclus; iclus++, dclus++) {
      for (int flip = -1; flip <= 1; flip += 2) {
	for (int rot = 0; rot < 3; rot++) {
	  int side = 0;
	  for (; side < 3; side++) {
	    if (fabs((*dseed)(side) - (*dclus)((3 + flip * side + rot) % 3)) >
		cutoff * 2)
	      break;
	  }
	  if (side == 3) {
	    rMatrix3d matseed;
	    rVector3d v1 = (*iseed)(1) - (*iseed)(0);
	    v1.normalize();
	    rVector3d v2 = (*iseed)(2) - (*iseed)(0);
	    v2.normalize();
	    rVector3d v3 = v1 ^ v2;
	    v3.normalize();
	    matseed.set_column(0, v1);
	    matseed.set_column(1, v2);
	    matseed.set_column(2, v3);
	    rMatrix3d matclus;
	    v1 = (*iclus)((3 + flip * 1 + rot) % 3) -
		 (*iclus)((3 + flip * 0 + rot) % 3);
	    v1.normalize();
	    v2 = (*iclus)((3 + flip * 2 + rot) % 3) -
		 (*iclus)((3 + flip * 0 + rot) % 3);
	    v2.normalize();
	    for (int mirror = -1; mirror <= 1; mirror += 2) {
	      rVector3d v3 = mirror * v1 ^ v2;
	      v3.normalize();
	      matclus.set_column(0, v1);
	      matclus.set_column(1, v2);
	      matclus.set_column(2, v3);
	      rMatrix3d m = matseed * (!matclus);
	      rMatrix3d mmt = m * (~m);
	      rMatrix3d invstrain;
	      pow(&invstrain, mmt, -0.5);
	      rMatrix3d matrot = invstrain * m;
	      rVector3d cmassseed =
		  ((*iseed)(0) + (*iseed)(1) + (*iseed)(2)) / 3;
	      rVector3d cmassclus =
		  ((*iclus)(0) + (*iclus)(1) + (*iclus)(2)) / 3;
	      rVector3d trans = cmassseed - matrot * cmassclus;
	      Array<rVector3d> newclus(clus.get_size());
	      Array<int> newsite(clus.get_size());
	      int i;
	      for (i = 0; i < clus.get_size(); i++) {
		rVector3d pos = matrot * clus(i) + trans;
		newclus(i) = pos;
		newsite(i) = which_atom(seed, pos, cutoff);
		int ii;
		for (ii = 0; ii < i; ii++) {
		  if (newsite(i) == newsite(ii)) break;
		}
		if (ii != i) break;
	      }
	      if (i == clus.get_size()) {
		LinkedListIterator<Array<int> > jsite(sitelist);
		for (; jsite; jsite++) {
		  if (newsite == (*jsite)) break;
		}
		if (!jsite) {
		  (*cluslist) << new Array<rVector3d>(newclus);
		  sitelist << new Array<int>(newsite);
		}
	      }
	    }
	  }
	}
      }
    }
  }
}

void find_periodicity(Array<rVector3d> *p_period, const AtomCluster &clus,
		      Real cutoff) {
  int invoidcutoff = clus.atom_pos.get_size() ? ? ? ? ;
  for (int i = 0; i < clus.atom_pos.get_size(); i++) {
    for (int j = 0; j < i; j++) {
      rVector v = clus.atom_pos(j) - clus.atom_pos(i);
      Array<int> newsite(clus.atom_pos.get_size());
      int invoid = 0;
      int k;
      for (k = 0; k < clus.atom_pos.get_size(); k++) {
	int l = which_atom(clus.atom_pos, clus.atom_pos(k) + v, cutoff);
	if (l == -1) {
	  invoid++;
	} else {
	  if (clus.atom_type(k) != clus.atom_type(l)) break;
	  int kk;
	  for (kk = 0; kk < k; ii++) {
	    if (l == newsite(kk)) break;
	  }
	  newsite(k) = l;
	  if (kk != k) break;
	}
      }
      if (k == clus.atom_pos.get_size()) {
	if (invoid < invoidcutoff) {
	  period_list << new rVector3d(v);
	}
      }
    }
  }
}

void read_clus(Array<rVector3d> *pclus, istream &file) {
  LinkedList<rVector3d> clus_list;
  while (!file.eof()) {
    rVector3d v(MAXFLOAT, MAXFLOAT, MAXFLOAT);
    file >> v;
    if (v(2) == MAXFLOAT) break;
    clus_list << new rVector3d(v);
  }
  LinkedList_to_Array(pclus, clus_list);
}

int main(void) {
  Real cutoff = 0.2;
  Array<rVector3d> seed;
  ifstream seedfile("seed.in");
  read_clus(&seed, seedfile);
  Array<rVector3d> clusorg;
  ifstream clusfile("clus.in");
  read_clus(&clusorg, clusfile);
  LinkedList<Array<rVector3d> > cluslist;
  find_overlapping_clusters(&cluslist, seed, clusorg, cutoff);
  LinkedListIterator<Array<rVector3d> > i(cluslist);
  cout.setf(ios::fixed);
  cout.precision(3);
  int cnt = 0;
  cout << cluslist.get_size() << endl;
  for (; i; i++) {
    Array<rVector3d> &clus = (*i);
    int overlap = 0;
    for (int j = 0; j < clus.get_size(); j++) {
      if (which_atom(seed, clus(j), cutoff) != -1) {
	overlap++;
      }
    }
    cout << overlap << endl;
    if (overlap < clus.get_size()) {
      ofstream file("clus.xyz");
      file << (seed.get_size() + clus.get_size() - overlap) << endl;
      file << "Title" << endl;
      for (int j = 0; j < seed.get_size(); j++) {
	if (which_atom(clus, seed(j), cutoff) == -1) {
	  file << "S " << seed(j) << endl;
	} else {
	  file << "O " << seed(j) << endl;
	}
      }
      for (int j = 0; j < clus.get_size(); j++) {
	if (which_atom(seed, clus(j), cutoff) == -1) {
	  file << "C " << clus(j) << endl;
	}
      }
      system("rasmol -xyz clus.xyz");
      cnt++;
      if (cnt == 5) break;
    }
  }
}
