#include "meshcalc.h"

#include <fstream>

#include "lstsqr.h"

void debug_mesh(const LinkedList<Array<int> > &polylist,
		const Array<Array<Real> > &pts) {
  ofstream file("mesh.out");
  LinkedListIterator<Array<int> > it(polylist);
  for (; it; it++) {
    ofstream last("lastpoly.out");
    for (int i = 0; i < it->get_size(); i++) {
      for (int ii = 0; ii < i; ii++) {
	for (int j = 0; j < pts((*it)(0)).get_size(); j++) {
	  file << pts((*it)(i))(j) << " ";
	  last << pts((*it)(i))(j) << " ";
	}
	file << endl;
	last << endl;
	for (int j = 0; j < pts((*it)(0)).get_size(); j++) {
	  file << pts((*it)(ii))(j) << " ";
	  last << pts((*it)(ii))(j) << " ";
	}
	file << endl;
	file << endl;
	file << endl;
	last << endl;
	last << endl;
	last << endl;
      }
    }
  }
}

void find_perp_to_vector_array(Array<Real> *pperp,
			       const Array<Array<Real> > &v) {
  int dim = v(0).get_size();
  pperp->resize(dim);
  Array2d<Real> B(dim, dim), iB;
  for (int i = 0; i < dim - 1; i++) {
    for (int j = 0; j < dim; j++) {
      B(i, j) = v(i)(j);
    }
  }
  for (int j = 0; j < dim; j++) {
    B(dim - 1, j) = uniform01() - 0.5;
  }
  invert_matrix(&iB, B);
  for (int i = 0; i < dim; i++) {
    (*pperp)(i) = iB(i, dim - 1);
  }
  normalize(pperp);
}

int calc_circumsphere(Array<Real> *pcenter, Real *pradius2,
		      const Array<Array<Real> > &pts) {
  int dim = pts(0).get_size();
  int nbpts = pts.get_size();
  Array2d<Real> A(dim, dim);
  Array<Real> v(dim);
  for (int i = 0; i < nbpts - 1; i++) {
    Array<Real> delta;
    diff(&delta, pts(i + 1), pts(i));
    set_row(&A, delta, i);
    v(i) = (norm2(pts(i + 1)) - norm2(pts(i))) / 2.;
  }
  if (nbpts < dim + 1) {
    Array2d<Real> B, iB;
    Array<Real> u(dim);
    B = A;
    for (int i = 0; i < dim; i++) {
      B(dim - 1, i) = uniform01() - 0.5;
    }
    invert_matrix(&iB, B);
    for (int i = 0; i < dim; i++) {
      u(i) = iB(i, dim - 1);
      A(dim - 1, i) = u(i);
    }
    v(dim - 1) = inner_product(u, pts(0));
  } else if (nbpts != dim + 1) {
    ERRORQUIT("wrong dimensions in calc_circumsphere");
  }
  solve_linsys(&A, &v);
  *pcenter = v;
  Array<Real> r;
  diff(&r, pts(0), *pcenter);
  *pradius2 = norm2(r);
  return 1;
}

void debug_side(const LinkedList<SideFlag> &edgelist,
		const Array<Array<Real> > &pts) {
  ofstream file("normal.out");
  ofstream filedone("done.out");
  LinkedListIterator<SideFlag> it(edgelist);
  for (; it; it++) {
    Array<Real> center;
    center = pts(it->ipts(0));
    for (int i = 1; i < it->ipts.get_size(); i++) {
      sum(&center, center, pts(it->ipts(i)));
    }
    product(&center, center, 1. / (Real)(it->ipts.get_size()));
    for (int j = 0; j < center.get_size(); j++) {
      filedone << center(j) << " ";
    }
    filedone << it->done << endl;
    if (it->forward.get_size() != 0) {
      Array<Real> endarrow;
      sum(&endarrow, center, it->forward);
      for (int j = 0; j < center.get_size(); j++) {
	file << center(j) << " ";
      }
      file << endl;
      for (int j = 0; j < center.get_size(); j++) {
	file << endarrow(j) << " ";
      }
      file << endl;
      file << endl;
      file << endl;
    }
  }
}

PointGridIndex::PointGridIndex(const Array<Array<Real> > &pts)
    : low(), boxsize(0), indices() {
  if (pts.get_size() == 0) return;
  int dim = pts(0).get_size();
  low = pts(0);
  Array<Real> high(pts(0));
  for (int i = 0; i < dim; i++) {
    for (int j = 0; j < pts.get_size(); j++) {
      if (pts(j)(i) < low(i)) {
	low(i) = pts(j)(i);
      }
      if (pts(j)(i) > high(i)) {
	high(i) = pts(j)(i);
      }
    }
  }
  Array<Real> range;
  diff(&range, high, low);
  boxsize =
      max(range) / pow((Real)pts.get_size(), 1. / (Real)dim) * 0.5;  // adjust;
  Array<int> size(dim);
  for (int i = 0; i < dim; i++) {
    size(i) = (int)ceil((high(i) + zero_tolerance - low(i)) / boxsize);
    if (size(i) == 0) {
      size(i) = 1;
    }
  }
  indices.resize(size);

  // cerr << "grid index" << endl;
  // cerr << low << endl;
  // cerr << boxsize << endl;
  for (int j = 0; j < pts.get_size(); j++) {
    Array<int> box(dim);
    for (int i = 0; i < dim; i++) {
      box(i) = (int)floor((pts(j)(i) - low(i)) / boxsize);
    }
    indices(box) << new int(j);
    // cerr << box << endl;
  }
  // cerr << "endgrid" << endl;
}

int PointFromCenterIteratorGridIndex::fill_buffer(int new_boxrad) {
  // cerr << new_boxrad << endl;
  LinkedList<RealIntPair> distindexlist;
  for (int i = cur_index; i < distindex.get_size(); i++) {
    distindexlist << new RealIntPair(distindex(i));
  }
  int dim = pts(0).get_size();
  Array<int> box_lo(dim);
  Array<int> box_hi(dim);
  int nobump = 0;
  for (int i = 0; i < dim; i++) {
    int uplim = grid_index.indices.get_size()(i) - 1;
    box_lo(i) = box_center(i) - new_boxrad;
    if (box_lo(i) <= 0) {
      box_lo(i) = 0;
    } else {
      nobump = 1;
    }
    if (box_lo(i) > uplim) {
      box_lo(i) = uplim;
    }
    box_hi(i) = box_center(i) + new_boxrad;
    if (box_hi(i) >= uplim) {
      box_hi(i) = uplim;
    } else {
      nobump = 1;
    }
    if (box_hi(i) < 0) {
      box_hi(i) = 0;
    }  //??;
    //      cerr << box_lo(i) << " " << box_hi(i) << endl;
  }
  //    cerr << "== " << cur_boxrad << " " << new_boxrad << endl;
  MultiDimIterator<Array<int> > b(box_lo, box_hi);
  for (; b; b++) {
    Array<int> db;
    diff(&db, (Array<int> &)b, box_center);
    for (int i = 0; i < db.get_size(); i++) {
      db(i) = fabs(db(i));
    }
    if (max(db) > cur_boxrad) {
      LinkedListIterator<int> it(grid_index.indices(b));
      for (; it; it++) {
	Array<Real> v;
	diff(&v, pts(*it), center);
	distindexlist << new RealIntPair(norm(v), *it);
      }
    }
  }
  LinkedList_to_Array(&distindex, distindexlist);
  sort_array(&distindex, OrderBy_r<RealIntPair>());
  cur_index = 0;
  cur_boxrad = new_boxrad;
  cur_rad = (Real)cur_boxrad * grid_index.boxsize;
  // cerr << "=== " << cur_rad << endl;
  return (nobump);
}

PointFromCenterIteratorGridIndex::PointFromCenterIteratorGridIndex(
    const Array<Array<Real> > &_pts, const PointGridIndex &_grid_index,
    const Array<Real> &_center)
    : pts(_pts),
      center(_center),
      box_center(_center.get_size()),
      distindex(0),
      cur_index(0),
      grid_index(_grid_index),
      cur_rad(0.),
      cur_boxrad(0) {
  if (pts.get_size() == 0) return;
  int dim = pts(0).get_size();
  for (int i = 0; i < dim; i++) {
    box_center(i) =
	(int)floor((center(i) - grid_index.low(i)) / grid_index.boxsize);
  }
  init();
}

void PointFromCenterIteratorGridIndex::init(void) {
  cur_boxrad = -1;
  fill_buffer(1);
  update_if_needed();
}

PointFromCenterIteratorGridIndex::operator void *() {
  return (void *)(cur_index < distindex.get_size() ? &cur_index : NULL);
}

int PointFromCenterIteratorGridIndex::operator()(void) {
  return distindex(cur_index).i;
}

void PointFromCenterIteratorGridIndex::update_if_needed(void) {
  while (1) {
    if (cur_index < distindex.get_size()) {
      // cerr << "<\n";
      if (distindex(cur_index).r <= cur_rad) {
	// cerr << "<=\n";
	break;
      } else {
	// cerr << "add " << distindex.get_size() << " " << cur_index << " " <<
	// center << " " << pts(distindex(cur_index).i) << " " <<
	// distindex(cur_index).r << " " << cur_rad << "\n";
	int addbox =
	    (int)ceil((distindex(cur_index).r - cur_rad) / grid_index.boxsize);
	if (!fill_buffer(cur_boxrad + max(1, addbox))) {
	  break;
	}
      }
    } else {
      // cerr << "+1\n";
      if (!fill_buffer(cur_boxrad + 1)) {
	break;
      }
    }
  }
}

void PointFromCenterIteratorGridIndex::operator++(int) {
  cur_index++;
  update_if_needed();
  /*
    while (1) {
    int good;
    if (cur_index<distindex.get_size()) {
    if (distindex(cur_index).r>cur_rad)  {good=0;}
    else {good=1;}
    }
    else {good=0;}
    if (!good) {
    if (!fill_buffer(cur_boxrad+1)) {break;}
    // cerr << "s=" << distindex.get_size() << endl;
    }
    else {break;}
    }
  */
}

PointFromCenterIterator::PointFromCenterIterator(
    const Array<Array<Real> > &_pts, const Array<Real> &_center)
    : pts(_pts), center(_center), distindex(_pts.get_size()), curindex(0) {
  for (int i = 0; i < pts.get_size(); i++) {
    Array<Real> d;
    diff(&d, pts(i), center);
    distindex(i).r = norm(d);
    distindex(i).i = i;
  }
  sort_array(&distindex, OrderBy_r<RealIntPair>());
}

void PointFromCenterIterator::init(void) { curindex = 0; }

PointFromCenterIterator::operator void *() {
  return (void *)(curindex < distindex.get_size() ? &curindex : NULL);
}

int PointFromCenterIterator::operator()(void) { return distindex(curindex).i; }

void PointFromCenterIterator::operator++(int) { curindex++; }

int find_point_to_add(const SideFlag &side, const Array<Array<Real> > &pts,
		      const PointGridIndex &grid_index, Real maxr = MAXFLOAT) {
  Array<Real> midside(pts(side.ipts(0)));
  for (int i = 1; i < side.ipts.get_size(); i++) {
    sum(&midside, midside, pts(side.ipts(i)));
  }
  product(&midside, midside, 1. / (Real)(side.ipts.get_size()));
  Array<Array<Real> > polytope(side.ipts.get_size() + 1);
  for (int i = 0; i < side.ipts.get_size(); i++) {
    polytope(i + 1) = pts(side.ipts(i));
    // cerr << "side[" << i << "]=" << side.ipts(i) << endl;
  }
  LinkedList<int> good_pts;
  // cerr << "midside=" << midside << endl;
  PointFromCenterIteratorGridIndex itpts(pts, grid_index, midside);
  for (; itpts; itpts++) {
    // cerr << "here\n";
    // cerr << "itpts=" << itpts() << endl;
    Array<Real> delta;
    diff(&delta, pts(itpts()), midside);
    if (norm(delta) > maxr) {
      break;
    }
    if (inner_product(side.forward, delta) > 0 &&
	!is_in_array(side.ipts, itpts())) {
      polytope(0) = pts(itpts());
      // cerr << polytope << endl;
      Array<Real> center;
      Real radius2;
      calc_circumsphere(&center, &radius2, polytope);
      int success = 0;
      // cerr << "center=" << center << endl;//DEBUG
      // cerr << "r= " << radius2 << endl;
      // cerr << "<pow: " << (radius2<pow(zero_tolerance,-2.)) << endl;
      if (radius2 < pow(zero_tolerance, -2.)) {
	// cerr << "do\n";
	PointFromCenterIteratorGridIndex jtpts(pts, grid_index, center);
	for (; jtpts; jtpts++) {
	  if (is_in_array(side.ipts, jtpts()) || jtpts() == itpts()) {
	    success = 1;
	    break;
	  }
	  Array<Real> delta;
	  diff(&delta, pts(jtpts()), center);
	  // cerr << norm(delta) << "<>" << sqrt(radius2) << endl;
	  if (norm2(delta) < radius2) break;
	}
      }
      if (success) {
	good_pts << new int(itpts());
	// cerr << "picked" << endl;
	// break;
      }
    }
  }
  // cerr << "good=" << good_pts.get_size() << endl;
  LinkedListIterator<int> ig(good_pts);
  int besti = -1;
  Real bestproj = 2.;
  for (; ig; ig++) {
    Array<Real> delta;
    diff(&delta, pts(*ig), midside);
    Real proj = fabs(inner_product(side.normal, delta) / norm(delta));
    // cerr << proj << endl;
    if (proj < bestproj) {
      besti = *ig;
      bestproj = proj;
    }
  }
  return besti;
  //  if (!itpts) {return -1;}
  //  return itpts();
}

void calc_forward_vector(Array<Real> *pforward, const Array<int> &sideipts,
			 int ptipts, const Array<Array<Real> > &pts) {
  int dim = pts(0).get_size();
  Array<Array<Real> > edge(dim - 1);
  for (int j = 0; j < sideipts.get_size() - 1; j++) {
    diff(&(edge(j)), pts(sideipts(j + 1)), pts(sideipts(0)));
  }
  // cerr << "edge=" << edge << endl;
  if (sideipts.get_size() < dim) {
    Array<Array<Real> > planeedge(dim - 1);
    for (int j = 0; j < dim - 1; j++) {
      // cerr << "j=" << j << endl;
      // cerr << sideipts(j) << endl;
      diff(&(planeedge(j)), pts(sideipts(j)), pts(ptipts));
    }
    find_perp_to_vector_array(&(edge(dim - 2)), planeedge);
  }
  // cerr << "here2\n";
  find_perp_to_vector_array(pforward, edge);
  // cerr << "here3\n";
  Array<Real> testv(dim);
  diff(&testv, pts(ptipts), pts(sideipts(0)));
  if (inner_product(testv, *pforward) > 0) {
    product(pforward, *pforward, -1.);
  }
}

void calc_normal_vector(Array<Real> *pnormal, const Array<int> &sideipts,
			int ptipts, const Array<Array<Real> > &pts) {
  int dim = pts(0).get_size();
  if (sideipts.get_size() < dim) {
    Array<Array<Real> > planeedge(dim - 1);
    for (int j = 0; j < dim - 1; j++) {
      diff(&(planeedge(j)), pts(sideipts(j)), pts(ptipts));
    }
    find_perp_to_vector_array(pnormal, planeedge);
    normalize(pnormal);
  } else {
    pnormal->resize(dim);
    zero_array(pnormal);
  }
}

int add_point_to_mesh(LinkedList<Array<int> > *ppolylist,
		      LinkedList<SideFlag> *psidelist, int pttoadd,
		      SideFlag *pside, const Array<Array<Real> > &pts) {
  int good = 1;
  LinkedList<SideFlag> newsidelist;
  for (int i = 0; i < pside->ipts.get_size(); i++) {
    SideFlag *pnewside = new SideFlag();
    pnewside->ipts = pside->ipts;
    pnewside->ipts(i) = pttoadd;
    sort_array(&(pnewside->ipts));
    pnewside->done = 1;
    calc_forward_vector(&pnewside->forward, pnewside->ipts, pside->ipts(i),
			pts);
    calc_normal_vector(&pnewside->normal, pnewside->ipts, pside->ipts(i), pts);
    LinkedListIterator<SideFlag> itexistside(*psidelist);
    for (; itexistside; itexistside++) {
      if (itexistside->ipts == pnewside->ipts) break;
    }
    newsidelist << pnewside;
    if (itexistside) {
      if (itexistside->done != 1 && itexistside->done != -1) {
	good = 0;
	// cerr << "Not adding side: done=" << itexistside->done << endl;
      }
    }
  }
  if (good) {
    LinkedListIterator<SideFlag> it(newsidelist);
    while (it) {
      LinkedListIterator<SideFlag> itexistside(*psidelist);
      for (; itexistside; itexistside++) {
	if (itexistside->ipts == it->ipts) break;
      }
      if (!itexistside) {
	// cerr << "side added\n";
	(*psidelist) << newsidelist.detach(it);
      } else {
	itexistside->done = 2;
	itexistside->forward.resize(0);
	it++;
      }
    }
    pside->done = 2;
    pside->forward.resize(0);
    Array<int> *ppoly = new Array<int>(pside->ipts.get_size() + 1);
    (*ppoly)(0) = pttoadd;
    for (int i = 0; i < pside->ipts.get_size(); i++) {
      (*ppoly)(i + 1) = pside->ipts(i);
    }
    // cerr << "POLY" << (*ppoly) << endl;
    sort_array(ppoly);
    (*ppolylist) << ppoly;
    return 1;
  } else {
    pside->done = -1;
    pside->forward.resize(0);
    return 0;
  }
}

void create_mesh(LinkedList<Array<int> > *ppolylist,
		 const Array<Array<Real> > &pts, int fulldim, Real maxr) {
  if (pts.get_size() == 0) {
    return;
  }
  if (pts.get_size() < pts(0).get_size() + fulldim) {
    return;
  }

  PointGridIndex grid_index(pts);
  LinkedList<SideFlag> sidelist;
  int dim = pts(0).get_size();
  if (pts.get_size() < dim + fulldim) return;

  Array<int> meshed(pts.get_size());
  zero_array(&meshed);
  int nbmeshed = 0;
  // cerr << "Start\n";
  while (1) {
    LinkedListIterator<Array<int> > ip(*ppolylist);
    for (; ip; ip++) {
      for (int j = 0; j < ip->get_size(); j++) {
	if (meshed((*ip)(j)) == 0) {
	  meshed((*ip)(j)) = 1;
	  nbmeshed++;
	}
      }
    }
    // cerr << "nbm= " << nbmeshed << "/" << pts.get_size() << endl;
    if (nbmeshed == pts.get_size()) break;

    int startpt = 0;
    while (meshed(startpt) == 1) {
      startpt++;
    }

    Array<int> *ppoly = new Array<int>(dim + fulldim);
    Array<Real> acc(pts(startpt));
    (*ppoly)(0) = startpt;
    meshed(startpt) = 1;
    nbmeshed++;
    // cerr << "m=" << endl << meshed << endl;
    int problem = 0;
    for (int n = 1; n < ppoly->get_size(); n++) {
      Array<Real> avg;
      product(&avg, acc, 1. / (Real)n);
      // cerr << "avg= " << avg << endl;
      PointFromCenterIteratorGridIndex itpts(pts, grid_index, avg);
      while (itpts) {
	//	if (meshed(itpts())==0) {
	int i;
	for (i = 0; i < n; i++) {
	  if ((*ppoly)(i) == itpts()) break;
	}
	if (i == n) break;
	// }
	itpts++;
      }
      if (itpts) {
	Array<Real> dr;
	diff(&dr, pts(itpts()), avg);
	if (norm(dr) <= maxr) {
	  (*ppoly)(n) = itpts();
	  sum(&acc, acc, pts(itpts()));
	  // meshed(itpts())=1;
	  // nbmeshed++;
	} else {
	  // cerr << "too close\n";
	  problem = 1;
	  break;
	}
      } else {
	// cerr << "no itpts\n";
	problem = 1;
	break;
      }
    }
    if (!problem) {
      Array<Real> center;
      Real radius2;
      Array<Array<Real> > mypts(ppoly->get_size());
      for (int i = 0; i < mypts.get_size(); i++) {
	mypts(i) = pts((*ppoly)(i));
      }
      calc_circumsphere(&center, &radius2, mypts);
      // cerr << "center2= " << center << endl;
      // cerr << fabs(radius2) << endl;
      // cerr <<  (fabs(radius2)>pow(zero_tolerance,-2.)) << endl;
      if (radius2 < pow(zero_tolerance, -2.)) {
	PointFromCenterIteratorGridIndex jtpts(pts, grid_index, center);
	if (jtpts) {
	  Array<Real> delta;
	  diff(&delta, pts(jtpts()), center);
	  // cerr << norm2(delta) << " " << radius2 << " " <<
	  // (norm2(delta)-radius2) << endl;
	  if (norm2(delta) < radius2) {
	    if (!is_in_array(*ppoly, jtpts())) {
	      problem = 1;
	      // cerr << "circum not ok\n";
	    }
	  }
	}
      } else {
	problem = 1;
	// cerr << "circle too big\n";
      }
    }

    if (problem) {
      delete ppoly;
    } else {
      //      cerr << "go!\n";
      (*ppolylist) << ppoly;

      // cerr << "poly= ";
      // cerr << *ppoly << endl;

      for (int n = 0; n < ppoly->get_size(); n++) {
	SideFlag *pside = new SideFlag();
	pside->ipts.resize(ppoly->get_size() - 1);
	for (int i = 0; i < ppoly->get_size() - 1; i++) {
	  pside->ipts(i) = (*ppoly)(i + (i >= n ? 1 : 0));
	}
	sort_array(&(pside->ipts));
	calc_forward_vector(&(pside->forward), pside->ipts, (*ppoly)(n), pts);
	calc_normal_vector(&(pside->normal), pside->ipts, (*ppoly)(n), pts);
	pside->done = 1;
	sidelist << pside;
	// cerr << "side= ";
	// cerr << pside->ipts << endl;
	// cerr << pside->forward << endl;
      }
      // debug_mesh(*ppolylist,pts);
      // debug_side(sidelist,pts);
      // cerr << "done 1st poly\n";

      int didsome;
      do {
	LinkedListIterator<SideFlag> itside(sidelist);
	didsome = 0;
	for (; itside; itside++) {
	  if (itside->done == 1) {
	    // cerr << itside->ipts << endl;
	    int pttoadd = find_point_to_add(*itside, pts, grid_index, maxr);
	    // cerr << "add=" << pttoadd << endl;
	    if (pttoadd == -1) {
	      itside->done = -1;
	    } else {
	      // cerr << "add= " << pttoadd << endl;
	      add_point_to_mesh(ppolylist, &sidelist, pttoadd, itside, pts);
	      // debug_mesh(*ppolylist,pts);
	      // debug_side(sidelist,pts);
	      Real tmp;
	      didsome = 1;
	    }
	  }
	}
      } while (didsome);
    }
  }
}

int remove_nearby(Array<Array<Real> > *ppts, const Array<Array<Real> > &orgpts,
		  Real cutoff) {
  int ret = 1;
  Real cutoff2 = sqr(cutoff);
  LinkedList<Array<Real> > list;
  for (int i = 0; i < orgpts.get_size(); i++) {
    LinkedListIterator<Array<Real> > j(list);
    for (; j; j++) {
      Array<Real> delta;
      diff(&delta, orgpts(i), *j);
      if (inner_product(delta, delta) < cutoff2) {
	ret = 0;
	break;
      }
    }
    if (!j) {
      list << new Array<Real>(orgpts(i));
    }
  }
  LinkedList_to_Array(ppts, list);
  return ret;
}

class TupleDoneInfo {
 public:
  Array<int> orgindex;
  int newindex;
  const Array<Real> *pprojpt;
  TupleDoneInfo(const Array<int> &_orgindex, int _newindex,
		const Array<Real> *_pprojpt)
      : orgindex(_orgindex), newindex(_newindex), pprojpt(_pprojpt) {}
};

void cross_section_simplexes(Array<Array<Real> > *pcutpts,
			     LinkedList<Array<int> > *pcutpoly,
			     const Array<Array<Real> > &axes,
			     const Array<Real> &org,
			     const Array<Array<Real> > &pts,
			     const LinkedList<Array<int> > &poly) {
  Array<Array<Real> > dummy;
  cross_section_simplexes(pcutpts, pcutpoly, NULL, axes, org, pts, poly, dummy);
}

void cross_section_simplexes(
    Array<Array<Real> > *pcutpts, LinkedList<Array<int> > *pcutpoly,
    Array<Array<Real> > *pcutcolors, const Array<Array<Real> > &axes,
    const Array<Real> &org, const Array<Array<Real> > &pts,
    const LinkedList<Array<int> > &poly, const Array<Array<Real> > &colors) {
  int dim = pts(0).get_size();
  int nax = axes.get_size();
  int nlam = dim - nax + 1;

  Array2d<Real> normal(dim - nax, dim);
  Array<Real> shift(dim - nax);

  Array<Array<Real> > basis;
  for (int i = 0; i < nax; i++) {
    build_ortho_basis(&basis, NULL, axes(i));
  }
  for (int i = 0; i < dim - nax;) {
    Array<Real> trialv(dim);
    for (int j = 0; j < dim; j++) {
      trialv(j) = uniform01() - 0.5;
    }
    if (build_ortho_basis(&basis, NULL, trialv)) {
      set_row(&normal, basis(basis.get_size() - 1), i);
      i++;
    }
  }
  product(&shift, normal, org);

  int curptsindex = 0;
  LinkedList<Array<Real> > cutptslist;
  LinkedList<Array<Real> > colorslist;
  Array<LinkedList<TupleDoneInfo> > tuple_done(pts.get_size());
  int icol = 0;
  LinkedListIterator<Array<int> > ipoly(poly);
  for (; ipoly; ipoly++, icol++) {
    LinkedList<int> polyindexlist;
    LinkedList<const Array<Real> *> polyptslist;
    NTupleIterator it(nlam, dim);
    for (; it; it++) {
      Array<int> ipts(nlam);
      for (int i = 0; i < ipts.get_size(); i++) {
	ipts(i) = (*ipoly)(((Array<int> &)it)(i));
      }
      sort_array(&ipts, TrivLessThan<int>());
      LinkedListIterator<TupleDoneInfo> ituple(tuple_done(ipts(0)));
      for (; ituple; ituple++) {
	int j;
	for (j = 1; j < ipts.get_size(); j++) {
	  if (ipts(j) != ituple->orgindex(j)) break;
	}
	if (j == ipts.get_size()) break;
      }
      if (ituple) {
	polyindexlist << new int(ituple->newindex);
	polyptslist << new (const Array<Real> *)(ituple->pprojpt);
      } else {
	Array<Real> c(nlam);
	Array<int> signflip(normal.get_size()(0));
	zero_array(&signflip);
	Array<Array<Real> > proj(ipts.get_size());
	for (int i = 0; i < proj.get_size(); i++) {
	  product(&proj(i), normal, pts(ipts(i)));
	  diff(&proj(i), proj(i), shift);
	  for (int j = 0; j < signflip.get_size(); j++) {
	    if (sgn(proj(i)(j)) != sgn(proj(0)(j))) {
	      signflip(j) = 1;
	    }
	  }
	}
	if (max(signflip) > 0) {
	  Array2d<Real> A(nlam, nlam);
	  zero_array(&c);
	  c(0) = 1.;
	  for (int i = 0; i < nlam; i++) {
	    A(0, i) = 1.;
	    for (int j = 1; j < nlam; j++) {
	      A(j, i) = proj(i)(j - 1);
	    }
	  }
	  if (solve_linsys(&A, &c)) {
	    if (min(c) >= 0 && max(c) <= 1) {
	      Array<Real> newpt(dim);
	      zero_array(&newpt);
	      for (int i = 0; i < nlam; i++) {
		Array<Real> tmp;
		product(&tmp, pts(ipts(i)), c(i));
		sum(&newpt, newpt, tmp);
	      }

	      Array<Real> *pnewprojpt = new Array<Real>(nax);
	      {
		Array<Real> shnewpt;
		diff(&shnewpt, newpt, org);
		for (int j = 0; j < nax; j++) {
		  (*pnewprojpt)(j) = inner_product(shnewpt, axes(j));
		}
	      }

	      cutptslist << pnewprojpt;
	      TupleDoneInfo *pnew_tuple =
		  new TupleDoneInfo(ipts, curptsindex, pnewprojpt);
	      curptsindex++;
	      tuple_done(ipts(0)) << pnew_tuple;
	      polyindexlist << new int(pnew_tuple->newindex);
	      polyptslist << new (const Array<Real> *)(pnew_tuple->pprojpt);
	    }
	  }
	}
      }
    }
    Array<int> polyindex;
    LinkedList_to_Array(&polyindex, polyindexlist);
    if (polyindex.get_size() > 0) {
      if (polyindex.get_size() == nax) {
	(*pcutpoly) << new Array<int>(polyindex);
	if (pcutcolors) {
	  colorslist << new Array<Real>(colors(icol % colors.get_size()));
	}
      }
      /*	else if (polyindex.get_size()==4 && nax==3) {
		for (int s=0; s<2; s++) {
		Array<int> *ppoly=new Array<int>(3);
		for (int i=0; i<3; i++) {
		(*ppoly)(i)=polyindex(s+i);
		}
		(*pcutpoly) << ppoly;
		if (pcutcolors) {
		colorslist << new Array<Real>(colors(icol % colors.get_size()));
		}
		}
		}*/
      else {
	LinkedList<Array<int> > subpolylist;
	Array<Array<Real> > subpts(polyindex.get_size());
	LinkedListIterator<const Array<Real> *> ipt(polyptslist);
	for (int j = 0; j < polyindex.get_size(); j++, ipt++) {
	  subpts(j) = **ipt;
	}
	create_mesh(&subpolylist, subpts, 0, MAXFLOAT);
	LinkedListIterator<Array<int> > isubpoly(subpolylist);
	for (; isubpoly; isubpoly++) {
	  Array<int> *pnewpoly = new Array<int>(isubpoly->get_size());
	  for (int j = 0; j < pnewpoly->get_size(); j++) {
	    (*pnewpoly)(j) = polyindex((*isubpoly)(j));
	  }
	  (*pcutpoly) << pnewpoly;
	  if (pcutcolors) {
	    colorslist << new Array<Real>(colors(icol % colors.get_size()));
	  }
	}
      }
    }
  }
  LinkedList_to_Array(pcutpts, cutptslist);
  if (pcutcolors) {
    LinkedList_to_Array(pcutcolors, colorslist);
  }
}

void simple_cross_section_simplexes(
    Array<Array<Real> > *pcutpts, LinkedList<Array<int> > *pcutpoly,
    Array<Array<Real> > *pcutcolors, const Array<int> &iaxes,
    const Array<Real> &org, const Array<Array<Real> > &pts,
    const LinkedList<Array<int> > &poly, const Array<Array<Real> > &colors) {
  int dim = pts(0).get_size();
  int nax = iaxes.get_size();
  int nlam = dim - nax + 1;
  int nnor = dim - nax;

  Array<int> inormal(nnor);
  int j = 0;
  for (int i = 0; i < inormal.get_size(); i++) {
    while (1) {
      int k;
      for (k = 0; k < nax; k++) {
	if (iaxes(k) == j) {
	  break;
	}
      }
      if (k == nax) break;
      j++;
    }
    inormal(i) = j;
  }
  Array<Real> shift(nnor);
  for (int i = 0; i < inormal.get_size(); i++) {
    shift(i) = org(inormal(i));
  }

  int curptsindex = 0;
  LinkedList<Array<Real> > cutptslist;
  LinkedList<Array<Real> > colorslist;
  Array<LinkedList<TupleDoneInfo> > tuple_done(pts.get_size());
  int icol = 0;
  LinkedListIterator<Array<int> > ipoly(poly);
  for (; ipoly; ipoly++, icol++) {
    LinkedList<int> polyindexlist;
    LinkedList<const Array<Real> *> polyptslist;
    NTupleIterator it(nlam, dim);
    // cerr << nlam << " " << dim << endl;
    for (; it; it++) {
      // cerr << "it=" << (Array<int> &)it << endl;
      Array<int> ipts(nlam);
      for (int i = 0; i < ipts.get_size(); i++) {
	ipts(i) = (*ipoly)(((Array<int> &)it)(i));
      }
      sort_array(&ipts, TrivLessThan<int>());
      LinkedListIterator<TupleDoneInfo> ituple(tuple_done(ipts(0)));
      for (; ituple; ituple++) {
	int j;
	for (j = 1; j < ipts.get_size(); j++) {
	  if (ipts(j) != ituple->orgindex(j)) break;
	}
	if (j == ipts.get_size()) break;
      }
      if (ituple) {
	polyindexlist << new int(ituple->newindex);
	polyptslist << new (const Array<Real> *)(ituple->pprojpt);
      } else {
	Array<Real> c(nlam);
	Array<int> signflip(nnor);
	zero_array(&signflip);
	Array<Array<Real> > proj(ipts.get_size());
	for (int i = 0; i < proj.get_size(); i++) {
	  proj(i).resize(nnor);
	  for (int j = 0; j < nnor; j++) {
	    proj(i)(j) = pts(ipts(i))(inormal(j));
	  }
	  diff(&proj(i), proj(i), shift);
	  for (int j = 0; j < signflip.get_size(); j++) {
	    if (sgn(proj(i)(j)) != sgn(proj(0)(j))) {
	      signflip(j) = 1;
	    }
	  }
	}
	if (max(signflip) > 0) {
	  Array2d<Real> A(nlam, nlam);
	  zero_array(&c);
	  c(0) = 1.;
	  for (int i = 0; i < nlam; i++) {
	    A(0, i) = 1.;
	    for (int j = 1; j < nlam; j++) {
	      A(j, i) = proj(i)(j - 1);
	    }
	  }
	  if (solve_linsys(&A, &c)) {
	    if (min(c) >= 0 && max(c) <= 1) {
	      Array<Real> newpt(dim);
	      zero_array(&newpt);
	      for (int i = 0; i < nlam; i++) {
		Array<Real> tmp;
		product(&tmp, pts(ipts(i)), c(i));
		sum(&newpt, newpt, tmp);
	      }

	      Array<Real> *pnewprojpt = new Array<Real>(nax);
	      {
		Array<Real> shnewpt;
		diff(&shnewpt, newpt, org);
		for (int j = 0; j < nax; j++) {
		  (*pnewprojpt)(j) = shnewpt(iaxes(j));
		}
	      }
	      cutptslist << pnewprojpt;
	      TupleDoneInfo *pnew_tuple =
		  new TupleDoneInfo(ipts, curptsindex, pnewprojpt);
	      curptsindex++;
	      tuple_done(ipts(0)) << pnew_tuple;
	      polyindexlist << new int(pnew_tuple->newindex);
	      polyptslist << new (const Array<Real> *)(pnew_tuple->pprojpt);
	    }
	  }
	}
      }
    }
    Array<int> polyindex;
    LinkedList_to_Array(&polyindex, polyindexlist);
    if (polyindex.get_size() > 0) {
      if (polyindex.get_size() == nax) {
	(*pcutpoly) << new Array<int>(polyindex);
	if (pcutcolors) {
	  colorslist << new Array<Real>(colors(icol % colors.get_size()));
	}
      }
      /*	else if (polyindex.get_size()==4 && nax==3) {
		for (int s=0; s<2; s++) {
		Array<int> *ppoly=new Array<int>(3);
		for (int i=0; i<3; i++) {
		(*ppoly)(i)=polyindex(s+i);
		}
		(*pcutpoly) << ppoly;
		if (pcutcolors) {
		colorslist << new Array<Real>(colors(icol % colors.get_size()));
		}
		}
		}*/
      else {
	LinkedList<Array<int> > subpolylist;
	Array<Array<Real> > subpts(polyindex.get_size());
	LinkedListIterator<const Array<Real> *> ipt(polyptslist);
	for (int j = 0; j < polyindex.get_size(); j++, ipt++) {
	  subpts(j) = **ipt;
	}
	// cerr << subpts << endl;
	// cerr << subpts.get_size() << " "; //DEBUG;
	create_mesh(&subpolylist, subpts, 0, MAXFLOAT);
	// cerr << subpolylist.get_size() << endl; //DEBUG;
	LinkedListIterator<Array<int> > isubpoly(subpolylist);
	for (; isubpoly; isubpoly++) {
	  Array<int> *pnewpoly = new Array<int>(isubpoly->get_size());
	  for (int j = 0; j < pnewpoly->get_size(); j++) {
	    (*pnewpoly)(j) = polyindex((*isubpoly)(j));
	  }
	  (*pcutpoly) << pnewpoly;
	  if (pcutcolors) {
	    colorslist << new Array<Real>(colors(icol % colors.get_size()));
	  }
	}
      }
    }
  }
  LinkedList_to_Array(pcutpts, cutptslist);
  if (pcutcolors) {
    LinkedList_to_Array(pcutcolors, colorslist);
  }
}

/*
void cross_section_simplexes(Array<Array<Real> > *pcutpts, const
Array<Array<Real> > &axes, const Array<Real> &org, const Array<Array<Real> >
&pts, const LinkedList<Array<int> > &poly) { int dim=pts(0).get_size(); int
nax=axes.get_size(); int nlam=dim-nax+1;

  Array2d<Real> normal(dim-nax,dim);
  Array<Real> shift(dim-nax);

  Array<Array<Real> > basis;
  for (int i=0; i<nax; i++) {
    build_ortho_basis(&basis,NULL,axes(i));
  }
  for (int i=0; i<dim-nax; ) {
    Array<Real> trialv(dim);
    for (int j=0; j<dim; j++) {trialv(j)=uniform01()-0.5;}
    if (build_ortho_basis(&basis,NULL,trialv)) {
      set_row(&normal,basis(basis.get_size()-1),i);
      i++;
    }
  }
  product(&shift,normal,org);

  LinkedList<Array<int> > tuple_list;
  LinkedListIterator<Array<int> > ipoly(poly);
  for ( ; ipoly; ipoly++) {
    Array<int> ipts(nlam);
    NTupleIterator it(nlam,dim);
    for ( ; it; it++) {
      for (int i=0; i<ipts.get_size(); i++) {
	ipts(i)=(*ipoly)(((Array<int> &)it)(i));
      }
      sort_array(&ipts,TrivLessThan<int>());
      LinkedListIterator<Array<int> > ituple(tuple_list);
      for ( ; ituple; ituple++) {
	if (ipts==*ituple) break;
      }
      if (!ituple) {
	tuple_list << new Array<int>(ipts);
      }
    }
  }

  //cerr << tuple_list << endl;

  LinkedList<Array<Real> > cutptslist;
  LinkedListIterator<Array<int> > ituple(tuple_list);
  for ( ; ituple; ituple++) {

    Array<Array<Real> > proj(ituple->get_size());
    for (int i=0; i<proj.get_size(); i++) {
      product(&proj(i),normal,pts((*ituple)(i)));
      diff(&proj(i),proj(i),shift);
    }

    Array2d<Real> A(nlam,nlam);
    Array<Real> c(nlam);
    zero_array(&c);
    c(0)=1.;
    for (int i=0; i<nlam; i++) {
      A(0,i)=1.;
      for (int j=1; j<nlam; j++) {
	A(j,i)=proj(i)(j-1);
      }
    }
    solve_linsys(&A,&c);
    if (min(c)>=0 && max(c)<=1) {
      Array<Real> *pnewpt=new Array<Real>(dim);
      zero_array(pnewpt);
      for (int i=0; i<nlam; i++) {
	Array<Real> tmp;
	product(&tmp,pts((*ituple)(i)),c(i));
	sum(pnewpt,*pnewpt,tmp);
      }
      cutptslist << pnewpt;
    }
  }

  pcutpts->resize(cutptslist.get_size());
  LinkedListIterator<Array<Real> > ipts(cutptslist);
  for (int i=0; i<pcutpts->get_size(); i++, ipts++) {
    (*pcutpts)(i).resize(nax);
    Array<Real> shpt;
    diff(&shpt,*ipts,org);
    for (int j=0; j<nax; j++) {
      (*pcutpts)(i)(j)=inner_product(shpt,axes(j));
    }
  }
}
*/

/*
void cross_section_simplex(Array<Array<Real> > *pcutpts, LinkedList<Array<int> >
*pcutpolylist, const Array2d<Real > &normal, const Array<Real> &shift, const
Array<Array<Real> > &pts, const Array<int> &poly) { int dim=pts(0).get_size();
  int nlam=normal.get_size()+1;

  Array<Array<Real> > proj(pts.get_size());
  for (int i=0; i<pts.get_size(); i++) {
    inner_product(&proj(i),normal,pts(poly(i)));
    diff(&proj(i),proj(i),shift);
  }

  LinkedList<Array<Real> > ptslist;

  Array2d<Real> A(nlam,nlam);
  Array<Real> c(nlam);
  NTupleIterator it(nlam,pts.get_size());
  for ( ; it; it++) {
    zero_array(&c);
    c(0)=1.;
    for (int i=0; i<nlam; i++) {
      A(0,i)=1.;
    }
    Array<int> &ita=(Array<int> &)it;
    for (int j=0; j<nlam; j++) {
      set_row(&A,proj(ita(j)),j+1);
    }
      solve_linsys(&A,&c);
    if (min(c)>=0 && max(c)<=1) {
      Array<Real> *pnewpt=new Array<Real>(dim);
      for (int i=0; i<nlam; i++) {
	Array<Real> tmp;
	product(&tmp,pts(ita(i)),c);
	sum(pnewpt,*pnewpt,tmp);
      }
      ptslist << pnewpts;
    }
  }
}
*/

