#include <fstream>
#include "getvalue.h"
#include "version.h"
#include "integer.h"
#include "plugin.h"

/*
template <class T>
ostream& operator<<(ostream& s, const LinkedList<T> &l) {
  s << l.get_size() << endl;
  LinkedListIterator<T> i(l);
  for (; i; i++) {
    s << *i << endl;
  }
  return s;
}
*/

template <class T>
T& first(LinkedList<T> &l) {
  LinkedListIterator<T> i(l);
  return *i;
}

/*
Real dirhausdorf(const Array<rVector3d> &a, const Array<rVector3d> &b) {
  Real mx=0.;
  for (int i=0; i<a.get_size(); i++) {
    Real mn=MAXFLOAT;
    for (int j=0; j<b.get_size(); j++) {
      mn=min(norm(a(i)-b(j)),mn);
    }
    mx=max(mx,mn);
  }
  return mx;
}

*/

Real norm2d(rVector3d x) {
  x(0)=0;
  return norm(x);
}

Real directed_hausdorf(const Array<int> &a, const Array<int> &b, const Array<rVector3d> &pts) {
  Real mx=0.;
  for (int i=0; i<a.get_size(); i++) {
    Real mn=MAXFLOAT;
    for (int j=0; j<b.get_size(); j++) {
      mn=min(norm2d(pts(a(i))-pts(b(j))),mn);
    }
    mx=max(mx,mn);
  }
  return mx;
}

Real hausdorf(const Array<int> &a, const Array<int> &b, const Array<rVector3d> &pts) {
  return max(directed_hausdorf(a,b,pts),directed_hausdorf(b,a,pts));
}

Real directed_hausdorf(int a, const Array<int> &b, const Array<rVector3d> &pts) {
  Real mn=MAXFLOAT;
  for (int j=0; j<b.get_size(); j++) {
    mn=min(norm2d(pts(a)-pts(b(j))),mn);
  }
  return mn;
}

template <class T>
void concat(Array<T> *pc, const Array<T> &a, const Array<T> &b) {
  pc->resize(a.get_size()+b.get_size());
  int j=0;
  for (int i=0; i<a.get_size(); i++,j++) {
    (*pc)(j)=a(i);
  }
  for (int i=0; i<b.get_size(); i++,j++) {
    (*pc)(j)=b(i);
  }
}

/*
int flip_if_needed(Array<rVector3d> *pa, const Array<rVector3d> &b) {
  Real nf=norm((*pa)(0)-b(0))+norm((*pa)(pa->get_size()-1)-b(b.get_size()-1));
  Real nr=norm((*pa)(0)-b(b.get_size()-1))+norm((*pa)(pa->get_size()-1)-b(0));
  if (nr<nf) {
    for (int i=0; i<pa->get_size()/2; i++) {
      swap(&((*pa)(i)),&((*pa)(pa->get_size()-1-i)));
    }
  }
}
*/

template <class T>
class Permutor {
  Array<T> org;
  Array<T> cur;
  MultiDimIterator<Array<int> > it;
public:
  Permutor(const Array<T> &_org): org(_org), cur(), it() {
    Array<int> m(org.get_size());
    for (int j=0; j<m.get_size(); j++) {m(j)=m.get_size()-j;}
    it.init(m);
  }
  operator void * () {return (void *)it;}
  operator const Array<T>& (void) {
    if (cur.get_size()==0) {
      cur.resize(org.get_size());
      Array<int> done(org.get_size());
      zero_array(&done);
      const Array<int> &p=it;
      for (int i=0; i<cur.get_size(); i++) {
	int j=0;
	int k=0;
	while (1) {
	  while (done(k)==1) {k++;}
	  if (j==p(i)) break;
	  k++;
	  j++;
	}
	cur(k)=org(i);
	done(k)=1;
      }
    }
    return cur;
  }
  void operator++ (int) {
    it++;
    cur.resize(0);
  }
};

template <class T>
class PermutorUnique {
  Permutor<T> p;
  LinkedList<Array<T> > list;
public:
  PermutorUnique(const Array<T> &_org): p(_org), list() {list << new Array<T>(p);}
  operator void * () {return (void *)p;}
  operator const Array<T>& (void) {
    return p;
  }
  void operator++ (int) {
    while (1) {
      p++;
      if (!p) break;
      const Array<T> &a=p;
      LinkedListIterator<Array<T> > i(list);
      for (; i; i++) {
	if (*i==a) break;
      }
      if (!i) {
	list << new Array<T>(a);
	break;
      }
    }
  }
};

void make_triangles(LinkedList<iVector3d> *ptriangles, const Array<int> &a0, const Array<int> &a1, int i1org, int i1end, const Array<rVector3d> &pts) {
  int di1=1;
  if (norm2d(pts(a1(i1org))-pts(a0(0)))>norm2d(pts(a1(i1end))-pts(a0(0)))) {
    di1=-1;
  }
  int i0=0;
  int i1=(di1==1 ? i1org : i1end);
  while (i0<a0.get_size() && i1>=i1org && i1<=i1end) {
    int inci0=1;
    int inci1=1;
    if (i0+1>=a0.get_size()) {inci0=0;}
    if (i1+di1<i1org || i1+di1>i1end) {inci1=0;}
    if (!inci0 && !inci1) break;
    if (inci0 && inci1) {
      rVector3d dx0=pts(a0(i0+1))-pts(a0(i0));
      rVector3d dx1=pts(a1(i1+di1))-pts(a1(i1));
      Real p0=(pts(a1(i1))-pts(a0(i0)))*dx0/norm2(dx0);
      Real p1=(pts(a0(i0))-pts(a1(i1)))*dx1/norm2(dx1);
      if (p0>=0 && p0<=1) {
	inci1=0;
      }
      else if (p1>=0 && p1<=1) {
	inci0=0;
      }
      else if (p0<p1) {
	inci0=0;
      }
      else {
	inci1=0;
      }
    }
    iVector3d *ptri=new iVector3d;
    if (inci0) {
      (*ptri)(0)=a0(i0);
      (*ptri)(1)=a1(i1);
      (*ptri)(2)=a0(i0+1);
      i0++;
    }
    else {
      (*ptri)(0)=a1(i1);
      (*ptri)(1)=a1(i1+di1);
      (*ptri)(2)=a0(i0);
      i1+=di1;
    }
    (*ptriangles) << ptri;
  }
}

void fix_triangle_normals_sub(int t, const rVector3d &dir, Array<int> *ptoflip, const Array<LinkedList<int> > &pts2tri, const Array<iVector3d> &triangles, const Array<rVector3d> &pts) {
  rVector3d u=(pts(triangles(t)(1))-pts(triangles(t)(0))) ^ (pts(triangles(t)(2))-pts(triangles(t)(1)));
  if (dir*u>=0) {(*ptoflip)(t)=1;} else {(*ptoflip)(t)=-1; u=-u;}
  for (int v=0; v<3; v++) {
    LinkedListIterator<int> i(pts2tri(triangles(t)(v)));
    for (; i; i++) {
      if ((*i)!=t) {
	int j;
	for (j=0; j<3; j++) {
	  if (triangles(*i)(j)==triangles(t)((v+1)%3)) break;
	}
	if (j<3) break;
      }
    }
    if (i) {
      if ((*ptoflip)(*i)==0) {
	fix_triangle_normals_sub(*i,u,ptoflip,pts2tri,triangles,pts);
      }
    }
  }
}

void fix_triangle_normals(LinkedList<iVector3d> *ptriangles, const Array<rVector3d> &pts) {
  Array<iVector3d> triangles;
  LinkedList_to_Array(&triangles,*ptriangles);
  Array<LinkedList<int> > pts2tri(pts.get_size());
  for (int t=0; t<triangles.get_size(); t++) {
    for (int i=0; i<3; i++) {
      pts2tri(triangles(t)(i)) << new int(t);
    }
  }
  Array<int> toflip(triangles.get_size());
  zero_array(&toflip);
  while (1) {
    int f;
    for (f=0; f<toflip.get_size(); f++) {
      if (toflip(f)==0) break;
    }
    if (f==toflip.get_size()) break;
    fix_triangle_normals_sub(f,rVector3d(0.,0.,0.),&toflip,pts2tri,triangles,pts);
  }
  LinkedListIterator<iVector3d> t(*ptriangles);
  for (int it=0; t; t++, it++) {
    if (toflip(it)==-1) {
      swap(&((*t)(1)),&((*t)(2)));
    }
  }
}

class TriangulationExporter {
public:
  virtual void write_onecolor(ostream &file, const Array<rVector3d> &pts, const LinkedList<iVector3d> &tri, int color) {}
};

template <> GenericPlugIn<TriangulationExporter> *GenericPlugIn<TriangulationExporter>::list=NULL;

class VTKTriangulationExporter: public TriangulationExporter {
  virtual void write_onecolor(ostream &file, const Array<rVector3d> &pts, const LinkedList<iVector3d> &tri, int color) {
    cout << "# vtk DataFile Version 3.0" << endl;
    cout << "vtk output" << endl;
    cout << "ASCII" << endl;
    cout << "DATASET POLYDATA" << endl;
    cout << "POINTS " << pts.get_size() << " float" << endl;
    for (int p=0; p<pts.get_size(); p++) {
      cout << pts(p) << endl;
    }
    cout << "POLYGONS " << tri.get_size() << " " << tri.get_size()*4 <<endl;
    LinkedListIterator<iVector3d> t(tri);
    for (; t; t++) {
      cout << 3 << " " << *t << endl; 
    }
    
    cout << "POINT_DATA " <<  pts.get_size() << endl;
    cout << "SCALARS mycolor int 1" << endl;
    cout << "LOOKUP_TABLE default" << endl;
    for (int p=0; p<pts.get_size(); p++) {
      cout << color << endl;
    }
  }
};

SpecificPlugIn<TriangulationExporter,VTKTriangulationExporter> VTKPlugIn("vtk");

const char *helpstring="Insert more help here";

int main(int argc, char *argv[]) {
  // parsing command line. See getvalue.hh for details;
  const char *tabfilename="tab.in";
  //  Real maxd=0.;
  Real sc=1000.;
  int color=0;
  int sigdig=5;
  const char *exporter_label="vtk";
  int nofixnormal=0;
  int backface=0;
  int dohelp=0;
  int dummy=0;
  AskStruct options[]={
    {"","TRIangulate PHase " MAPS_VERSION ", by Axel van de Walle",TITLEVAL,NULL},
    {"-t","Input table file",STRINGVAL,&tabfilename},
    //    {"-md","Max distance between points of the same surface",REALVAL,&maxd},
    {"-sc","Scale for the unique axis (1st)",REALVAL,&sc},
    {"-col","Color",INTVAL,&color},
    {"-nfn","do Not fix Normal vector",BOOLVAL,&nofixnormal},
    {"-bf","generate back face",BOOLVAL,&backface},
    {"-sig","Number of significant digits printed (default: 5)",INTVAL,&sigdig},
    {"-ef","Export format (default: vtk)",STRINGVAL,&exporter_label},
    {"-h","Display more help",BOOLVAL,&dohelp},
    {"-d","Use all default values",BOOLVAL,&dummy}
  };
  if (!get_values(argc,argv,countof(options),options)) {
    display_help(countof(options),options);
    return 1;
  }
  if (dohelp) {
    cout << helpstring;
    return 1;
  }

  cout.setf(ios::fixed);
  cout.precision(sigdig);

  if (!check_plug_in(TriangulationExporter(),exporter_label)) {
    cerr << exporter_label << endl;
    ERRORQUIT("Unknown export format");
  }
  TriangulationExporter *pexporter=GenericPlugIn<TriangulationExporter>::create(exporter_label);

  // this needs to be cleaned up;
  Array<rVector3d> pts;
  LinkedList<Array<int> > strips;
  {
    Array<Array<Real> > rawtab;
    ifstream tabfile(tabfilename);
    if (!tabfile) ERRORQUIT("Unable to open table file");
    read_table(&rawtab,tabfile,1);
    LinkedList<rVector3d> ptslist;
    LinkedList<int> curstrip;
    int j=0;
    for (int i=0; i<=rawtab.get_size(); i++) {
      if (i==rawtab.get_size() || rawtab(i).get_size()==0) {
	if (curstrip.get_size()>0) {
	  Array<int> *pstrip=new Array<int>;
	  LinkedList_to_Array(pstrip,curstrip);
	  strips << pstrip;
	  curstrip.delete_all();
	}
      }
      else {
	rVector3d *pv=new rVector3d(rawtab(i)(0)/sc,rawtab(i)(1)+rawtab(i)(2)/2.,rawtab(i)(2)*sqrt(3.)/2);
	ptslist << pv;
	curstrip << new int(j);
	j++;
      }
    }
    LinkedList_to_Array(&pts, ptslist);
  }
  LinkedList<iVector3d> triangles;

  Real cursl0=MAXFLOAT;
  for (int i=1; i<pts.get_size(); i++) {cursl0=min(pts(i)(0),cursl0);}
  while (1) {
    Real cursl1=MAXFLOAT;
    for (int i=1; i<pts.get_size(); i++) {
      if (pts(i)(0)>cursl0+zero_tolerance) {cursl1=min(cursl1,pts(i)(0));}
    }
    if (cursl1==MAXFLOAT) break;
    
    int nsl0=0;
    int nsl1=0;
    LinkedListIterator<Array<int> > s(strips);
    for (; s; s++) {
      if (near_zero(pts((*s)(0))(0)-cursl0)) nsl0++;
      if (near_zero(pts((*s)(0))(0)-cursl1)) nsl1++;
    }
    Array<Array<int> *> pa0(nsl0);
    Array<Array<int> *> pa1(nsl1);
    int i0=0;
    int i1=0;
    for (s.init(strips); s; s++) {
      if (near_zero(pts((*s)(0))(0)-cursl0)) {pa0(i0)=s; i0++;}
    }
    for (s.init(strips); s; s++) {
      if (near_zero(pts((*s)(0))(0)-cursl1)) {pa1(i1)=s; i1++;}
    }
    if (pa0.get_size()<pa1.get_size()) {
      swap(&pa0,&pa1);
    }
    Real bestd=MAXFLOAT;
    Array<int> bestmrg;
    Array<Array<int>* > bestperm;
    Array<int> merge(pa0.get_size());
    zero_array(&merge);
    for (int i=0; i<pa0.get_size()-pa1.get_size();i++) {merge(i)=1;}
    PermutorUnique<int> mrg(merge);
    for (; mrg; mrg++) {
      const Array<int> &amrg=mrg;
      Permutor<Array<int> * > perm(pa0);
      for (; perm; perm++) {
	const Array<Array<int> *> &aperm=perm;
	Real d=0;
	int s0=0;
	for (int s1=0; s1<pa1.get_size(); s1++, s0++) {
	  if (amrg(s1)==1) {
	    Array<int> comb;
	    concat(&comb,*aperm(s0),*aperm(s0+1));
	    d+=hausdorf(comb,*pa1(s1),pts);
	    s0++;
	  }
	  else {
	    d+=hausdorf(*aperm(s0),*pa1(s1),pts);
	  }
	}
	if (d<bestd) {
	  bestmrg=amrg;
	  bestperm=aperm;
	  bestd=d;
	}
      }
    }
    int s0=0;
    for (int s1=0; s1<pa1.get_size(); s1++, s0++) {
      if (bestmrg(s1)==1) {
	Array<Real> diffd(pa1(s1)->get_size());
	for (int i1=0; i1<pa1(s1)->get_size(); i1++) {
	  diffd(i1)=directed_hausdorf((*pa1(s1))(i1),*bestperm(s0),pts)-directed_hausdorf((*pa1(s1))(i1),*bestperm(s0+1),pts);
	}
	int i1split=-1;
	int hlf=diffd.get_size()/2;
	for (int di=0; di<hlf; di++) {
	  for (int sgn=-1; sgn<=1; sgn+=2) {
	    int i2=hlf+sgn*di-(sgn==-1);
	    if (i2+1<diffd.get_size()) {
	      if ((diffd(i2)>0)!=(diffd(i2+1)>0)) {
		i1split=i2;
		break;
	      }
	      //	      if (near_zero(diffd(i2)) || near_zero(diffd(i2+1))) {
	      //		i1split=i2;
	      //		break;
	      //	      }
	    }
	  }
	  if (i1split!=-1) break;
	}
	if (i1split==-1) {
	  make_triangles(&triangles,*bestperm(s0+(diffd(0)>0)),*pa1(s1),0,pa1(s1)->get_size()-1,pts);
	}
	else {
	  int shft=(fabs(diffd(i1split+1))<fabs(diffd(i1split)));
	  make_triangles(&triangles,*bestperm(s0+(diffd(i1split)>0)),*pa1(s1),0,i1split+shft,pts);
	  make_triangles(&triangles,*bestperm(s0+(diffd(i1split+1)>0)),*pa1(s1),i1split+shft,pa1(s1)->get_size()-1,pts);
	}
	s0++;
      }
      else {
	make_triangles(&triangles,*bestperm(s0),*pa1(s1),0,pa1(s1)->get_size()-1,pts);
      }
    }
    cursl0=cursl1;
  }

  if (!nofixnormal) {fix_triangle_normals(&triangles,pts);}
  if (backface) {
    LinkedListIterator<iVector3d> t(triangles);
    for (; t; t++) {
      triangles.add(new iVector3d((*t)(0),(*t)(2),(*t)(1)),t);
      t++;
    }
  }

  pexporter->write_onecolor(cout,pts,triangles,color);
  delete pexporter;

}
