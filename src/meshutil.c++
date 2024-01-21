#include "meshutil.h"
#include "getvalue.h"

void MeshExporter::write(ostream &file, const Array<rVector3d> &pts, const LinkedList<iVector3d> &poly, const Array<Array<Real> > &ptsdata) {
  LinkedList<Array<int> > apoly;
  LinkedListFixedVector_to_LinkedListArray(&apoly,poly);
  write(file,pts,apoly,ptsdata);
}

void MeshExporter::write(ostream &file, const Array<Array<Real> > &pts, const LinkedList<Array<int> > &poly, const Array<Array<Real> > &ptsdata) {
  Array<rVector3d> fpts;
  ArrayArray_to_ArrayFixedVector(&fpts,pts);
  write(file,fpts,poly,ptsdata);
}

void BuiltInMeshExporter::write(ostream &file, const Array<rVector3d> &pts, const LinkedList<Array<int> > &poly, const Array<Array<Real> > &ptsdata) {
  Array<Array<Real> > apts;
  ArrayFixedVector_to_ArrayArray(&apts,pts);
  write(file,apts,poly,ptsdata);
}

void BuiltInMeshExporter::write(ostream &file, const Array<Array<Real> > &pts, const LinkedList<Array<int> > &poly, const Array<Array<Real> > &ptsdata) {
  file << pts;
  file << poly;
  file << ptsdata;
}

void VTKMeshExporter::write(ostream &file, const Array<rVector3d> &pts, const LinkedList<Array<int> > &poly, const Array<Array<Real> > &ptsdata) {
  if (pts.get_size()==0 || poly.get_size()==0) return;
  file << "# vtk DataFile Version 3.0" << endl;
  file << "vtk output" << endl;
  file << "ASCII" << endl;
  file << "DATASET POLYDATA" << endl;
  file << "POINTS " << pts.get_size() << " float" << endl;
  for (int p=0; p<pts.get_size(); p++) {
    file << pts(p) << endl;
  }
  LinkedListIterator<Array<int> > t(poly);
  int polytype=t->get_size();
  if (polytype==2) {
    file << "LINES ";
  }
  else {
    file << "POLYGONS ";
  }
  file  << poly.get_size() << " " << poly.get_size()*(polytype+1) <<endl;
  for (; t; t++) {
    file << t->get_size();
    for (int i=0; i<t->get_size(); i++) {
      file << " " << (*t)(i);
    }
    file << endl; 
  }
  if (ptsdata.get_size()>0) {
    int dimdata=ptsdata(0).get_size();
    if (dimdata==3 || dimdata==4) {
      file << "CELL_DATA " <<  poly.get_size() << endl;
      file << "SCALARS mycolor unsigned_char " << 4 << endl;
      file << "LOOKUP_TABLE default" << endl;
      for (int p=0; p<poly.get_size(); p++) {
	for (int j=0; j<4; j++) {
	  file << (j<dimdata ? (int)(ptsdata(p%ptsdata.get_size())(j)) : 255);
	  if (j<3) {file << " ";}
	}
	file << endl;
      }
    }
    else {
      file << "CELL_DATA " <<  poly.get_size() << endl;
      file << "SCALARS mycolor float " << dimdata << endl;
      file << "LOOKUP_TABLE default" << endl;
      for (int p=0; p<poly.get_size(); p++) {
	for (int j=0; j<dimdata; j++) {
	  file << ptsdata(p%ptsdata.get_size())(j);
	  if (j<dimdata-1) {file << " ";}
	}
	file << endl;
      }
    }
  }
}

template <> GenericPlugIn<MeshExporter> *GenericPlugIn<MeshExporter>::list=NULL;

SpecificPlugIn<MeshExporter,VTKMeshExporter> VTKPlugIn("vtk");
SpecificPlugIn<MeshExporter,BuiltInMeshExporter> BuiltInPlugIn("bi");

void combine_mesh(LinkedList<rVector3d> *ppts, LinkedList<Array<int> > *ppoly, const Array<rVector3d> &newpts, const LinkedList<Array<int> > &newpoly) {
  int shift=ppts->get_size();
  LinkedListIterator<Array<int> > it(newpoly);
  for (; it; it++) {
    Array<int> *p=new Array<int>(*it);
    (*ppoly) << p;
    for (int i=0; i<p->get_size(); i++) {
      (*p)(i)+=shift;
    }
  }
  (*ppts) << newpts;
}

void combine_mesh(LinkedList<Array<Real> > *ppts, LinkedList<Array<int> > *ppoly, const Array<Array<Real> > &newpts, const LinkedList<Array<int> > &newpoly) {
  int shift=ppts->get_size();
  LinkedListIterator<Array<int> > it(newpoly);
  for (; it; it++) {
    Array<int> *p=new Array<int>(*it);
    (*ppoly) << p;
    for (int i=0; i<p->get_size(); i++) {
      (*p)(i)+=shift;
    }
  }
  (*ppts) << newpts;
}

void transform(Array<rVector3d> *ptpts, const Array<rVector3d> pts, const rMatrix3d &op, const rVector3d &trans) {
  ptpts->resize(pts.get_size());
  for (int i=0; i<pts.get_size(); i++) {
    (*ptpts)(i)=op*pts(i)+trans;
  }
}

void PolyFont3D::init(ifstream &file) {
  AutoString dummy;
  get_string(&dummy,file,"\n");
  while (!file.eof()) {
    get_string(&dummy,file," ");
    if (file.eof()) break;
    skip_delim(file," ");
    char c;
    file.get(c);
    get_string(&dummy,file," ");
    file >> charwidth(c);
    get_string(&dummy,file," ");
    int numpts;
    file >> numpts;
    charpts(c).resize(numpts);
    for (int i=0; i<numpts; i++) {
      file >> charpts(c)(i);
    }
    get_string(&dummy,file," ");
    int numpoly;
    file >> numpoly;
    for (int i=0; i<numpoly; i++) {
      Array<int> *ppoly=new Array<int>();
      file >> (*ppoly);
      charpoly(c) << ppoly;
    }
    get_string(&dummy,file,"\n");
  }
  charwidth(' ')=0.5;
}

void PolyFont3D::write(LinkedList<rVector3d> *plistpts, LinkedList<Array<int> > *plistpoly, const rVector3d &where, const rVector3d &right, const rVector3d &up, char *string) {
  rMatrix3d rot;
  rot.set_column(0,right);
  rot.set_column(1,up);
  rot.set_column(2,right^up);
  rVector3d t(where);
  for (int i=0; i<strlen(string); i++) {
    char c=string[i];
    Array<rVector3d> tpts;
    transform(&tpts,charpts(c),rot,t);
    combine_mesh(plistpts,plistpoly,tpts,charpoly(c));
    t+=right*(charwidth(c)+spacing);
  }
}

Real PolyFont3D::get_length(char *string, int includelast) {
  Real len=0;
  for (int i=0; i<strlen(string); i++) {
    len+=charwidth(string[i]);
  }
  len+=spacing*(Real)(strlen(string)-1+includelast);
  return len;
}
