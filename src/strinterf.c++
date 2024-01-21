#include "strinterf.h"
#include "parse.h"

int num_strain_dim(int n) {
  n=n%4;
  return (n*(n+1)/2);
}

void Array6_to_Matrix3d(rMatrix3d *pm, const Array<Real> &x, int dimstrain) {
  pm->zero();
  int uaxis=dimstrain/4;
  dimstrain=dimstrain%4;
  iVector3d sh;
  if (dimstrain==3) {
    sh=iVector3d(0,3,5);
  }
  else if (dimstrain==2) {
    sh=iVector3d(0,2,0);
  }
  else {
    sh=iVector3d(0,0,0);
  }
  for (int i=0; i<dimstrain; i++) {
    for (int j=0; j<dimstrain; j++) {
      int ii=i;
      int jj=j;
      if (j<i) {swap(&ii,&jj);}
      (*pm)((i+uaxis)%3,(j+uaxis)%3)=x(sh(ii)+jj-ii);
    }
  }
}

/*
void Array6_to_Matrix3d(rMatrix3d *pm, const Array<Real> &x, int dimstrain) {
  if (dimstrain==3) {
    for (int i=0; i<3; i++) {
      for (int j=0; j<3; j++) {
	int ii=i;
	int jj=j;
	if (j<i) {swap(&ii,&jj);}
	(*pm)(i,j)=x((2*3+1-ii)*ii/2+jj-ii);
      }
    }
  }
  else {
    for (int i=0; i<3; i++) {
      for (int j=0; j<3; j++) {
	if (i<2 && j<2) {
	  int ii=i;
	  int jj=j;
	  if (j<i) {swap(&ii,&jj);}
	  (*pm)(i,j)=x(ii+jj);
	}
	else {
	  (*pm)(i,j)=0.;
	}
      }
    }
  }
}
*/

void Matrix3d_to_Array6(Array<Real> *px, const rMatrix3d &m, int dimstrain) {
  if (px->get_size()==0) {px->resize(num_strain_dim(dimstrain));}
  int uaxis=dimstrain/4;
  dimstrain=dimstrain%4;
  int k=0;
  for (int i=0; i<dimstrain; i++) {
    for (int j=i; j<dimstrain; j++) {
      (*px)(k)=m((i+uaxis)%3,(j+uaxis)%3);
      k++;
    }
  }
}

/*
void Matrix3d_to_Array6(Array<Real> *px, const rMatrix3d &m, int dimstrain) {
  if (dimstrain==3) {
    if (px->get_size()==0) {px->resize(6);}
    int k=0;
    for (int i=0; i<3; i++) {
      for (int j=i; j<3; j++) {
	(*px)(k)=m(i,j);
	k++;
      }
    }
  }
  else {
    if (px->get_size()==0) {px->resize(3);}
    int k=0;
    for (int i=0; i<2; i++) {
      for (int j=i; j<2; j++) {
	(*px)(k)=m(i,j);
	k++;
      }
    }
  }
}
*/

OptimizedStructure::OptimizedStructure(void): FunctionWithGrad(), grad(), val(0), axes(), str(), atom_label(), sumpos() {stressscale=0; strainscale=0; sleeptime=5; dimstrain=3;}

OptimizedStructure::~OptimizedStructure(void) {ofstream file("stop");}

void OptimizedStructure::init(Array<Real> *px, const rMatrix3d &_axes, const Structure &_str, const Array<AutoString> &_atom_label, Real forcefact, int _dimstrain) {
  axes=_axes;
  str=_str;
  dimstrain=_dimstrain;
  Real omega=det(str.cell)/(Real)(str.atom_pos.get_size());
  stressscale=pow(omega,2./3.)*forcefact;
  strainscale=pow(omega,1./3.)*(Real)(str.atom_pos.get_size())/forcefact;
  atom_label=_atom_label;
  px->resize(num_strain_dim(dimstrain)+3*(str.atom_pos.get_size()-1));
  rMatrix3d id;
  id.zero();//
  Matrix3d_to_Array6(px,id,dimstrain);
  sumpos=rVector3d(0.,0.,0.);
  for (int at=0; at<str.atom_pos.get_size(); at++) {
    rVector3d v=str.atom_pos(at);
    sumpos+=v;
    if (at>0) {
      for (int j=0; j<3; j++) {
	(*px)(num_strain_dim(dimstrain)+3*(at-1)+j)=v(j);
      }
    }
  }
}

void OptimizedStructure::vect_to_str(Structure *pstr, const Array<Real> &x) {
  pstr->atom_type=str.atom_type;
  rMatrix3d strainI;
  Array6_to_Matrix3d(&strainI,x,dimstrain);
  strainI=strainI/strainscale;
  rMatrix3d id;//
  id.identity();
  strainI=id+strainI;
  pstr->cell=strainI*str.cell;
  pstr->atom_pos.resize(str.atom_pos.get_size());
  rVector3d v0(0.,0.,0.);
  rMatrix3d icell=!str.cell;
  for (int at=1; at<pstr->atom_pos.get_size(); at++) {
    rVector3d v;
    for (int j=0; j<3; j++) {
      v(j)=x(num_strain_dim(dimstrain)+3*(at-1)+j);
    }
    pstr->atom_pos(at)=pstr->cell*icell*v;
    v0+=v;
  }
  pstr->atom_pos(0)=pstr->cell*icell*(sumpos-v0);
}

void OptimizedStructure::set_arg(const Array<Real> &x) {
  Real stressconv=1e8/1.6e-19*1e-30; // kB to eV/A^3;
  Structure curstr;
  cout << "x= ";
  tracev(x);
  vect_to_str(&curstr,x);
  {
    ofstream file("str.out");
    file.setf(ios::fixed);
    file.precision(5);
    write_structure(curstr,atom_label,axes,file,0);
  }
  {ofstream flag("busy");}
  while (file_exists("busy")) {
    sleep(sleeptime);//
  }
  {
    ifstream efile("energy");
    efile >> val;
  }
  cout << "read energy\n";
  grad.resize(num_strain_dim(dimstrain)+3*(str.atom_pos.get_size()-1));
  rMatrix3d strainI;
  {
    ifstream sfile("stress.out");
    rMatrix3d rawstress;
    rawstress(2,2)=MAXFLOAT;
    for (int i=0; i<3; i++) {
      for (int j=0; j<3; j++) {
	sfile >> rawstress(i,j);
      }
    }
    if (rawstress(2,2)==MAXFLOAT) {
      ERRORQUIT("Error reading stress.out file.");
    }
    /*
    for (int i=0; i<3; i++) {
      for (int j=0; j<3; j++) {
	if (i!=j) {rawstress(i,j)*=2.;}
      }
    }
    */
    cout << "read stress\n";
    rMatrix3d stress;
    rMatrix3d id;
    id.identity();//
    Array6_to_Matrix3d(&strainI,x,dimstrain);
    strainI=strainI/strainscale;
    strainI=id+strainI;
    stress=-stressscale*stressconv*rawstress*(~ ! strainI);
    //cout << "stress= " << stress << endl;
    for (int i=0; i<3; i++) {
      for (int j=i+1; j<3; j++) {
	stress(i,j)+=stress(j,i);
	stress(j,i)=0.;
      }
    }
    Matrix3d_to_Array6(&grad,stress,dimstrain);
  }
  {
    ifstream rfile("str_relax.out");
    Structure str_relax;
    rMatrix3d axes_relax;
    parse_structure_file(&(str_relax.cell), &(str_relax.atom_pos), &(str_relax.atom_type), atom_label, rfile, &axes_relax);
    Array<int> copyfrom;
    Array<iVector3d> cellshift;
    reorder_atoms(&str_relax,str,&copyfrom,&cellshift);
    cout << "read str_relax.out\n";
    
    ifstream ffile("force.out");
    Array<rVector3d> rawforce(str.atom_pos.get_size());
    rawforce(rawforce.get_size()-1)(2)=MAXFLOAT;
    for (int i=0; i<rawforce.get_size(); i++) {
      ffile >> rawforce(i);
    }
    cout << "read forces\n";
    if (rawforce(rawforce.get_size()-1)(2)==MAXFLOAT) {
      ERRORQUIT("Error reading force.out file.");
    }
    
    for (int at=1; at<rawforce.get_size(); at++) {
      rVector3d force=(~ strainI)*(rawforce(copyfrom(at))-rawforce(copyfrom(0)));
      for (int j=0; j<3; j++) {
	grad(num_strain_dim(dimstrain)+3*(at-1)+j)=-force(j);
      }
    }
  }
  cout << "done parsing output files\n";
}

void SmallEpicycle::init(const Array<Real> _c) {
  c=_c;
  pfunc->set_arg(c);
  pfunc->get_grad(&g0);
}

void SmallEpicycle::set_arg(const Array<Real> &_s) {
  cout << "Epicycle" << endl;
  Array<Real> x;
  s=_s;
  sum(&x,c,s);
  //    tracev(c);
  cout << "epidx= ";
  tracev(s);
  cout << endl;
  pfunc->set_arg(x);
  e=pfunc->get_val();
  Array<Real> g;
  pfunc->get_grad(&g);
  diff(&df,g,g0);
  product(&df,df,1./norm(s));
}
