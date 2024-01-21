#include <fstream>
#include "parse.h"
#include "getvalue.h"
#include "meshutil.h"
#include "version.h"
#include "meshcalc.h"

void calc_outshift(Array<rVector3d> *pshift, const Array<rVector3d> &vertex) {
    rVector3d middle(0.,0.,0.);
    for (int i=0; i<vertex.get_size(); i++) {
      middle+=vertex(i);
    }
    middle=middle*(1./(Real)(vertex.get_size()));
    pshift->resize(vertex.get_size());
    for (int i=0; i<vertex.get_size(); i++) {
      (*pshift)(i)=vertex(i)-middle;
      (*pshift)(i).normalize();
    }
}

void make_simplex(Array<Array<Real> > *psimplex, int dim) {
  psimplex->resize(dim+1);
  Array<Real> accum(dim);
  zero_array(&accum);
  (*psimplex)(0)=accum;
  for (int i=1; i<dim+1; i++) {
    Array<Real> &mid=(*psimplex)(i);
    mid=accum;
    product(&mid,mid,1./(Real)(i));
    Real d=sqrt(1.-norm2(mid));
    mid(i-1)+=d;
    sum(&accum,accum,mid);
  }
}

// write extra help as plain text in *.hlp
// extern const char *helpstring;
const char *helpstring="";

int main(int argc, char *argv[]) {
  // parsing command line. See getvalue.hh for details;
  int dotetra=0;
  int dotri=0;
  int hidim=3;
  int dotrans=0;
  Real color=0;
  Array<Real> rgb;
  const char *allvertexlabels="";
  const char *font_file="";
  const char *exporter_label="vtk";
  int docut=0;
  const char *cutfilename="cut.in";
  Real zmax=1000.;
  Real zmin=0.;
  Real zscale=1000.;
  int sigdig=5;
  int dohelp=0;
  int dummy=0;
  AskStruct options[]={
    {"","MaKe Axes " MAPS_VERSION ", by Axel van de Walle",TITLEVAL,NULL},
    {"-h","Display more help",BOOLVAL,&dohelp},
    {"-trans","Transform to coordinate system",BOOLVAL,&dotrans},
    {"-tri","Gibbs triangle",BOOLVAL,&dotri},
    {"-tetra","Gibbs tetrahedron",BOOLVAL,&dotetra},
    {"-hidim","Dimension",INTVAL,&hidim},
    {"-col","Color",REALVAL,&color},
    {"-rgb","RGB Color",ARRAYRVAL,&rgb},
    {"-vl","Vertices labels",STRINGVAL,&allvertexlabels},
    {"-zmax","",REALVAL,&zmax},
    {"-zmin","",REALVAL,&zmin},
    {"-zscale","",REALVAL,&zscale},
    {"-cs","Cross-Section",BOOLVAL,&docut},
    {"-cf","Cross section input File",STRINGVAL,&cutfilename},
    {"-font","font file",STRINGVAL,&font_file},
    {"-ef","Export format (default: vtk)",STRINGVAL,&exporter_label},
    {"-sig","Number of significant digits printed (Default: 5)",INTVAL,&sigdig},
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

  if (!check_plug_in(MeshExporter(),exporter_label)) {
    cerr << exporter_label << endl;
    ERRORQUIT("Unknown export format");
  }

  if (dotrans) {
    Array<Array<Real> > rawpts;
    read_table(&rawpts,cin,1);
    Array2d<Real> tmat;
    Array<Real> shft;
    if (hidim<3) {
      ERRORQUIT("Need dim>=3");
    }
    else if (hidim==3) {
      if (dotri) {
	tmat.resize(3,4);
	shft.resize(4);
	zero_array(&tmat);
	zero_array(&shft);
	tmat(0,2)=1.;
	tmat(0,3)=0.5;
	tmat(1,3)=sqrt(3.)/2.;
	tmat(2,0)=1./zscale;
	shft(0)=-zmin;
      }
      else if (dotetra) {
	tmat.resize(3,5);
	shft.resize(5);
	zero_array(&tmat);
	zero_array(&shft);
	tmat(0,2)=1.;
	tmat(0,3)=0.5;
	tmat(1,3)=sqrt(3.)/2.;
	tmat(0,4)=0.5;
	tmat(1,4)=1./sqrt(12.);
	tmat(2,4)=sqrt(2./3.);
      }
      else {
	ERRORQUIT("Must specify -tri or -tetra");
      }
    }
    else {
      if (dotri) {
	tmat.resize(hidim,hidim+1);
	shft.resize(hidim+1);
	zero_array(&tmat);
	zero_array(&shft);
	Array<Array<Real> > axes;
	make_simplex(&axes,hidim-1);
	for (int i=1; i<axes.get_size(); i++) {
	  for (int j=0; j<axes(i).get_size(); j++) {
	    tmat(j,1+i)=axes(i)(j);
	  }
	}
	tmat(hidim-1,0)=1./zscale;
	shft(0)=-zmin;
      }
      else if (dotetra) {
	tmat.resize(hidim,hidim+2);
	shft.resize(hidim+2);
	zero_array(&tmat);
	zero_array(&shft);
	Array<Array<Real> > axes;
	make_simplex(&axes,hidim);
	for (int i=1; i<axes.get_size(); i++) {
	  for (int j=0; j<axes(i).get_size(); j++) {
	    tmat(j,1+i)=axes(i)(j);
	  }
	}
      }
      else {
	ERRORQUIT("Must specify -tri (i.e. include temperature axis) or -tetra (no temperature axis)");
      }
    }
    for (int l=0; l<rawpts.get_size(); l++) {
      if (rawpts(l).get_size()==0) {
	cout << endl;
      }
      else {
	Array<Real> tpt,tmp;
	sum(&tmp,rawpts(l),shft);
	product(&tpt,tmat,tmp);
	for (int c=0; c<tpt.get_size(); c++) {
	  if (c>0) {cout << " ";}
	  cout << tpt(c);
	}
	cout << endl;
      }
    }
    exit(0);
  }

  MeshExporter *pexporter=GenericPlugIn<MeshExporter>::create(exporter_label);

  cout.setf(ios::fixed);
  cout.precision(sigdig);

  PolyFont3D font;
  Array<AutoString> labels(4);
  Real fontsize=0.05;
  if (strlen(allvertexlabels)>0) {
    AutoString aallvertexlabels(allvertexlabels);
    char *begstr=aallvertexlabels;
    char *endstr=begstr+strlen(begstr);
    for (int l=0; l<labels.get_size(); l++) {
      char *seek=begstr;
      while (seek<endstr && *seek!=',') {seek++;}
      *seek=0;
      labels(l).set(begstr);
      begstr=seek;
      if (begstr<endstr) {begstr++;}
    }

    AutoString fontfilename(font_file);
    if (strlen(font_file)==0) {
      get_atat_root(&fontfilename);
      fontfilename+="/data/font.dat";
    }
    ifstream fontfile(fontfilename);
    if (!fontfile) ERRORQUIT("Unable to open font file.");
    font.init(fontfile);
  }
  Array<Array<Real> > colors(1);
  if (rgb.get_size()>0) {
    colors(0)=rgb;
  }
  else {
    colors(0).resize(1);
    colors(0)(0)=color;
  }

  if (hidim==3) {
    if (dotri) {
      Array<rVector3d> vertex(3);
      vertex(0)=rVector3d(0.,0.,0.);
      vertex(1)=rVector3d(1.,0.,0.);
      vertex(2)=rVector3d(0.5,sqrt(3.)/2.,0.);
      if (strlen(allvertexlabels)>0) {
	Array<rVector3d> outshift;
	calc_outshift(&outshift,vertex);
	LinkedList<rVector3d> listpts;
	LinkedList<Array<int> > listpoly;
	rVector3d right=rVector3d(1.,0.,0.)*fontsize;
	rVector3d up=rVector3d(0.,0.,1.)*fontsize;
	for (int i=0; i<3; i++) {
	  rVector3d where=vertex(i)+rVector3d(0.,0.,(zmax-zmin)/zscale/2.0)+outshift(i)*fontsize*2.-right*font.get_length(labels(i))/2.;
	  font.write(&listpts,&listpoly,where,right,up,labels(i));
	}
	Array<rVector3d> apts;
	LinkedList_to_Array(&apts,listpts);
	pexporter->write(cout,apts,listpoly,colors);
      }
      else {
	Array<rVector3d> pts(6);
	for (int i=0; i<3; i++) {
	  pts(i)=vertex(i);
	  pts(3+i)=vertex(i)+rVector3d(0.,0.,(zmax-zmin)/zscale);
	}
	LinkedList<Array<int> > polys;
	for (int i=0; i<3; i++) {
	  Array<int> *ppoly=new Array<int>(2);
	  (*ppoly)(0)=i;
	  (*ppoly)(1)=(i+1)%3;
	  polys << ppoly;
	  ppoly=new Array<int>(2);
	  (*ppoly)(0)=i+3;
	  (*ppoly)(1)=((i+1)%3)+3;
	  polys << ppoly;
	  ppoly=new Array<int>(2);
	  (*ppoly)(0)=i;
	  (*ppoly)(1)=i+3;
	  polys << ppoly;
	}    
	pexporter->write(cout,pts,polys,colors);
      }
    }
    if (dotetra) {
      Array<rVector3d> vertex(4);
      vertex(0)=rVector3d(0.,0.,0.);
      vertex(1)=rVector3d(1.,0.,0.);
      vertex(2)=rVector3d(0.5,sqrt(3.)/2.,0.);
      vertex(3)=rVector3d(0.5,1./sqrt(12.),sqrt(2./3.));
      if (strlen(allvertexlabels)>0) {
	Array<rVector3d> outshift;
	calc_outshift(&outshift,vertex);
	LinkedList<rVector3d> listpts;
	LinkedList<Array<int> > listpoly;
	rVector3d right=rVector3d(1.,0.,0.)*fontsize;
	rVector3d up=rVector3d(0.,0.,1.)*fontsize;
	for (int i=0; i<4; i++) {
	  rVector3d where=vertex(i)+outshift(i)*fontsize*2.-right*font.get_length(labels(i))/2.;
	  font.write(&listpts,&listpoly,where,right,up,labels(i));
	}
	Array<rVector3d> apts;
	LinkedList_to_Array(&apts,listpts);
	pexporter->write(cout,apts,listpoly,colors);
      }
      else {
	LinkedList<Array<int> > polys;
	for (int i=1; i<4; i++) {
	  for (int j=0; j<i; j++) {
	    Array<int> *ppoly=new Array<int>(2);
	    (*ppoly)(0)=i;
	    (*ppoly)(1)=j;
	    polys << ppoly;
	  }
	}
	pexporter->write(cout,vertex,polys,colors);
      }
    }
  }
  else {
    Array<Array<Real> > pts;
    LinkedList<Array<int> > polys;
    if (dotri) {
      Array<Array<Real> > axes;
      make_simplex(&axes,hidim-1);
      pts.resize(hidim*2);
      for (int i=0; i<axes.get_size(); i++) {
	pts(i).resize(hidim);
	extract_elements(&(pts(i)),0,axes(i),0,hidim-1);
	pts(i)(hidim-1)=0.;
	pts(i+axes.get_size())=pts(i);
	pts(i+axes.get_size())(hidim-1)=(zmax-zmin)/zscale;
      }
      for (int omit=0; omit<hidim; omit++) {
	Array<int> sideface(2*(hidim-1));
	for (int i=0; i<hidim-1; i++) {
	  sideface(i)=i+(i>=omit ? 1 : 0);
	  sideface(i+hidim-1)=i+(i>=omit ? 1 : 0)+hidim;
	}
	for (int subsim=0; subsim<(hidim-1); subsim++) {
	  Array<int> *ppoly=new Array<int>(hidim);
	  for (int i=0; i<hidim; i++) {
	    (*ppoly)(i)=sideface(subsim+i);
	  }
	  polys << ppoly;
	}
      }
      for (int f=0; f<2; f++) {
	Array<int> *ppoly=new Array<int>(hidim);
	for (int i=0; i<hidim; i++) {
	  (*ppoly)(i)=i+f*hidim;
	}
	polys << ppoly;
      }
    }
    else if (dotetra) {
      make_simplex(&pts,hidim);
      for (int omit=0; omit<hidim+1; omit++) {
	Array<int> *ppoly=new Array<int>(hidim);
	for (int i=0; i<hidim; i++) {
	  (*ppoly)(i)=i+(i>=omit ? 1 : 0);
	}
	polys << ppoly;
      }
    }
    if (docut) {
      Array<Array<Real> > cutpts;
      LinkedList<Array<int> > cutpolys;

      Array<Real> org;
      Array<Array<Real> > axes;
      {
	ifstream cutfile(cutfilename);
	if (!cutfile) ERRORQUIT("Unable to open input file for cross section.");
	Array<Array<Real> > rawcut;
	read_table(&rawcut,cutfile,0);
	org=rawcut(0);
	axes.resize(rawcut.get_size()-1);
	for (int j=0; j<axes.get_size(); j++) {
	  diff(&(axes(j)),rawcut(j+1),org);
	}
      }
      cross_section_simplexes(&cutpts,&cutpolys, axes,org,pts,polys);
      pexporter->write(cout,cutpts,cutpolys,colors);
    }
    else {
      pexporter->write(cout,pts,polys,colors);
    }
  }
}
