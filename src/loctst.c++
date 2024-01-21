#include <fstream>
#include "parse.h"
#include "getvalue.h"
#include "version.h"

// write extra help as plain text in *.hlp
//extern char *helpstring;
char *helpstring="";

int main(int argc, char *argv[]) {
  // parsing command line. See getvalue.hh for details;
  char *wellstrfilename="wel_str.out";
  char *saddlestrfilename="sad_str.out";
  Real r=0;
  rVector3d center(0.,0.,0.);
  Array<Real> centera;
  FixedVector_to_Array(&centera,center);
  int sigdig=5;
  int dohelp=0;
  int dummy=0;
  AskStruct options[]={
    {"","Skeleton for atat-like codes " MAPS_VERSION ", by Axel van de Walle",TITLEVAL,NULL},
    {"-h","Display more help",BOOLVAL,&dohelp},
    {"-ws","Input file defining the well structure",STRINGVAL,&wellstrfilename},
    {"-ss","Input file defining the saddle structure",STRINGVAL,&saddlestrfilename},
    {"-c","Center",ARRAYRVAL,&centera},
    {"-r","Radius",REALVAL,&r},
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

  cout.setf(ios::fixed);
  cout.precision(sigdig);

  // parsing structure files. See parse.hh for detail;
  Structure str[2];
  Array<AutoString> label;
  rMatrix3d axes;
  {
    ifstream wellfile(wellstrfilename);
    if (!wellfile) ERRORQUIT("Unable to open well structure file");
    Array<Arrayint> labellookup;
    parse_lattice_file(&(str[0].cell), &(str[0].atom_pos), &(str[0].atom_type), &labellookup, &label, wellfile, &axes);
    wrap_inside_cell(&(str[0].atom_pos),str[0].atom_pos,str[0].cell);
    fix_atom_type(&str[0],labellookup);
  }

  {
    rMatrix3d axes1;
    ifstream saddlefile(saddlestrfilename);
    if (!saddlefile) ERRORQUIT("Unable to open saddle structure file.");
    parse_structure_file(&(str[1].cell), &(str[1].atom_pos), &(str[1].atom_type), label, saddlefile, &axes1);
    wrap_inside_cell(&(str[1].atom_pos),str[1].atom_pos,str[1].cell);
    reorder_atoms(&str[1],str[0]);
  }

  Array_to_FixedVector(&center,centera);
  center=axes*center;


  int fixsh=label.get_size();
  Array<AutoString> label_f(fixsh*2);
  for (int i=0; i<fixsh; i++) {
    label_f(i)=label(i);
    label_f(i)+=AutoString("_T");
    label_f(i+fixsh)=label(i);
    label_f(i+fixsh)+=AutoString("_F");
  }

  for (int ws=0; ws<2; ws++) {
    rMatrix3d icell=!(str[ws].cell);
    for (int at=0; at<str[ws].atom_pos.get_size(); at++) {
      Real d=norm(str[ws].cell*cylinder(icell*(str[ws].atom_pos(at)-center)));
      if (d<=r) {
	str[ws].atom_type(at)+=fixsh;
      }
    }
    write_structure(str[ws],label_f,axes,cout,0);
  }

}
