#include <fstream>
#include "parse.h"
#include "getvalue.h"
#include "version.h"

// write extra help as plain text
char *helpstring="More help!\n";

int main(int argc, char *argv[]) {
  // parsing command line. See getvalue.hh for details;
  char *strfilename="str.in";
  Real r=0;
  int n=0;
  int b=0;
  Array<Real> v;
  int sigdig=5;
  int dohelp=0;
  int dummy=0;
  AskStruct options[]={
    {"","Skeleton for atat-like codes " MAPS_VERSION ", by Axel van de Walle",TITLEVAL,NULL},
    {"-h","Display more help",BOOLVAL,&dohelp},
    {"-s","Input file defining a structure (Default: str.in)",STRINGVAL,&strfilename},
    {"-r","A real parameter",REALVAL,&r},
    {"-n","An integer parameter",INTVAL,&n},
    {"-b","A boolean parameter",BOOLVAL,&b},
    {"-v","A vector of reals",ARRAYRVAL,&v},
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

  //output what we have parsed
  cout.setf(ios::fixed);
  cout.precision(sigdig);
  cout << r << endl;
  cout << n << endl;
  cout << b << endl;
  cout << v << endl;

  // parsing structure file. See parse.hh for detail;
  Structure str; // contains cell and atom positions in cartesian and atom types as integers;
  Array<AutoString> label; //these integers will point to strings in this array;
  rMatrix3d axes; // coordinate system the user wants to use (multiply coord by (!axes) to get them into user coordinates);
  {
    Array<Arrayint> labellookup;
    ifstream strfile(strfilename);
    if (!strfile) ERRORQUIT("Unable to open structure file");
    parse_lattice_file(&str.cell, &str.atom_pos, &str.atom_type, &labellookup, &label, strfile, &axes);
    wrap_inside_cell(&str.atom_pos,str.atom_pos,str.cell);
    fix_atom_type(&str, labellookup);   // str.atom_type converted to indices into label;
  }
  write_structure(str,label,axes,cout);
}
