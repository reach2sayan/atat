#include <fstream>
#include "parse.h"
#include "getvalue.h"
#include "version.h"

// write extra help as plain text in *help.hlp
extern char *helpstring;

int main(int argc, char *argv[]) {
  // parsing command line. See getvalue.hh for details;
  char *latfilename="lat.in";
  char *strfilename="str.out";
  Real r=0;
  int n=0;
  int b=0;
  int sigdig=5;
  int dohelp=0;
  int dummy=0;
  AskStruct options[]={
    {"","Skeleton for atat-like codes " MAPS_VERSION ", by Axel van de Walle",TITLEVAL,NULL},
    {"-h","Display more help",BOOLVAL,&dohelp},
    {"-l","Input file defining the lattice (Default: lat.in)",STRINGVAL,&latfilename},
    {"-s","Input file defining the structure (Default: str.out)",STRINGVAL,&strfilename},
    {"-r","A real parameter",REALVAL,&r},
    {"-n","An integer parameter",INTVAL,&n},
    {"-b","A boolean parameter",BOOLVAL,&b},
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
  cout << r << endl;
  cout << n << endl;
  cout << b << endl;

  // parsing lattice and structure files. See parse.hh for detail;
  Structure lat;
  Array<Arrayint> labellookup;
  Array<AutoString> label;
  rMatrix3d axes;
  {
    ifstream latfile(latfilename);
    if (!latfile) ERRORQUIT("Unable to open lattice file");
    parse_lattice_file(&lat.cell, &lat.atom_pos, &lat.atom_type, &labellookup, &label, latfile, &axes);
    wrap_inside_cell(&lat.atom_pos,lat.atom_pos,lat.cell);
  }

  Structure str;
  {
    ifstream strfile(strfilename);
    if (!strfile) ERRORQUIT("Unable to open structure file");
    parse_structure_file(&str.cell, &str.atom_pos, &str.atom_type, label, strfile, NULL);
    wrap_inside_cell(&str.atom_pos,str.atom_pos,str.cell);
  }

  // str.atom_type are indices into label;
  write_structure(str,label,axes,cout,0);

  // str.atom_type converted to indices into labellookup;
  fix_atom_type(&str, lat,labellookup,0);
  write_structure(str,lat,labellookup,label,axes,cout,0);

  // if space group needed;
  SpaceGroup spacegroup;
  spacegroup.cell=lat.cell;
  find_spacegroup(&spacegroup.point_op,&spacegroup.trans,lat.cell,lat.atom_pos,lat.atom_type);
}
