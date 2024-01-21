#include <fstream>
#include "parse.h"
#include "getvalue.h"
#include "version.h"

extern const char *helpstring;

int main(int argc, char *argv[]) {
  // parsing command line. See getvalue.hh for details;
  const char *latfilename="-";
  Real r=0;
  int doerror=0;
  int dohelp=0;
  AskStruct options[]={
    {"","Test if MINimum DISTance between two atoms is less than a given threshold, version " MAPS_VERSION ", by Axel van de Walle",TITLEVAL,NULL},
    {"-s","Input file defining the structure (Default: stdin)",STRINGVAL,&latfilename},
    {"-r","Threshold radius",REALVAL,&r},
    {"-e","Create error file atoms are too close",BOOLVAL,&doerror},
    {"-h","Display more help",BOOLVAL,&dohelp}
  };
  if (!get_values(argc,argv,countof(options),options)) {
    display_help(countof(options),options);
    return 1;
  }
  if (dohelp) {
    cout << helpstring;
    return 1;
  }

  // parsing lattice and structure files. See parse.hh for detail;
  Structure lat;
  Array<Arrayint> labellookup;
  Array<AutoString> label;
  rMatrix3d axes;
  {
    if (strcmp(latfilename,"-")==0) {
      parse_lattice_file(&lat.cell, &lat.atom_pos, &lat.atom_type, &labellookup, &label, cin, &axes);
    }
    else {
      ifstream latfile(latfilename);
      if (!latfile) ERRORQUIT("Unable to open lattice file");
      parse_lattice_file(&lat.cell, &lat.atom_pos, &lat.atom_type, &labellookup, &label, latfile, &axes);
    }
    wrap_inside_cell(&lat.atom_pos,lat.atom_pos,lat.cell);
  }

  AtomPairIterator pair(lat.cell,lat.atom_pos);
  if (pair.length()<r) {
    if (doerror) {
      ofstream file("error");
    }
    return 1;
  }
  else {
    return 0;
  }
}
