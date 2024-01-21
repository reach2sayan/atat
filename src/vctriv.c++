#include <fstream>

#include "getvalue.h"
#include "normal.h"
#include "parse.h"
#include "version.h"

char *helpstring = "Insert more help here";

int main(int argc, char *argv[]) {
  // parsing command line. See getvalue.hh for details;
  char *latfilename = "lat.in";
  int n = 0;
  int s = 100;
  Real v = 1.;
  int sigdig = 5;
  int dohelp = 0;
  int dummy = 0;
  AskStruct options[] = {
      {"", "Skeleton for atat-like codes" MAPS_VERSION ", by Axel van de Walle",
       TITLEVAL, NULL},
      {"-l", "Input file defining the lattice (defaults: lat.in)", STRINGVAL,
       &latfilename},
      {"-n", "An integer parameter", INTVAL, &n},
      {"-s", "An integer parameter", INTVAL, &s},
      {"-v", "An real parameter", REALVAL, &v},
      {"-sig", "Number of significant digits printed (default: 5)", INTVAL,
       &sigdig},
      {"-h", "Display more help", BOOLVAL, &dohelp},
      {"-d", "Use all default values", BOOLVAL, &dummy}};
  if (!get_values(argc, argv, countof(options), options)) {
    display_help(countof(options), options);
    return 1;
  }
  if (dohelp) {
    cout << helpstring;
    return 1;
  }

  cout.setf(ios::fixed);
  cout.precision(sigdig);

  // parsing lattice and structure files. See parse.hh for detail;
  Structure lat;
  Array<Arrayint> labellookup;
  Array<AutoString> label;
  rMatrix3d axes;
  {
    ifstream latfile(latfilename);
    if (!latfile) ERRORQUIT("Unable to open lattice file");
    parse_lattice_file(&lat.cell, &lat.atom_pos, &lat.atom_type, &labellookup,
		       &label, latfile, &axes);
    //    wrap_inside_cell(&lat.atom_pos,lat.atom_pos,lat.cell);
  }
  rndseed(0);
  rVector3d x = lat.atom_pos(n);
  //  cout << lat.atom_pos << endl;
  for (int i = 0; i < s; i++) {
    rVector3d nx;
    while (1) {
      nx = x + v * rVector3d(normal01(), normal01(), normal01());
      AtomPairIterator nn(lat.cell, nx, lat.atom_pos);
      if (norm(lat.atom_pos(n) - nn(1)) < zero_tolerance) break;
    }
    cout << nx << endl;
    x = nx;
  }
}
