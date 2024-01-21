#include <fstream>

#include "getvalue.h"
#include "parse.h"
#include "version.h"

char *helpstring = "Insert more help here";

int main(int argc, char *argv[]) {
  // parsing command line. See getvalue.hh for details;
  char *strfilename = "str.out";
  char *procfilename = "";
  Real r = 10;
  int nbneigh = 0;
  int sigdig = 5;
  int nbu = 40;
  int nbv = 20;
  Real slab = 5;
  Real bufwid = 0.;
  int dohelp = 0;
  int dummy = 0;
  AskStruct options[] = {
      {"", "Inverse WULFF construction " MAPS_VERSION ", by Axel van de Walle",
       TITLEVAL, NULL},
      {"-s", "Input file defining the structure (default: str.out)", STRINGVAL,
       &strfilename},
      {"-r", "radius", REALVAL, &r},
      {"-nn", "max number of neighbor for solute removal", INTVAL, &nbneigh},
      {"-nu", "number of horizontal subdivisions", INTVAL, &nbu},
      {"-nv", "number of vertical   subdivisions", INTVAL, &nbv},
      {"-pf",
       "Output file containing the processed structure (default: no output)",
       STRINGVAL, &procfilename},
      {"-st", "Slab Thickness", REALVAL, &slab},
      {"-bw", "Buffer width", REALVAL, &bufwid},
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
  Structure str;
  Array<Arrayint> labellookup;
  Array<AutoString> label;
  rMatrix3d axes;
  {
    ifstream strfile(strfilename);
    if (!strfile) ERRORQUIT("Unable to open structure file");
    parse_lattice_file(&str.cell, &str.atom_pos, &str.atom_type, &labellookup,
		       &label, strfile, &axes);
  }

  int nbcut = 0;
  for (int at = 0; at < str.atom_pos.get_size(); at++) {
    int nn = 0;
    for (int at2 = 0; at2 < str.atom_pos.get_size(); at2++) {
      if (norm(str.atom_pos(at) - str.atom_pos(at2)) < r) {
	nn++;
      }
    }
    nn--;
    //    cerr << nn << endl;
    if (nn <= nbneigh) {
      str.atom_type(at) = -1;
      nbcut++;
    }
  }
  cerr << "Solutes deleted= " << nbcut << endl;

  if (bufwid > 0) {
    Array<rVector3d> shift(3);
    rMatrix3d reclat = ~(!str.cell);
    for (int d = 0; d < 3; d++) {
      cerr << d << endl;
      rVector3d u = reclat.get_column(d);
      u.normalize();
      for (int e = 0; e < 2; e++) {
	rVector3d refpt = str.cell.get_column(d) * (Real)e;
	cerr << refpt << endl;
	int at;
	for (at = 0; at < str.atom_pos.get_size(); at++) {
	  if (str.atom_type(at) >= 0) {
	    Real d = fabs((str.atom_pos(at) - refpt) * u);
	    if (d <= bufwid) break;
	  }
	}
	if (at != str.atom_pos.get_size()) {
	  for (Real z = 0; z < fabs(str.cell.get_column(d) * u); z += bufwid) {
	    int at2;
	    for (at2 = 0; at2 < str.atom_pos.get_size(); at2++) {
	      if (str.atom_type(at2) >= 0) {
		Real d = fabs(str.atom_pos(at) * u - z);
		if (d <= bufwid) break;
	      }
	    }
	    if (at2 == str.atom_pos.get_size()) {
	      shift(d) = str.cell.get_column(d);
	      shift(d).normalize();
	      shift(d) = shift(d) * (-z);
	      cerr << "Shift along dimension " << d << " = " << shift(d)
		   << endl;
	      break;
	    }
	  }
	  break;
	} else {
	  cerr << "No shift along dimension " << d << endl;
	  shift(d) = rVector3d(0., 0., 0.);
	}
      }
    }
    rMatrix3d icell = !str.cell;
    for (int at = 0; at < str.atom_pos.get_size(); at++) {
      rVector3d v = str.atom_pos(at) + shift(0) + shift(1) + shift(2);
      str.atom_pos(at) = str.cell * (mod1(icell * v));
    }
  }
  if (strlen(procfilename) > 0) {
    ofstream file(procfilename);
    rMatrix3d iaxes = !axes;
    write_axes(axes, file, 0);
    rMatrix3d frac_cell = iaxes * str.cell;
    for (int i = 0; i < 3; i++) {
      file << frac_cell.get_column(i) << endl;
    }
    for (int i = 0; i < str.atom_pos.get_size(); i++) {
      file << (iaxes * str.atom_pos(i)) << " ";
      if (str.atom_type(i) == -1) {
	file << "Vac" << endl;
      } else {
	file << label(labellookup(str.atom_type(i))(0)) << endl;
      }
    }
  }

  rVector3d org(0., 0., 0.);
  int nbat = 0;
  for (int at = 0; at < str.atom_pos.get_size(); at++) {
    if (str.atom_type(at) >= 0) {
      nbat++;
      org += str.atom_pos(at);
    }
  }
  org = org / (Real)nbat;

  Real du = 2. * M_PI / (Real)nbu;
  Real dv = M_PI / (Real)nbv;
  for (Real u = 0; u < 2. * M_PI + du / 2.; u += du) {
    for (Real v = 0; v < M_PI + dv / 2; v += dv) {
      rVector3d dir(cos(u) * sin(v), sin(u) * sin(v), cos(v));
      Real maxp = 0;
      for (int at = 0; at < str.atom_pos.get_size(); at++) {
	if (str.atom_type(at) >= 0) {
	  maxp = max(maxp, (str.atom_pos(at) - org) * dir);
	}
      }
      int nbslab = 0;
      for (int at = 0; at < str.atom_pos.get_size(); at++) {
	if (str.atom_type(at) >= 0) {
	  Real p = (str.atom_pos(at) - org) * dir;
	  if (maxp - p <= slab) {
	    nbslab++;
	  }
	}
      }
      cout << u << "\t" << v << "\t" << maxp << "\t" << nbslab << endl;
    }
    cout << endl;
  }
}
