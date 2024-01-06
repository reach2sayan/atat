#include "getvalue.h"
#include "parse.h"
#include "version.h"
#include "xtalutil.h"
#include <fstream.h>

char *helpstring = "Insert more help here";

int equivalent_mod_cell(const rVector3d &a, const rVector3d &b,
                        const rMatrix3d &inv_cell, Real tol) {
  rVector3d frac_a = inv_cell * a;
  rVector3d frac_b = inv_cell * b;
  rVector3d delta = frac_a - frac_b;
  for (int i = 0; i < 3; i++) {
    delta(i) = cylinder_norm(delta(i));
  }
  return (norm(delta) < tol ? 1 : 0);
}

int which_atom(const Array<rVector3d> &atom_pos, const rVector3d &pos,
               const rMatrix3d &inv_cell, Real tol) {
  for (int i = 0; i < atom_pos.get_size(); i++) {
    if (equivalent_mod_cell(atom_pos(i), pos, inv_cell, tol))
      return i;
  }
  return -1;
}

Real angle(const rVector3d &a, const rVector3d &b) {
  return acos(min(max(cos_angle(a, b), -1.), 1.));
}

int struct_compare(Array<int> *pmatches, const Structure &str1,
                   const Structure &str2_, Real ftol, Real vtol, Real ltol) {
  Structure str2;

  Real ratio = fabs(det(str1.cell) / det(str2_.cell));
  int lonratio = max((int)round(ratio * (1 - vtol)), 1);
  int hinratio = max((int)round(ratio * (1 + vtol)), 1);
  for (int nratio = lonratio; nratio <= hinratio; nratio++) {
    Real scale = pow(ratio / nratio, 1. / 3.);
    str2.cell = str2_.cell * scale;
    str2.atom_pos.resize(str2_.atom_pos.get_size());
    str2.atom_type = str2_.atom_type;
    for (int at = 0; at < str2_.atom_type.get_size(); at++) {
      str2.atom_pos(at) = str2_.atom_pos(at) * scale;
    }

    rMatrix3d istr2cell = (!str2.cell);
    pmatches->resize(str2.atom_pos.get_size());
    zero_array(pmatches);

    rVector3d len1;
    for (int i = 0; i < 3; i++) {
      len1(i) = norm(str1.cell.get_column(i));
    }
    LatticePointIterator lp0(str2.cell);
    for (; norm(lp0) <= len1(0) * (1. + 2. * ltol); lp0++) {
      //          cout << "0:" << (rVector3d)lp0 << endl;
      if (fabs(norm(lp0) - len1(0)) < ltol * len1(0)) {
        LatticePointIterator lp1(str2.cell);
        for (; norm(lp1) <= len1(1) * (1. + 2. * ltol); lp1++) {
          //	      cout << "1:" << (rVector3d)lp1 << endl;
          if (fabs(norm(lp1) - len1(1)) < ltol * len1(1)) {
            if (fabs(angle(str1.cell.get_column(0), str1.cell.get_column(1)) -
                     angle((rVector3d)lp0, (rVector3d)lp1)) < ltol) {
              LatticePointIterator lp2(str2.cell);
              for (; norm(lp2) <= len1(2) * (1. + 2 * ltol); lp2++) {
                //	            cout << "2:" << (rVector3d)lp2 << endl;
                if (fabs(norm(lp2) - len1(2)) < ltol * len1(2)) {
                  if (fabs(angle(str1.cell.get_column(1),
                                 str1.cell.get_column(2)) -
                           angle((rVector3d)lp1, (rVector3d)lp2)) < ltol &&
                      fabs(angle(str1.cell.get_column(0),
                                 str1.cell.get_column(2)) -
                           angle((rVector3d)lp0, (rVector3d)lp2)) < ltol &&
                      det(str1.cell) * ((rVector3d)lp0 *
                                        ((rVector3d)lp1 ^ (rVector3d)lp2)) >
                          0.) {
                    Structure rstr1;
                    rstr1.cell.set_column(0, lp0);
                    rstr1.cell.set_column(1, lp1);
                    rstr1.cell.set_column(2, lp2);
                    rstr1.atom_pos.resize(str1.atom_pos.get_size());
                    rstr1.atom_type.resize(str1.atom_pos.get_size());
                    rMatrix3d rot = rstr1.cell * (!str1.cell);
                    for (int at = 0; at < str1.atom_type.get_size(); at++) {
                      rstr1.atom_pos(at) = rot * str1.atom_pos(at);
                      rstr1.atom_type(at) = str1.atom_type(at);
                    }
                    for (int at2 = 0; at2 < str2.atom_type.get_size(); at2++) {
                      if (str2.atom_type(at2) == rstr1.atom_type(0)) {
                        zero_array(pmatches);
                        rVector3d t = str2.atom_pos(at2) - rstr1.atom_pos(0);
                        int at1;
                        for (at1 = 0; at1 < rstr1.atom_type.get_size(); at1++) {
                          int at =
                              which_atom(str2.atom_pos, rstr1.atom_pos(at1) + t,
                                         istr2cell, ftol);
                          if (at == -1)
                            break;
                          if (str2.atom_type(at) != rstr1.atom_type(at1))
                            break;
                          (*pmatches)(at)++;
                        }
                        if (at1 == rstr1.atom_type.get_size()) {
                          int at;
                          for (at = 0; at < pmatches->get_size(); at++) {
                            if ((*pmatches)(at) == nratio)
                              break;
                          }
                          if (at < pmatches->get_size()) {
                            return 1;
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  return 0;
}

int main(int argc, char *argv[]) {
  // parsing command line. See getvalue.hh for details;
  char *str1filename = "str1.out";
  char *str2filename = "str2.out";
  Real ftol = 1e-2;
  Real vtol = 0.25;
  Real ltol = 0.25;

  int sigdig = 5;
  int dohelp = 0;
  int dummy = 0;
  AskStruct options[] = {
      {"", "CoMPare STRuctures " MAPS_VERSION ", by Axel van de Walle",
       TITLEVAL, NULL},
      {"-s1", "Input file defining the lattice (defaults: str1.out)", STRINGVAL,
       &str1filename},
      {"-s2", "Input file defining the structure (default: str2.out)",
       STRINGVAL, &str2filename},
      {"-ft", "Fractional Tolerance", REALVAL, &ftol},
      {"-vt", "Volume Tolerance", REALVAL, &vtol},
      {"-lt", "Length Tolerance", REALVAL, &ltol},
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

  Array<AutoString> label;
  Structure str1;
  {
    Array<Arrayint> labellookup;
    rMatrix3d axes;
    ifstream latfile(str1filename);
    if (!latfile)
      ERRORQUIT("Unable to open lattice file");
    parse_lattice_file(&str1.cell, &str1.atom_pos, &str1.atom_type,
                       &labellookup, &label, latfile, &axes);
    wrap_inside_cell(&str1.atom_pos, str1.atom_pos, str1.cell);
    for (int at = 0; at < str1.atom_type.get_size(); at++) {
      str1.atom_type(at) = labellookup(str1.atom_type(at))(0);
    }
  }

  Structure str2;
  {
    ifstream strfile(str2filename);
    if (!strfile)
      ERRORQUIT("Unable to open structure file");
    parse_structure_file(&str2.cell, &str2.atom_pos, &str2.atom_type, label,
                         strfile, NULL);
    wrap_inside_cell(&str2.atom_pos, str2.atom_pos, str2.cell);
  }
  Array<int> matches;
  if (struct_compare(&matches, str1, str2, ftol, vtol, ltol)) {
    cout << matches << endl
         << str1.atom_type.get_size() << endl
         << str2.atom_type.get_size() << endl;
  } else {
    cout << "No_match" << endl;
  }
}
