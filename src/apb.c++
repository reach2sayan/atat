/*
 * =====================================================================================
 *
 *       Filename:  apb.cc
 *
 *    Description:  This code creates an APB for a given structure and slip vector.
 *    							The APB energy is calculated by calling memc2 with clusters and
 *    							ECIs from a prior cluster expansion.
 *
 *         Author:  Ruoshi Sun and Axel van de Walle
 *   Organization:  School of Engineering, Brown University
 *
 * =====================================================================================
 */


#include <fstream>
#include <string>
#include "apb.h"
#include "getvalue.h"
#include "parse.h"
#include "version.h"

#define MAXMULTIPLET 6

extern const char *helpstring;

const double eV_mJ=16021.77;

	int
main ( int argc, char *argv[] )
{
  char *latfilename="lat.in";
  char *strfilename="str.out";
  char *ecifilename="eci.out";
  char *apbfilename="str_apb.out";
  char *gammafilename="gamma_apb.out";
  char *clusterfilename="clusters.out";
  char *gsfilename="gs_str.out";
  rVector3d slipvec(0, 0, 0);       // APB slip vector
  int dohelp=0;
  int dofileonly=0;
  int domc=0;
  int neq=-1;
  int nav=-1;
  Real T=0;
  int dummy=0;
  AskStruct options[]={
    {"", "\033[1;4mA\033[0mnti\033[1;4mP\033[0mhase \033[1;4mB\033[0moundary " MAPS_VERSION ", by Ruoshi Sun and Axel van de Walle\nAutomating impurity-enhanced antiphase boundary energy calculations from ab initio Monte Carlo,\n Calphad 53, 20 (2016)\n\n\033[1mFile options:\033[0m", TITLEVAL, NULL}, 
    {"-l",  "Input file: lattice   (Default: lat.in)", STRINGVAL, &latfilename}, 
    {"-s",  "          : structure (Default: str.out)", STRINGVAL, &strfilename}, 
    {"-o",  "Output file: APB structure (Default: str_apb.out)", STRINGVAL, &apbfilename}, 
    {"-og", "           : energies      (Default: gamma_apb.out)\n\n\033[1mAPB options:\033[0m", STRINGVAL, &gammafilename}, 
    {"-f",  "Generate APB structure file and exit; do not compute APB energy", BOOLVAL, &dofileonly}, 
    {"-sx", "APB slip vector: x-component", REALVAL, &slipvec(0)}, 
    {"-sy", "               : y-component", REALVAL, &slipvec(1)}, 
    {"-sz", "               : z-component (Default: 0)\n\n\033[1mMonte Carlo options:\033[0m", REALVAL, &slipvec(2)}, 
    {"-mc", "Run Monte Carlo", BOOLVAL, &domc},
    {"-eq", "Number of: equilibration passes (Default: -1)", INTVAL, &neq},
    {"-n",  "         : averaging passes     (Default: -1)", INTVAL, &nav},
    {"-T",  "Temperature (Default: 0)\n\n\033[1mOther options:\033[0m", REALVAL, &T},
    {"-d", "Use all default values", BOOLVAL, &dummy},
    {"-h",  "Display more help", BOOLVAL, &dohelp}
  };
  if ( !get_values(argc, argv, countof(options), options) ) {
    display_help(countof(options), options);
    return 1;
  }
  if (dohelp) {
    cout << helpstring;
    exit(1);
  }

  /* ----------  begin parsing ---------- */
  Structure lat;
  Array<Arrayint> labellookup;
  Array<AutoString> label;
  rMatrix3d axes;
  ifstream latfile(latfilename);
  if (!latfile) ERRORQUIT("Unable to open lattice file");
  parse_lattice_file(&lat.cell, &lat.atom_pos, &lat.atom_type, &labellookup, &label, latfile, &axes);
  wrap_inside_cell(&lat.atom_pos,lat.atom_pos,lat.cell);

  Structure str;
  ifstream strfile(strfilename);
  if (!strfile) ERRORQUIT("Unable to open structure file");
  parse_structure_file(&str.cell, &str.atom_pos, &str.atom_type, label, strfile, NULL);
  wrap_inside_cell(&str.atom_pos,str.atom_pos,str.cell);
  /* ----------  end of parsing ---------- */

  // write APB structure file
  ofstream apbfile(apbfilename);
  write_apb_structure ( str, label, axes, apbfile, slipvec );

  if (dofileonly) exit(1);

  // check if clusters.out, eci.out, gs_str.out exist
  ifstream clusterfile(clusterfilename);
  if (!clusterfile) ERRORQUIT("Unable to open clusters.out");
  ifstream ecifile(ecifilename);
  if (!ecifile) ERRORQUIT("Unable to open eci.out");
  ifstream gsfile(gsfilename);
  if (!gsfile) ERRORQUIT("Unable to open gs_str.out");

  // number of atoms
  int n_atom = str.atom_pos.get_size();
  // field number of "E_mc"
  int nf=0;

  // get cross-sectional area
  Real area = get_area(str);

  int n_type;
  n_type=label.get_size();

  // get APB enegy using clusters and ECIs
  apb_energy(strfilename, apbfilename, gammafilename, area, n_atom, n_type, nf, false, 0);

  if (domc) {
    // fix concentration in conccons.in
    char *concfilename="conccons.in";
    ofstream concfile(concfilename);
    int type;
    Array<Real> n_atom_type(n_type);
    for ( int i=0; i!=n_type; ++i ) {
      n_atom_type[i]=0.;
    }
    for ( int i=0; i!=n_atom; ++i ) {
      type=str.atom_type[i];
      ++n_atom_type[type];
    }
    for ( int i=0; i!=n_type; ++i ) {
      concfile << "1.0*" << label[i] << " = " << n_atom_type[i]/n_atom << endl;
    }
    
    char *controlfilename="control.in";
    
    // write control.in
    write_control(T, nav, n_type, controlfilename);

    // run MC
    char memc2_command[200];
    snprintf(memc2_command, 200, "memc2 -cf=%s -is=%s -n=%d -eq=%d -opss=str0000.out -keV -q",
             controlfilename,strfilename, nav, neq);
    system(memc2_command);

    // calculate APB energy
    for ( int i=0; i!=nav; ++i ) {
      char strf[30], apbf[30];
      sprintf(strf, "str%04d.out", i);
      sprintf(apbf, "str%04d_apb.out", i);
      cout << strf << apbf << endl;
      ifstream strfile(strf);
      ofstream apbfile(apbf);
      parse_structure_file(&str.cell, &str.atom_pos, &str.atom_type, label, strfile, NULL);
      wrap_inside_cell(&str.atom_pos,str.atom_pos,str.cell);
      write_apb_structure(str, label, axes, apbfile, slipvec);
      apb_energy(strf, apbf, gammafilename, area, n_atom, n_type, nf, domc, i);
    }
  }

  return 0;
}				/* ----------  end of function main  ---------- */

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  write_apb_structure
 *  Description:  Write APB structure file. Duplicate atoms along T3 and shift them by
 *  							APB slip vector.
 * =====================================================================================
 */
  void
write_apb_structure ( const Structure &str, const Array<AutoString> &atom_label,
                      const rMatrix3d &axes, ostream &file, const rVector3d &slipvec )
{
	file.setf(ios::fixed);
	file.precision(6);

	// write axes
	write_axes(axes, file, 0);

	// write translation vectors
	rMatrix3d iaxes = !axes;
	rMatrix3d frac_cell = iaxes*str.cell;
	rVector3d T3 = frac_cell.get_column(2);    // copy T3 to write atoms later
	for ( int i=0; i!=3; ++i ) {
		frac_cell(i, 2) *= 2;                    // double T3
	}
	for ( int i=0; i!=3; ++i ) {
		file << frac_cell.get_column(i) << endl; // T1, T2 remain the same
	}

	// write atoms
	for ( int i=0; i!=str.atom_pos.get_size(); ++i ) {
		file << iaxes*str.atom_pos(i) << " " << atom_label(str.atom_type(i)) << endl;
	  // copy atom to r+T3, apply slip vector
		file << iaxes*str.atom_pos(i) + T3 + slipvec << " "
			   << atom_label(str.atom_type(i)) << endl;
	}
}		/* -----  end of function write_apb_structure  ----- */

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  get_area
 *  Description:  The cross-sectional area of the supercell is given by ||T1 x T2||.
 * =====================================================================================
 */
	Real
get_area ( const Structure &str )
{
	rVector3d T1 = str.cell.get_column(0);
	rVector3d T2 = str.cell.get_column(1);
	return norm(T1^T2);
}		/* -----  end of function get_area  ----- */


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  write_control
 *  Description:  Write control.in
 * =====================================================================================
 */
void
write_control ( const Real &T, const int &n, const int &n_type, const char *controlfilename)
{
  Real T1=T+1;
  ofstream controlfile(controlfilename);
  controlfile << T;
  for (int i = 0; i != n_type; ++i) {
    controlfile << " 0";
  }
  controlfile << endl << T1;
  for (int i = 0; i != n_type; ++i) {
    controlfile << " 0";
  }
  controlfile << " " << n << endl;
}
/* -----  end of function write_control  ----- */

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  get_energy
 *  Description:  Obtain the total energy (eV/atom) of a structure using ECIs via memc2.
 * =====================================================================================
 */
Real
get_energy ( const char *strfilename, int &nf, const int &n_type, const char *mcfilename)
{
  
  char *controlfilename="energy_control.in";
  // write control.in
  write_control(0,1, n_type, controlfilename);

  char memc2_command[200];
  snprintf(memc2_command, 200, "memc2 -cf=%s -is=%s -o=%s -n=0 -eq=0 -g2c -q", controlfilename, strfilename, mcfilename);
  system(memc2_command);

  ifstream mcfile(mcfilename);
  if (!mcfile) ERRORQUIT("Unable to open monte carlo output file");

  // find field number of "E_mc"
  if (nf==0) {
    char *mcheaderfilename="mcheader.out";
    ifstream mcheaderfile(mcheaderfilename);
    if (!mcheaderfile) ERRORQUIT("Unable to open mcheader.out");

    bool found=false;
    string line;
    while(getline(mcheaderfile, line)) {
      nf++;
      if (line.find("E_mc") != string::npos) {
        found=true;
        break;
      }
    }
    if (!found) ERRORQUIT("Cannot find E_mc in mcheader.out");
  }

  // read energy
  Real energy;
  string line;
  for (int i=0; i!=nf; ++i) {
    mcfile >> energy;
  }
  return energy;
}		/* -----  end of function get_energy  ----- */

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  apb_energy
 *  Description:  Compute APB energy, defined as
 *                    2N*e1 - 2N*e0   (e1-e0)*N
 *                    ------------- = ---------
 *                         2*A            A
 *                where e0 and e1 are total energies in eV/atom.
 * =====================================================================================
 */
	void
apb_energy ( const char *strfilename, const char *apbfilename, const char *gammafilename,
	     const Real &area, const int &n_atom, const int &n_type, int &nf, bool domc, const int &i )
{
  ofstream gammafile;
  if (domc && i!=0) {
    gammafile.open(gammafilename, ofstream::out | ofstream::app);
  } else {
    gammafile.open(gammafilename);
  }
  Real e0 = get_energy(strfilename, nf, n_type, "energy_mc.out");
  Real e1 = get_energy(apbfilename, nf, n_type, "energy_mc_apb.out");
  Real gamma = (e1-e0)/area * n_atom * eV_mJ; // convert to mJ/m^2
  if (domc) {
    gammafile << setw(4) << setfill('0') << i << " ";
  }
  gammafile << e0 << " " << e1 << " " << gamma << endl;
}

/* -----  end of function apb_energy  ----- */
