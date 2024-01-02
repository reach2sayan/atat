#ifndef __APB_H__
#define __APB_H__

#include "stringo.h"
#include "xtalutil.h"
#include "vectmac.h"

void write_apb_structure ( const Structure &str, const Array<AutoString> &atom_label, const rMatrix3d &axes,
		                       ostream &file, const rVector3d &slipvec );

Real get_area ( const Structure &str );

void write_control ( const Real &T, const int &n, const int &n_type, const char *controlfilename);

Real get_energy ( const char *strfilename, int &nf, const int &n_type, const char *mcfilename);

void apb_energy ( const char *strfilename, const char *apbfilename, const char *gammafilename, 
                  const Real &area, const int &n_atom, const int &n_type, int &nf, bool append, const int &i );

#endif
