#include <fstream>
#include "parse.h"
#include "getvalue.h"
#include "version.h"
#include "xtalutil.h"

void read_atom_clus(Array<rVector3d> *patom_pos, Array<AutoString> *plabel, istream &file) {
  char *delim=" \t";
  LinkedList<rVector3d> latom_pos;
  LinkedList<AutoString> llabel;
  while (1) {
    rVector3d pos(MAXFLOAT,MAXFLOAT,MAXFLOAT);
    file >> pos;
    if (file.eof() || pos(0)==MAXFLOAT || !file) break;
    latom_pos << new rVector3d(pos);
    char buf[MAX_LINE_LEN];
    file.get(buf,MAX_LINE_LEN-1);
    char *atom_begin=buf;
    char *string_end=buf+strlen(buf);
    while (strchr(delim,*atom_begin) && atom_begin<string_end) {atom_begin++;}
    char *atom_end=atom_begin;
    while (!strchr(delim,*atom_end) && atom_begin<string_end) {atom_end++;}
    *atom_end=0;
    llabel << new AutoString(atom_begin);
  }
  LinkedList_to_Array(patom_pos,latom_pos);
  LinkedList_to_Array(plabel,llabel);
}

// write extra help as plain text in *.hlp
char *helpstring="";

int main(int argc, char *argv[]) {
  // parsing command line. See getvalue.hh for details;
  char *latfilename="lat.in";
  char *clusaxesfilename="";
  const char *clusfilename[2]={"",""};
  int dohelp=0;
  int dummy=0;
  AskStruct options[]={
    {"","Skeleton for atat-like codes " MAPS_VERSION ", by Axel van de Walle",TITLEVAL,NULL},
    {"-h","Display more help",BOOLVAL,&dohelp},
    {"-l","Input file defining the lattice (Default: lat.in)",STRINGVAL,&latfilename},
    {"-a","Axes input file",STRINGVAL,&clusaxesfilename},
    {"-c1","Cluster 1 file",STRINGVAL,&(clusfilename[0])},
    {"-c2","Cluster 2 file",STRINGVAL,&(clusfilename[1])},
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
  SpaceGroup spacegroup;
  spacegroup.cell=lat.cell;
  find_spacegroup(&spacegroup.point_op,&spacegroup.trans,lat.cell,lat.atom_pos,lat.atom_type);

  rMatrix3d clusaxes=axes;
  if (strlen(clusaxesfilename)>0) {
    ifstream clusaxesfile(clusaxesfilename);
    if (!clusaxesfile) ERRORQUIT("Unable to read axes.");
    read_cell(&clusaxes,clusaxesfile);
  }

  Array<Array<rVector3d> > cluspos(2);
  Array<Array<AutoString> > cluslabel(2);
  for (int i=0; i<2; i++) {
    ifstream clusfile(clusfilename[i]);
    if (!clusfile) ERRORQUIT("Unable to read cluster.");
    read_atom_clus(&(cluspos(i)),&(cluslabel(i)),clusfile);
    for (int at=0; at<cluspos(i).get_size(); at++) {
      cluspos(i)(at)=clusaxes*cluspos(i)(at);
    }
  }

  /*
  for (int i=0; i<2; i++) {
    cerr << i << endl;
    cerr << cluspos(i) << endl;
    cerr << cluslabel(i) << endl;
  }
  */

  rMatrix3d invcell=!lat.cell;
  int op;
  for (op=0; op<spacegroup.point_op.get_size(); op++) {
    Array<rVector3d> transpos0;
    apply_symmetry(&transpos0, spacegroup.point_op(op),spacegroup.trans(op),cluspos(0));
    rVector3d shift;
    if (equivalent_mod_cell(cluspos(1),transpos0,invcell,&shift)) {
      // cerr << op << endl;
      // cerr << transpos0 << endl;
      int at;
      for (at=0; at<transpos0.get_size(); at++) {
	// cerr << "--" << endl;
	// cerr << cluslabel(1)(which_atom(cluspos(1),transpos0(at)-shift)) << endl;
	// cerr << cluslabel(0)(at) << endl;
	if (strcmp(cluslabel(1)(which_atom(cluspos(1),transpos0(at)-shift)),cluslabel(0)(at)) != 0) break;
      }
      if (at==transpos0.get_size()) break;
    }
  }
  return (op==spacegroup.point_op.get_size());
}
