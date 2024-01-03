#include <fstream>
#include "clus_str.h"
#include "parse.h"
#include "getvalue.h"
#include "version.h"
#include "calccorr.h"

extern const char *helpstring;

void generate_permutation(Array<int> *pperm, int n) {
  pperm->resize(n);
  for (int i=0; i<n; i++) {(*pperm)(i)=-1;}
  for (int d=0; d<n; d++) {
    int i=0;
    int r=random(n-d);
    while (1) {
      while ((*pperm)(i)!=-1) {i++;}
      if (r==0) break;
      i++;
      r--;
    }
    (*pperm)(i)=d;
  }
}

void calc_all_corr(Array<Real> *pcorr, const Structure &str, const Array<Array<MultiCluster> > &eqclus, const rMatrix3d unitcell, const Array<Array<Array<Real> > > &corrfunc) {
  pcorr->resize(eqclus.get_size());
  for (int t=0; t<eqclus.get_size(); t++) {
    (*pcorr)(t)=calc_correlation(str,eqclus(t),unitcell,corrfunc);
  }
}

/*
Real calc_objective_func(const Array<Real> &corr, const Array<Real> &tcorr, Real mysqstol, Real weightdist, const Array<Real> &diam, const Array<int> &nbpt) {
  Array<Real> dcorr(tcorr.get_size());
  for (int t=0; t<tcorr.get_size(); t++) {
    dcorr(t)=fabs(corr(t)-tcorr(t));
  }
  Real objdist=MAXFLOAT;
  Real objdev=0.;
  for (int t=0; t<tcorr.get_size(); t++) {
    if (dcorr(t)>mysqstol) {
      objdist=min(diam(t),objdist);
      objdev+=dcorr(t);
    }
  }
  if (objdist==MAXFLOAT) {return -MAXFLOAT;}
  return -objdist*weightdist+objdev;
}
*/

Real calc_objective_func(const Array<Real> &corr, const Array<Real> &tcorr, Real mysqstol, Real weightdist, Real weightnbpt, Real weightdecay, const Array<Real> &diam, const Array<int> &nbpt) {
  Array<Real> dcorr(tcorr.get_size());
  Array<Real> maxdist(max(nbpt)-1);
  zero_array(&maxdist);
  for (int t=0; t<tcorr.get_size(); t++) {
    dcorr(t)=fabs(corr(t)-tcorr(t));
    maxdist(nbpt(t)-2)=max(maxdist(nbpt(t)-2),diam(t));
  }
  Real d0=MAXFLOAT;
  for (int t=0; t<tcorr.get_size(); t++) {
    d0=min(diam(t),d0);
  }
  for (int p=0; p<maxdist.get_size(); p++) {
    maxdist(p)+=d0;
  }
  for (int t=0; t<tcorr.get_size(); t++) {
    if (dcorr(t)>mysqstol) {
      maxdist(nbpt(t)-2)=min(diam(t),maxdist(nbpt(t)-2));
    }
  }
  Real d1=maxdist(0);
  for (int p=1; p<maxdist.get_size(); p++) {
    maxdist(p)=min(maxdist(p),maxdist(p-1));
    d1=min(maxdist(p),d1);
  }
  Real objdev=0.;
  Real den=0.;
  for (int t=0; t<tcorr.get_size(); t++) {
    if (diam(t)>=d1-zero_tolerance) {
      Real w=exp(-weightdecay*diam(t)/d0)*pow(weightnbpt,nbpt(t)-2);
      objdev+=dcorr(t)*w;
      den+=w;
    }
  }
  if (near_zero(objdev)) {return -MAXFLOAT;}
  objdev/=den;

  Real obj=objdev;
  for (int p=0; p<maxdist.get_size(); p++) {
    obj-=weightdist*pow(weightnbpt,p)*maxdist(p)/d0;
  }
  return obj;
}

void open_numbered_file(ofstream &file, const char *prefix, int num, const char *suffix) {
  ostringstream filename;
  if (num>-1) {
    filename << prefix << num << suffix << '\0';
  }
  else {
    filename << prefix << suffix << '\0';
  }
  file.open(filename.str().c_str());
  if (!file) {
    cerr << "Unable to open " << filename.str().c_str() << endl;
    ERRORQUIT("Aborting");
  }
}

class SupercellData {
public:
  SupercellData(void): str(), nbcomp(), corr(), obj() {};
  Structure str;
  Array<int> nbcomp;
  Array<int> typebeg;
  Array<int> typeend;
  int nbactive;
  Array<Real> corr;
  Real obj;
public:
  void operator=(const SupercellData &mc) {
    str=mc.str;
    nbcomp=mc.nbcomp;
    typebeg=mc.typebeg;
    typeend=mc.typeend;
    nbactive=mc.nbactive;
    corr=mc.corr;
    obj=mc.obj;
  }
};

#define MAXMULTIPLET 6

int main(int argc, char *argv[]) {
  int maxvol=0;
  Array<Real> maxd(MAXMULTIPLET+1);
  zero_array(&maxd);
  const char *rndstrfilename="rndstr.in";
  const char *clusterfilename="clusters.out";
  const char *tcorrfilename="";
  const char *paramfilename="sqsparam.in";
  int readcell=0;
  int ip=-1;
  int findbest=0;
  Real mysqstol=zero_tolerance;
  Real weightdist=1.;
  Real weightnbpt=1.;
  Real weightdecay=0.;
  Real T=1.;
  int seed=0;
  int maxtic=10000;
  int do2d=0;
  int sigdig=6;
  const char *corrfunc_label="trigo";
  int dohelp=0;
  AskStruct options[]={
    {"","Monte Carlo generator of Special Quasirandom Structures " MAPS_VERSION ", by Axel van de Walle",TITLEVAL,NULL},
    {"-n","nb of atom/unit cell",INTVAL,&maxvol},
    {"","If -n is not specified, generate clusters of the following sizes:",TITLEVAL,NULL},
    {"-2"," Maximum distance between two points within a pair",REALVAL,&maxd(2)},
    {"-3"," Maximum distance between two points within a triplet",REALVAL,&maxd(3)},
    {"-4"," Maximum distance between two points within a quadruplet",REALVAL,&maxd(4)},
    {"-5"," Maximum distance between two points within a quintuplet",REALVAL,&maxd(5)},
    {"-6"," Maximum distance between two points within a sextuplet",REALVAL,&maxd(6)},
    {"-l","Input file defining the random structure (default: rndstr.in)",STRINGVAL,&rndstrfilename},
    {"-cf","Input file defining the clusters (default: clusters.out)",STRINGVAL,&clusterfilename},
    {"-tcf","Input file defining the target multibody correlations (default: internally calculated values for fully disordered state)",STRINGVAL,&tcorrfilename},
    {"-tol","Tolerance for matching correlations (default: 1e-3)",REALVAL,&mysqstol},
    {"-wr","Weight assigned to range of perfect correlation match in objective function (default 1)",REALVAL,&weightdist},
    {"-wn","Multiplicative decrease in weight per additional point in cluster (default 1)",REALVAL,&weightnbpt},
    {"-wd","Exponent of decay in weight as function of cluster diameter (default 0)",REALVAL,&weightdecay},
    {"-T","Temperature (default 1)",REALVAL,&T},
    {"-pf","Input file defining the optimization parameters (default: sqsparam.in)",STRINGVAL,&paramfilename},
    {"-rc","Read supercells from file sqscell.out (default: generate internally and write to sqscell.out)",BOOLVAL,&readcell},
    {"-ip","Index of current process (for parallel operation)",INTVAL,&ip},
    {"-best","Collect best SQS among the outputs of prior parallel runs",BOOLVAL,&findbest},
    {"-crf","Select correlation functions (default: trigo)",STRINGVAL,&corrfunc_label},
    {"-sd","Seed for random number generation (default: use clock)",INTVAL,&seed},
    {"-rt","Read parameter file every (rt) step (default:10000)",INTVAL,&maxtic},
    {"-2d","Generate only supercells in the plane of a,b axes",BOOLVAL,&do2d},
    {"-sig","Number of significant digits to print in output files (default: 6)",INTVAL,&sigdig},
    {"-h","Display more help",BOOLVAL,&dohelp}
  };
  if (!get_values(argc,argv,countof(options),options)) {
    display_help(countof(options),options);
    if (argc==1) {cerr << endl << "Perhaps missing -n=[nb of atom/cell] option?" << endl;}
    return 1;
  }
  if (dohelp) {
    cout << helpstring;
    return 1;
  }

  if (findbest) {
    {
      ofstream script("findbestsqs");
      script << "#!/bin/sh" << endl 
	     << "stem=$(grep -n Objective bestcorr?*.out | sort -n -k 2 | head -1 | sed 's/:.*$//g')" << endl
	     << "if [ \"${stem}\" == \"\" ]; then" << endl
	     << "  echo No output yet." << endl
	     << "  exit 1" << endl
	     << "fi" << endl
	     << "cp $stem bestcorr.out" << endl
	     << "cp $(echo $stem | sed 's/corr/sqs/g') bestsqs.out" << endl
	     << "cat bestcorr.out" << endl;
    }
    system("chmod +x ./findbestsqs ; sh -c ./findbestsqs");
    exit(0);
  }
  if (maxvol==0 && !readcell) {
    cout << "Generating clusters..." << endl;
    ostringstream cmd;
    cmd << "corrdump -clus -noe -nop -ro -l=" << rndstrfilename << " ";
    for (int i=2; i<=MAXMULTIPLET; i++) {
      if (maxd(i)!=0.) {
	cmd << "-" << i << "=" << maxd(i) << " ";
      }
    }
    cmd << "; getclus" << '\0';
    //    cout << cmd.str().c_str() << endl;
    system(cmd.str().c_str());
    exit(0);
  }

  if (maxvol==0 && !readcell) {ERRORQUIT("Please specify -n option.");}

  cerr.setf(ios::fixed);
  cerr.precision(sigdig);

  // read in lattice (see parse.h);
  Structure ulat;
  Array<Array<Real> > uatomprob;
  Array<Arrayint> site_type_list;
  Array<AutoString> atom_label;
  rMatrix3d axes;
  {
    ifstream file(rndstrfilename);
    if (!file) ERRORQUIT("Unable to open random structure file.");
    parse_rndstr_file(&ulat.cell, &ulat.atom_pos, &ulat.atom_type, &uatomprob, &site_type_list, &atom_label, file, &axes);
  }
  rMatrix3d inv_cell=!ulat.cell;
  rMatrix3d inv_axes=!axes;

  SpaceGroup spacegroup;
  find_spacegroup(&spacegroup.point_op,&spacegroup.trans,ulat.cell,ulat.atom_pos,ulat.atom_type);
  spacegroup.cell=ulat.cell;

  // initialize a table of correlation functions;
  if (!check_plug_in(CorrFuncTable(),corrfunc_label)) {
    ERRORQUIT("Aborting");
  }
  CorrFuncTable *pcorrfunc=GenericPlugIn<CorrFuncTable>::create(corrfunc_label);
  pcorrfunc->init_from_site_type_list(site_type_list);

  // group sites that are equivalent by symmetry;
  Structure lat;
  Array<int> sym_to_type;
  Array<Array<Real> > atomprob;
  Array<Array<Real> > sym_type_prob;
  {
    lat.cell=ulat.cell;
    lat.atom_pos.resize(ulat.atom_pos.get_size());
    lat.atom_type.resize(ulat.atom_type.get_size());
    atomprob.resize(ulat.atom_type.get_size());
    ClusterBank pt_clus(ulat.cell,ulat.atom_pos,1,spacegroup);
    sym_to_type.resize(pt_clus.get_cluster_list().get_size());
    sym_type_prob.resize( sym_to_type.get_size());
    LinkedListIterator<Cluster> i_pt_clus(pt_clus.get_cluster_list());
    int at=0;
    for (int symt=0; i_pt_clus; i_pt_clus++,symt++) {
      int uat=which_atom(ulat.atom_pos,(*i_pt_clus)(0),inv_cell);
      sym_to_type(symt)=ulat.atom_type(uat);
      Array<Real> cur_atomprob;
      sym_type_prob(symt)=uatomprob(uat);
      Array<Cluster> pt_clus_eq;
      find_equivalent_clusters(&pt_clus_eq, *i_pt_clus, spacegroup.cell, spacegroup.point_op, spacegroup.trans);
      for (int eq=0; eq<pt_clus_eq.get_size(); eq++) {
	lat.atom_pos(at)=pt_clus_eq(eq)(0);
	lat.atom_type(at)=symt;
	Array<Real> dprob;
	atomprob(at)=sym_type_prob(symt);
	int uat=which_atom(ulat.atom_pos,pt_clus_eq(eq)(0),inv_cell);
	diff(&dprob,sym_type_prob(symt),uatomprob(uat));
	if (norm(dprob)>zero_tolerance) {
	  cerr << "Different occupations assigned to symmetrically equivalent site:" << endl;
	  cerr << inv_axes*(*i_pt_clus)(0) << endl << "and" << endl << inv_axes*pt_clus_eq(eq)(0) << endl;
	  ERRORQUIT("Aborting.");
	}
	at++;
      }
    }
/*    int uat=0;
    while (uat<ulat.atom_type.get_size()) {
      if (site_type_list(ulat.atom_type(uat)).get_size()==1) break;
      uat++;
    }
    for (; at<lat.atom_type.get_size(); at++,uat++) {
      lat.atom_pos(at)=ulat.atom_pos(uat);
      lat.atom_type(at)=sym_to_type.get_size()-1;
      atomprob(at).resize(1);
      atomprob(at)(0)=1.;
    }
    sym_to_type(sym_to_type.get_size()-1)=-1;
    for (int st=0; st<site_type_list.get_size(); st++) {
      if (site_type_list(st).get_size()==1) {
	sym_to_type(sym_to_type.get_size()-1)=st;
	sym_type_prob(sym_to_type.get_size()-1).resize(1);
	sym_type_prob(sym_to_type.get_size()-1)(0)=1.;
	break;
      }
    }
*/ 
 }
  if (ip<=0) {
    ofstream groupfile("rndstrgrp.out");
    groupfile.setf(ios::fixed);
    groupfile.precision(sigdig);
    write_axes(axes,groupfile,0);
    rMatrix3d frac_cell=inv_axes*lat.cell;
    for (int i=0; i<3; i++) {
      groupfile << frac_cell.get_column(i) << endl;
    }
    int last_type=-1;
    for (int i=0; i<lat.atom_pos.get_size(); i++) {
      if (last_type!=lat.atom_type(i)) {groupfile << endl;}
      last_type=lat.atom_type(i);
      groupfile << (inv_axes*lat.atom_pos(i)) << " ";
      //      int site_type=ulat.atom_type(which_atom(ulat.atom_pos,lat.atom_pos(i),inv_cell));
      int site_type=sym_to_type(lat.atom_type(i));
      for (int j=0; j<site_type_list(site_type).get_size(); j++) {
	if (j!=0) {
	  groupfile  << ",";
	}
	groupfile << atom_label(site_type_list(site_type)(j)) << "=" << sym_type_prob(lat.atom_type(i))(j);
      }
      groupfile << endl;
    }
  }

  //  TRACEIT(ulat.atom_pos);
  //  TRACEIT(ulat.atom_type);
  //  TRACEIT(lat.atom_pos);
  //  TRACEIT(lat.atom_type);

  Array<Array<Real> > pointcorr(lat.atom_type.get_size());
  for (int at=0; at<lat.atom_type.get_size(); at++) {
    int comp=atomprob(at).get_size();
    pointcorr(at).resize(comp-1);
    for (int func=0; func<comp-1; func++) {
      Real corr=0.;
      Real totprob=0;
      for (int p=0; p<atomprob(at).get_size(); p++) {
	corr+=atomprob(at)(p)*(*pcorrfunc)(comp-2)(func)(p);
	if (atomprob(at)(p)<0) {ERRORQUIT("Negative occupations!");}
	if (atomprob(at).get_size()>1 && near_zero(atomprob(at)(p)-1)) {
	  ERRORQUIT("Please drop irrelevant species when site is fully occupied with one specie");
	}
	totprob+=atomprob(at)(p);
      }
      if (!near_zero(totprob-1.)) ERRORQUIT("Probabilities must sum to 1.");
      pointcorr(at)(func)=corr;
    }
  }
  // cerr << pointcorr << endl;
  LinkedList<MultiCluster> clusterlist;
  {
    ifstream clusterfile(clusterfilename);
    if (!clusterfile) ERRORQUIT("Unable to open cluster file.  Please generate it with\n mcsqs -2=[max pair diameter] -3=[max triplet diameter] etc.");
    read_clusters_and_eci(&clusterlist, NULL, clusterfile, clusterfile, axes);
  }

  ifstream tcorrfile;
  if (strlen(tcorrfilename)>0) {
    tcorrfile.open(tcorrfilename);
    if ( !tcorrfile ) {ERRORQUIT("Unable to open correlation file");}
  }
  
  LinkedList<Real> tcorr_list;
  LinkedListIterator<MultiCluster> icluster(clusterlist);
  for (; icluster; ) {
    if (icluster->clus.get_size()<=1) {
      delete clusterlist.detach(icluster);
    }
    else {
      Real corr;
      if (strlen(tcorrfilename)>0) {
	tcorrfile >> corr;
      }
      else {
	corr=1.;
	for (int i=0; i<icluster->clus.get_size(); i++) {
	  int at=which_atom(lat.atom_pos,icluster->clus(i),inv_cell);
	  if (at==-1) {ERRORQUIT("Cluster file inconsistent with random structure file.");}
	  corr*=pointcorr(at)(icluster->func(i));
	}
      }
      tcorr_list << new Real(corr);
      icluster++;
    }
  }
  Array<Real> tcorr;
  LinkedList_to_Array(&tcorr,tcorr_list);
  //  cerr << tcorr;

  cout.setf(ios::fixed);
  cout.precision(sigdig);

  Array<Array<MultiCluster> > eqclus(clusterlist.get_size());
  Array<Real> diam(clusterlist.get_size());
  Array<int> nbpt(clusterlist.get_size());
  LinkedListIterator<MultiCluster> ic(clusterlist);
  for (int t=0; ic; t++, ic++) {
    find_equivalent_clusters(&eqclus(t),*ic,lat.cell,spacegroup.point_op,spacegroup.trans);
    diam(t)=get_length_quick(ic->clus);
    nbpt(t)=ic->clus.get_size();
  }

  ofstream logfile;
  open_numbered_file(logfile, "mcsqs",ip,".log");
  logfile.setf(ios::fixed);
  logfile.precision(sigdig);


  Array<rMatrix3d> supercell;
  int v=(int)(maxvol/MAX(1,lat.atom_pos.get_size()));
  if (readcell) {
    logfile << "Reading supercells..." << endl <<  flush;
    ifstream cellfile("sqscell.out");
    if (!cellfile) {
      ERRORQUIT("Unable to read sqscell.out");
    }
    int nc=0;
    cellfile >> nc;
    supercell.resize(nc);
    for (int i=0; i<supercell.get_size(); i++) {
      rMatrix3d tmpcell;
      read_cell(&tmpcell,cellfile);
      supercell(i)=axes*tmpcell;
    }
  }
  else {
    logfile << "Generating supercells..." << endl <<  flush;
    rMatrix3d iaxes=!axes;
    Array<rMatrix3d> pointgroup;
    pointgroup_from_spacegroup(&pointgroup, spacegroup.point_op);
    if (do2d) {
      find_supercells_2D(&supercell, v, v, lat.cell, pointgroup);
    }
    else {
      find_supercells(&supercell, v, v, lat.cell, pointgroup);
    }
    if (ip<=0) {
      logfile << "Writing supercells..." << endl <<  flush;
      ofstream cellfile("sqscell.out");
      cellfile.setf(ios::fixed);
      cellfile.precision(sigdig);
      cellfile << supercell.get_size() << endl << endl;
      for (int i=0; i<supercell.get_size(); i++) {
	//cout << supercell(i) << " ";
	//cout << shorten_lattice_vectors(supercell(i)) << endl;
	supercell(i)=find_symmetric_cell(supercell(i));
	write_axes(iaxes*supercell(i),cellfile,0);
	cellfile << endl;
      }
    }
    else {
      logfile << "Not writing supercells (process -ip=0 is doing that)." << endl <<  flush;
    }
  }

  rndseed(seed);

  Array<SupercellData> mc(supercell.get_size());
  SupercellData best;
  best.obj=MAXFLOAT;
  int cc=0;

  logfile << "Initializing random supercells..." << endl << flush;
  for (int c=0; c<supercell.get_size(); c++) {
    mc(c).str.cell=supercell(c);
    find_all_atom_in_supercell_ordered(&(mc(c).str.atom_pos),
		 &(mc(c).str.atom_type),lat.atom_pos,
		 lat.atom_type,
		 lat.cell, mc(c).str.cell);
    Array<int> curtype=mc(c).str.atom_type;
    int curbeg=0;
    int curend=0;
    mc(c).nbcomp.resize(curtype.get_size());
    mc(c).typebeg.resize(curtype.get_size());
    mc(c).typeend.resize(curtype.get_size());
    for (int at=0; at<curtype.get_size(); at++) {
      mc(c).nbcomp(at)=site_type_list(sym_to_type(curtype(at))).get_size();
      if (at==curend) {
	curbeg=curend;
	while (curend<curtype.get_size()) {
	  if (curtype(at)!=curtype(curend)) break;
	  curend++;
	}
      }
      mc(c).typebeg(at)=curbeg;
      mc(c).typeend(at)=curend;
    }
    mc(c).nbactive=curtype.get_size();
    for (int at=0; at<mc(c).str.atom_pos.get_size(); ) {
      int nbat=mc(c).typeend(at)-mc(c).typebeg(at);
      Array<int> perm;
      generate_permutation(&perm,nbat);
      if (mc(c).nbcomp(at)==1) {
	if (mc(c).nbactive==curtype.get_size()) {
	  mc(c).nbactive=at;
	}
      }
      int at2=0;
      for (int t=0; t<mc(c).nbcomp(at); t++) {
	Real rnum=(Real)nbat*sym_type_prob(curtype(at))(t);
	int inum=(int)round(rnum);
	//cerr << rnum << " " << inum << endl;
	if (!near_zero((Real)inum-rnum)) ERRORQUIT("Impossible to match point correlations due to incompatible supercell size.");
	for (int i=0; i<inum; i++) {
	  mc(c).str.atom_type(at+perm(at2))=t;
	  at2++;
	}
      }
      at=mc(c).typeend(at);
    }
    calc_all_corr(&(mc(c).corr),mc(c).str,eqclus,lat.cell,*pcorrfunc);
    mc(c).obj=calc_objective_func(mc(c).corr,tcorr,mysqstol,weightdist,weightnbpt,weightdecay,diam,nbpt);
    if (mc(c).obj<best.obj) {
      best=mc(c);
      cc=c;
    }
  }
  logfile << "Initialization done." << endl;

  int tic=0;
  Real obj=best.obj;
  best.obj=MAXFLOAT;
  while (1) {
    if (obj<best.obj) {
      // cerr << "Best" << endl;
      best=mc(cc);
      ofstream strfile;
      open_numbered_file(strfile, "bestsqs",ip,".out");
      strfile.setf(ios::fixed);
      strfile.precision(sigdig);
      write_structure(best.str,ulat,site_type_list,atom_label,axes,strfile);
      ofstream corrfile;
      open_numbered_file(corrfile, "bestcorr",ip,".out");
      corrfile.setf(ios::fixed);
      corrfile.precision(sigdig);
      for (int t=0; t<best.corr.get_size(); t++) {
	corrfile << nbpt(t) << "\t" << diam(t) << "\t" << best.corr(t) << "\t" << tcorr(t) << "\t" << (best.corr(t)-tcorr(t)) << endl;
      }
      if (best.obj==-MAXFLOAT) {
	corrfile << "Objective_function= Perfect_match" << endl; 
	logfile  << "Objective_function= Perfect_match" << endl << flush;
      }
      else {
        corrfile << "Objective_function= " << best.obj << endl;
	logfile  << "Objective_function= " << best.obj << endl << flush;
      }
      logfile << "Correlations_mismatch= ";
      for (int t=0; t<best.corr.get_size(); t++) {
	logfile << "\t" << (best.corr(t)-tcorr(t));
      }
      logfile << endl << flush;
      if (best.obj==-MAXFLOAT) {
	break;
      }
    }

    if (tic==0) {
      ifstream tempfile(paramfilename);
      Real newweightdist=1.;
      Real newweightnbpt=1.;
      Real newweightdecay=0.;
      Real newT=MAXFLOAT;
      if (tempfile) {
	tempfile >> newweightdist >> newweightnbpt >> newweightdecay >> newT;
	if (newT==MAXFLOAT) {newT=newweightnbpt; newweightdecay=0; newweightnbpt=1;}
	if (newweightdist!=weightdist || newweightnbpt!=weightnbpt || newweightdecay!=weightdecay || newT!=T) {
	  logfile << "New parameters read in: -wr=" << newweightdist << " -wp=" << newweightnbpt << " -wd=" << newweightdecay << " -T=" << newT << endl;
	}
	T=newT;
	if (newweightdist!=weightdist || newweightnbpt!=weightnbpt  || newweightdecay!=weightdecay) {
	  weightdist=newweightdist;
	  weightnbpt=newweightnbpt;
	  weightdecay=newweightdecay;
	  best.obj=MAXFLOAT;
	  for (int c=0; c<mc.get_size(); c++) {
	    mc(c).obj=calc_objective_func(mc(c).corr,tcorr,mysqstol,weightdist,weightnbpt,weightdecay,diam,nbpt);
	    if (mc(c).obj<best.obj) {
	      best=mc(c);
	      cc=c;
	    }
	  }
	  obj=best.obj;
	  best.obj=MAXFLOAT;
	}
      }
      if (file_exists("stopsqs")) {
	unlink("stopsqs");
	logfile << "Stopped" << endl; 
	break;
      }
      tic=maxtic;
    }
    int newcc=random(mc.get_size());
    int at1,at2;
    do {
      at1=random(mc(newcc).nbactive);
      int nbat=mc(newcc).typeend(at1)-mc(newcc).typebeg(at1);
      at2=mc(newcc).typebeg(at1) + ((1+random(nbat-1))%nbat);
    } while (mc(newcc).str.atom_type(at1)==mc(newcc).str.atom_type(at2));
    swap(&(mc(newcc).str.atom_type(at1)),&(mc(newcc).str.atom_type(at2)));
    Array<Real> newcorr;
    calc_all_corr(&newcorr,mc(newcc).str,eqclus,lat.cell,*pcorrfunc);
    Real newobj=calc_objective_func(newcorr,tcorr,mysqstol,weightdist,weightnbpt,weightdecay,diam,nbpt);
    //  cerr << newcc << " " << best.obj << " " << obj << " " << newobj << " ";
    //    for (int t=0; t<newcorr.get_size(); t++) {cerr << newcorr(t) << " ";}
    if (uniform01() < exp((obj-newobj)/T) ) {
      mc(newcc).corr=newcorr;
      mc(newcc).obj=newobj;
      cc=newcc;
      obj=newobj;
      //  cerr << "A";
    }
    else {
      swap(&(mc(newcc).str.atom_type(at1)),&(mc(newcc).str.atom_type(at2)));
      //  cerr << "r"; 
    }
    //  cerr << endl << flush;
    tic--;
  }
}
