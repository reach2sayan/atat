#include "CVMDataHolder.h"

#include <ifopt/ipopt_solver.h>

#include "CVMLogger.h"
#include "parse.h"

ostream &operator<<(ostream &s, const MultiCluster &multiclus) {
  s << multiclus.clus.get_size() << endl;
  for (int i = 0; i < multiclus.clus.get_size(); i++) {
    s << multiclus.clus(i) << " " << multiclus.func(i) << " "
      << multiclus.site_type(i) << endl;
  }
  return s;
}

void CVMOptimizerDataHolder::log() const {
  for (const auto logger : loggers) logger->log();
}

void CVMOptimizerDataHolder::addInfo(const VectorXd &optcorr, double optfe,
				     double ordfe, double disordfe, double T) {
  double conv = GetConvFactor(
      UNITTYPE);  // evperatom_kJpermole / lat.atom_pos.get_size();
  cvminfo.opt_correlations.push_back(optcorr);
  cvminfo.opt_fe.push_back(optfe * conv);
  cvminfo.ordered_fe.push_back(ordfe * conv);
  cvminfo.disordered_fe.push_back(disordfe * conv);
  cvminfo.temperature.push_back(T);
}

double CVMOptimizerDataHolder::GetConvFactor(const UnitType type) const {
  double retval = 1.;
  switch (type) {
    case EVPERATOM:
      retval = 1. / lat.atom_pos.get_size();
      break;
    case KJPERMOLE:
      retval = evperatom_kJpermole / lat.atom_pos.get_size();
      break;
  }
  return retval;
}

void CVMOptimizerDataHolder::PrettyPrintRho(const VectorXd &corrs, ostream &out,
					    const int sigdig) {
  out.precision(sigdig);
  const VectorXd rhos = get_vmatrix() * corrs;
  int rhocount = 0;
  LinkedListIterator<LinkedList<MultiCluster>> iconfiglist(configlistlist);
  for (; iconfiglist; iconfiglist++) {
    LinkedListIterator<MultiCluster> iequiconfig(*iconfiglist);
    for (; iequiconfig; iequiconfig++) {
      out << rhos[rhocount] << " ";
      rhocount++;
    }
    out << endl;
  }
}

void CVMOptimizerDataHolder::saveClusterInformation(
    const string &clusterfname, const string &clustermultfname,
    const string &configfname, const string &configmultfname,
    const string &vmatrixfname, const string &kbfname) {
  constexpr int sigdig = 5;

  ofstream clusterfile(clusterfname);
  clusterfile.setf(ios::fixed);
  clusterfile.precision(sigdig);
  LinkedListIterator<MultiCluster> icluster(multicluslist);

  for (; icluster; icluster++) {
    int mult = calc_multiplicity(*icluster, spacegroup.cell,
				 spacegroup.point_op, spacegroup.trans);

    clusterfile << mult << endl;
    clusterfile << get_length_slow(icluster->clus) << endl;
    clusterfile << icluster->clus.get_size() << endl;
    for (int i = 0; i < icluster->clus.get_size(); i++)
      clusterfile << ((!axes) * (icluster->clus(i))) << " "
		  << icluster->site_type(i) << " " << icluster->func(i) << endl;

    clusterfile << endl;
  }

  ofstream clusmultfile(clustermultfname);
  clusmultfile << multicluslist.get_size() << endl;
  LinkedListIterator<MultiCluster> icluster2(multicluslist);
  for (; icluster2; icluster2++) {
    int mult = calc_multiplicity(*icluster2, spacegroup.cell,
				 spacegroup.point_op, spacegroup.trans);
    clusmultfile << mult << endl;
  }

  ofstream configfile(configfname);
  configfile.setf(ios::fixed);
  configfile.precision(sigdig);
  configfile << configlistlist;

  ofstream configmultfile(configmultfname);
  configmultfile.setf(ios::fixed);
  configmultfile.precision(sigdig);
  configmultfile << multlistlist;

  ofstream vmatfile(vmatrixfname);
  vmatfile.setf(ios::fixed);
  vmatfile.precision(sigdig);
  vmatfile << vmatlist;

  ofstream configkbfile(kbfname);
  configkbfile.setf(ios::fixed);
  configkbfile.precision(sigdig);
  configkbfile << kbcoef;
}

void CVMOptimizerDataHolder::PrintOptimizationVectorsAndMatrices() {
  const VectorXd &mult_eci = get_mults_eci();
  const VectorXd &multconfig_kb = get_multconfig_kb();
  const MatrixXd &vmatrix = get_vmatrix();

  cout << "multiplicity * eci" << endl << mult_eci << endl;
  cout << "config multiplicity * kb" << endl << multconfig_kb << endl;
  cout << "vmatrix" << endl << vmatrix << endl;

  cout << "ECIs" << endl;
  LinkedListIterator<double> ieci(ecilist);
  for (; ieci; ieci++) cout << *ieci << endl;

  cout << "No. of Clusters: " << mult_eci.size() << endl;
  cout << "No. of single-point clusters: " << get_num_point_clusters() << endl;
  cout << "No. of configs: " << multconfig_kb.size() << endl;
}

void CVMOptimizerDataHolder::find_subclusters(
    LinkedList<MultiCluster> *pcluslist, const MultiCluster &maxclus,
    const SpaceGroup &spacegroup) {
  Array<int> maxselect(maxclus.clus.get_size());
  for (int i = 0; i < maxselect.get_size(); i++)
    maxselect(i) = maxclus.site_type(i) > 1 ? 2 : 1;
  MultiDimIterator<Array<int>> select(maxselect);

  for (; select; select++) {
    Array<int> &sel = (Array<int> &)select;
    int nbpt = sum(sel);
    MultiCluster *pclus = new MultiCluster(nbpt);
    int j = 0;
    for (int i = 0; i < maxselect.get_size(); i++) {
      if (sel(i) == 1) {
	pclus->clus(j) = maxclus.clus(i);
	pclus->func(j) = 0;
	pclus->site_type(j) = maxclus.site_type(i);  // nb of species;
	j++;
      }
    }

    LinkedListIterator<MultiCluster> iclus(*pcluslist);
    for (; iclus; iclus++)
      if (equivalent_by_symmetry(*pclus, *iclus, spacegroup.cell,
				 spacegroup.point_op, spacegroup.trans))
	break;

    if (!iclus) {
      iclus.init(*pcluslist);
      for (; iclus; iclus++) {
	if (iclus->clus.get_size() > pclus->clus.get_size() ||
	    (iclus->clus.get_size() == pclus->clus.get_size() &&
	     get_length_slow(iclus->clus) > get_length_slow(pclus->clus)))
	  break;
      }
      pcluslist->add(pclus, iclus);
    } else {
      delete pclus;
    }
  }
}

void CVMOptimizerDataHolder::generate_functions_on_clusters(
    LinkedList<MultiCluster> *pmulticluslist,
    const LinkedList<MultiCluster> &cluslist, const SpaceGroup &spacegroup) {
  LinkedListIterator<MultiCluster> begclus(*pmulticluslist);
  LinkedListIterator<MultiCluster> iclus(cluslist);
  for (; iclus; iclus++) {
    Array<int> occup(iclus->clus.get_size());
    for (int i = 0; i < occup.get_size(); i++) {
      occup(i) = iclus->site_type(i) - 1;
    }
    MultiDimIterator<Array<int>> conf(occup);
    for (; conf; conf++) {
      MultiCluster *pclus = new MultiCluster;
      pclus->clus = iclus->clus;
      pclus->func = conf;
      pclus->site_type.resize(occup.get_size());
      for (int i = 0; i < occup.get_size(); i++) {
	pclus->site_type(i) = occup(i) - 1;
      }	 // func type: 0=binary, 1=ternary;
      LinkedListIterator<MultiCluster> jclus(begclus);
      for (; jclus; jclus++) {
	if (equivalent_by_symmetry(*pclus, *jclus, spacegroup.cell,
				   spacegroup.point_op, spacegroup.trans)) {
	  break;
	}
      }
      if (jclus) {
	delete pclus;
      } else {
	(*pmulticluslist) << pclus;
      }
    }
    while (begclus) {
      begclus++;
    }
  }
}

void CVMOptimizerDataHolder::generate_config_on_clusters(
    LinkedList<LinkedList<MultiCluster>> *pconfiglistlist,
    LinkedList<LinkedList<double>> *pmultlistlist,
    const LinkedList<MultiCluster> &cluslist, const SpaceGroup &spacegroup) {
  LinkedListIterator<MultiCluster> iclus(cluslist);
  for (; iclus; iclus++) {
    LinkedList<MultiCluster> *pconfiglist = new LinkedList<MultiCluster>;
    LinkedList<double> *pmultlist = new LinkedList<double>;
    SpaceGroup subgroup;
    find_clusters_symmetry(&subgroup, *iclus, spacegroup);

    Array<int> occup(iclus->clus.get_size());
    for (int i = 0; i < occup.get_size(); i++) {
      occup(i) = iclus->site_type(i);
    }
    MultiDimIterator<Array<int>> conf(occup);
    for (; conf; conf++) {
      MultiCluster *pclus = new MultiCluster;
      pclus->clus = iclus->clus;
      pclus->func = conf;
      pclus->site_type = occup;	 // nb of species;
      LinkedListIterator<MultiCluster> jclus(*pconfiglist);
      LinkedListIterator<double> jmult(*pmultlist);
      for (; jclus; jclus++, jmult++) {
	if (equivalent_by_symmetry(*pclus, *jclus, subgroup.point_op,
				   subgroup.trans)) {
	  break;
	}
      }
      if (jclus) {
	(*jmult) += 1.;
	delete pclus;
      } else {
	(*pconfiglist) << pclus;
	(*pmultlist) << new double(1.);
      }
    }
    (*pconfiglistlist) << pconfiglist;
    (*pmultlistlist) << pmultlist;
  }
}

void CVMOptimizerDataHolder::calc_v_matrix(
    LinkedList<Array2d<double>> *pvmatlist,
    const LinkedList<LinkedList<MultiCluster>> &configlistlist,
    const LinkedList<LinkedList<double>> &configmultlistlist,
    const LinkedList<MultiCluster> &multicluslist, const SpaceGroup &spacegroup,
    const CorrFuncTable &corrfunc) {
  LinkedListIterator<LinkedList<MultiCluster>> iconfiglist(configlistlist);
  LinkedListIterator<LinkedList<double>> iconfigmultlist(
      configmultlistlist);  //

  for (; iconfiglist; iconfiglist++, iconfigmultlist++) {
    Array2d<double> *pvmat =
	new Array2d<double>(iconfiglist->get_size(), multicluslist.get_size());
    Array2d<double> ximat(multicluslist.get_size(),
			  iconfiglist->get_size());  // debug

    LinkedListIterator<MultiCluster> iclus(multicluslist);
    for (int j = 0; iclus; iclus++, j++) {
      double prefac = 1.;
      {
	LinkedListIterator<MultiCluster> iconf(*iconfiglist);
	for (int p = 0; p < iconf->site_type.get_size(); p++) {
	  int at = which_atom(iclus->clus, iconf->clus(p));
	  if (at == -1) {
	    prefac *= (double)(iconf->site_type(p));
	  } else {
	    double sum2 = 0.;
	    for (int q = 0;
		 q < corrfunc(iclus->site_type(at))(iclus->func(at)).get_size();
		 q++) {
	      sum2 += sqr(corrfunc(iclus->site_type(at))(iclus->func(at))(q));
	    }
	    prefac *= sum2;
	  }
	}
      }
      // cout << prefac << " ";

      Array<MultiCluster> equiv_clusters;
      find_equivalent_clusters(&equiv_clusters, *iclus, spacegroup.cell,
			       spacegroup.point_op, spacegroup.trans);
      LinkedListIterator<MultiCluster> iconf(*iconfiglist);
      LinkedListIterator<double> iconfigmult(*iconfigmultlist);	 // debug
      for (int i = 0; iconf; iconf++, i++, iconfigmult++) {
	int mult;
	double corrmult = calc_correlation_atomcluster(
	    iconf->clus, iconf->func,
	    const_cast<Array<MultiCluster> &>(equiv_clusters), spacegroup.cell,
	    corrfunc, &mult);
	(*pvmat)(i, j) = corrmult / prefac;
	ximat(j, i) = (mult == 0 ? 0 : corrmult / (double)mult) *
		      (*iconfigmult);  // debug
      }
    }
    (*pvmatlist) << pvmat;
  }
}

void CVMOptimizerDataHolder::calc_kikuchi_barker(
    Array<double> *pkbcoef, const LinkedList<MultiCluster> &cluslist,
    const SpaceGroup &spacegroup) {
  Array<double> &kbcoef = *pkbcoef;
  kbcoef.resize(cluslist.get_size());
  kbcoef(0) = 0.;
  for (int curclus = kbcoef.get_size() - 1; curclus > 0; curclus--) {
    kbcoef(curclus) = 1.;
    LinkedListIterator<MultiCluster> iclus =
	LinkedListIterator<MultiCluster>(cluslist) + curclus;
    Array<MultiCluster> candi_clus;
    Array<int> which_clus;
    find_clusters_overlapping_site(&candi_clus, &which_clus, iclus->clus(0),
				   spacegroup, cluslist);
    for (int c = 0; c < candi_clus.get_size(); c++) {
      if (candi_clus(c).clus.get_size() > iclus->clus.get_size()) {
	int i;
	for (i = 0; i < iclus->clus.get_size(); i++) {
	  if (which_atom(candi_clus(c).clus, iclus->clus(i)) == -1) break;
	}
	if (i == iclus->clus.get_size()) {
	  kbcoef(curclus) -= kbcoef(which_clus(c));
	}
      }
    }
  }
}

MatrixXd CVMOptimizerDataHolder::get_vmatrix() const {
  LinkedListIterator<Array2d<double>> ivmat(vmatlist);
  int num_clusters = get_num_clusters();
  int num_configs = get_num_configs();
  MatrixXd vmatrix = MatrixXd::Zero(num_configs, num_clusters);
  int row_position = 0;
  for (; ivmat; ivmat++) {
    iVector2d size = (*ivmat).get_size();
    assert(size[1] == num_clusters);
    int num_configs_in_cluster = size[0];
    MatrixXd current_vmat =
	MatrixXd::Zero(num_configs_in_cluster, num_clusters);
    for (int i = 0; i < num_configs_in_cluster; i++)
      for (int j = 0; j < num_clusters; j++)
	current_vmat(i, j) = (*ivmat)(i, j);

    vmatrix.block(row_position, 0, num_configs_in_cluster, num_clusters) =
	current_vmat;
    row_position += num_configs_in_cluster;
  }

  return vmatrix;
}

VectorXd CVMOptimizerDataHolder::get_mults_eci() const {
  int num_clusters = get_num_clusters();
  VectorXd mult_eci = VectorXd::Zero(num_clusters);
  LinkedListIterator<MultiCluster> icluster(multicluslist);
  LinkedListIterator<Real> ieci(ecilist);

  int count = 0;
  for (; icluster && ieci; icluster++, ieci++, count++) {
    mult_eci(count) = calc_multiplicity(*icluster, spacegroup.cell,
					spacegroup.point_op, spacegroup.trans) *
		      (*ieci);
  }

  return mult_eci;
}

int CVMOptimizerDataHolder::get_num_point_clusters() const {
  int num_point_clusters = 0;
  LinkedListIterator<MultiCluster> icluster(multicluslist);
  for (; icluster; icluster++) {
    if (icluster->clus.get_size() == 1) num_point_clusters++;
  }

  return num_point_clusters;
}

VectorXd CVMOptimizerDataHolder::get_multconfig_kb() const {
  int num_configs = get_num_configs();
  VectorXd multconfig_kb = VectorXd::Zero(num_configs);
  LinkedListIterator<LinkedList<double>> iconfigmult(multlistlist);
  int num_mults = 0, kb_count = 0, mult_count = 0;
  for (; iconfigmult; iconfigmult++, kb_count++) {
    num_mults += (*iconfigmult).get_size();
    LinkedListIterator<double> imult(*iconfigmult);
    for (; imult; imult++, mult_count++)
      multconfig_kb(mult_count) = (*imult) * kbcoef(kb_count);
  }

  return multconfig_kb;
}

int CVMOptimizerDataHolder::get_num_configs() const {
  LinkedListIterator<LinkedList<double>> iconfig(multlistlist);
  int num_configs = 0;
  for (; iconfig; iconfig++) num_configs += (*iconfig).get_size();

  return num_configs;
}

CVMOptimizerDataHolder::CVMOptimizerDataHolder(const string &maxclusfilename,
					       const string &latfilename,
					       const string &ecifilename,
					       const string &strfilename) {
  Array<Arrayint> labellookup;
  Array<AutoString> label;
  constexpr const char *corrfunc_label = "trigo";

  ifstream latfile(latfilename);
  if (!latfile) ERRORQUIT("Unable to open lattice file");
  parse_lattice_file(&lat.cell, &lat.atom_pos, &lat.atom_type, &labellookup,
		     &label, latfile, &axes);
  wrap_inside_cell(&lat.atom_pos, lat.atom_pos, lat.cell);

  rMatrix3d icell = !lat.cell;
  spacegroup.cell = lat.cell;
  find_spacegroup(&spacegroup.point_op, &spacegroup.trans, lat.cell,
		  lat.atom_pos, lat.atom_type);
  if (contains_pure_translations(spacegroup.point_op, spacegroup.trans)) {
    cerr << "Warning: unit cell is not primitive." << endl;
  }
  LinkedList<MultiCluster> cluslist;

  ifstream maxclusfile(maxclusfilename);
  if (!maxclusfile) ERRORQUIT("Unable to open maximal cluster file");

  while (true) {
    LinkedList<rVector3d> maxcluslist;
    while (true) {
      rVector3d v(MAXFLOAT, MAXFLOAT, MAXFLOAT);
      maxclusfile >> v;
      if (maxclusfile.eof() || v(0) == MAXFLOAT || !maxclusfile) break;
      maxcluslist << new rVector3d(axes * v);
    }
    skip_to_next_structure(maxclusfile);
    MultiCluster maxclus;
    LinkedList_to_Array(&maxclus.clus, maxcluslist);

    if (maxclus.clus.get_size() == 0) break;

    maxclus.site_type.resize(maxclus.clus.get_size());
    for (int i = 0; i < maxclus.clus.get_size(); i++)
      maxclus.site_type(i) =
	  labellookup(
	      lat.atom_type(which_atom(lat.atom_pos, maxclus.clus(i), icell)))
	      .get_size();
    find_subclusters(&cluslist, maxclus, spacegroup);
  }

  generate_functions_on_clusters(&multicluslist, cluslist, spacegroup);
  int num_clusters = multicluslist.get_size();

  Array<Real> eci;
  ifstream ecifile(ecifilename);
  int i = 0;
  if (!ecifile) ERRORQUIT("Unable to open ECI file.");
  while (skip_delim(ecifile)) {
    Real e;
    ecifile >> e;
    if (i == 0) {
      ++i;
      continue;
    }
    ecilist << new Real(e);
  }
  if (multicluslist.get_size() != ecilist.get_size())
    ERRORQUIT("Number of ECI does not match number of clusters.");

  ecifile.close();
  maxclusfile.close();
  latfile.close();

  generate_config_on_clusters(&configlistlist, &multlistlist, cluslist,
			      spacegroup);
  calc_kikuchi_barker(&kbcoef, cluslist, spacegroup);

  if (!check_plug_in(CorrFuncTable(), corrfunc_label)) ERRORQUIT("Aborting");
  CorrFuncTable *pcorrfunc_ =
      GenericPlugIn<CorrFuncTable>::create(corrfunc_label);
  pcorrfunc_->init_from_site_type_list(labellookup);
  pcorrfunc = pcorrfunc_;

  calc_v_matrix(&vmatlist, configlistlist, multlistlist, multicluslist,
		spacegroup, *pcorrfunc);
  get_structure_from_file(str, strfilename, lat, label, labellookup);

  LinkedListIterator<MultiCluster> icluster(multicluslist);
  for (; icluster; icluster++) {
    Array<MultiCluster> *pmulticlus = new Array<MultiCluster>;
    find_equivalent_clusters(pmulticlus, *icluster, spacegroup.cell,
			     spacegroup.point_op, spacegroup.trans);
    equivcluslist << pmulticlus;
  }

  cvminfo.disordered_correlation = VectorXd::Zero(multicluslist.get_size());
  cvminfo.ordered_correlation = VectorXd::Zero(multicluslist.get_size());

  UNITTYPE = KJPERMOLE;
}

void CVMOptimizerDataHolder::get_structure_from_file(
    Structure &str, const string &strfilename, const Structure &lat,
    const Array<AutoString> &label, const Array<Arrayint> &labellookup) {
  ifstream strfile(strfilename);
  if (!strfile) ERRORQUIT("Unable to open structure file");
  int strnum = 1;
  while (!strfile.eof()) {
    parse_structure_file(&str.cell, &str.atom_pos, &str.atom_type, label,
			 strfile);
    skip_to_next_structure(strfile);
    wrap_inside_cell(&str.atom_pos, str.atom_pos, str.cell);

    rMatrix3d supercell = (!lat.cell) * str.cell;
    rMatrix3d rounded_supercell = to_real(to_int(supercell));
    rMatrix3d transfo = lat.cell * rounded_supercell * (!str.cell);
    str.cell = transfo * str.cell;
    for (int i = 0; i < str.atom_pos.get_size(); i++)
      str.atom_pos(i) = transfo * str.atom_pos(i);
    fix_atom_type(&str, lat, labellookup, 0);
  }
}

VectorXd CVMOptimizerDataHolder::GetSampleCorrelation() const {
  static Structure mystructure = str;
  VectorXd corrs = VectorXd::Zero(get_num_clusters());
  int ieci = 0;

  LinkedListIterator<Array<MultiCluster>> icluster(equivcluslist);
  mystructure.atom_pos.shuffle();
  for (; icluster; icluster++, ieci++) {
    double rho =
	calc_correlation(mystructure, *icluster, spacegroup.cell, *pcorrfunc);
    corrs(ieci) = rho;
  }
  corrs(0) = 1.;
  return corrs;
}

VectorXd CVMOptimizerDataHolder::GetDisorderedCorrelation() const {
  VectorXd corrs = VectorXd::Zero(get_num_clusters());
  LinkedList<Real> pointcorr;
  LinkedListIterator<Array<MultiCluster>> ipcluster(equivcluslist);

  int ieci = 0;
  for (; ipcluster; ipcluster++) {
    if ((*ipcluster)(0).clus.get_size() > 1) break;
    if ((*ipcluster)(0).clus.get_size() == 1) {
      pointcorr << new Real(
	  calc_correlation(str, *ipcluster, spacegroup.cell, *pcorrfunc));
    }
  }

  ieci = 0;
  LinkedListIterator<MultiCluster> icluster(multicluslist);
  LinkedListIterator<Array<MultiCluster>> ieqcluster(equivcluslist);
  for (; icluster; icluster++, ieqcluster++, ieci++) {
    double rho = 1.;
    for (int s = 0; s < icluster->clus.get_size(); s++) {
      LinkedListIterator<MultiCluster> ipcluster(multicluslist);
      while (ipcluster->clus.get_size() == 0) ipcluster++;
      LinkedListIterator<Real> ipcorr(pointcorr);
      while (true) {
	if (ipcluster->site_type(0) == icluster->site_type(s) &&
	    ipcluster->func(0) == icluster->func(s) &&
	    equivalent_by_symmetry(ipcluster->clus(0), icluster->clus(s),
				   spacegroup.cell, spacegroup.point_op,
				   spacegroup.trans))
	  break;
	ipcluster++;
	ipcorr++;
      }
      rho *= (*ipcorr);
      corrs(ieci) = rho;
    }
  }
  corrs(0) = 1.;
  return corrs;
}
