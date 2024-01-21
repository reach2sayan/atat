#ifndef __CVMDATAHOLDER_HPP__
#define __CVMDATAHOLDER_HPP__

#include <Eigen/Dense>
#include <list>
#include <memory>
#include <vector>

#include "clus_str.h"
#include "findsym.h"
#include "linklist.h"
#include "xtalutil.h"

using VectorXd = Eigen::Matrix<double, Eigen::Dynamic, 1>;
using MatrixXd = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>;
constexpr double evperatom_Jpermole = 96491.5666370759;
enum UnitType { EVPERATOM, JPERMOLE };

template <class T>
T sum(const Array<T> &a) {
  T s = 0;
  for (int i = 0; i < a.get_size(); i++) s += a(i);
  return s;
}

struct CVMInfoPacket {
  VectorXd disordered_correlation;
  VectorXd ordered_correlation;
  std::vector<VectorXd, Eigen::aligned_allocator<VectorXd>> opt_correlations;
  std::vector<double> opt_fe;
  std::vector<double> ordered_fe;
  std::vector<double> disordered_fe;
  std::vector<double> temperature;
};

struct CVMTemperatureFunction {
  Array2d<Real> mat;
  Array<Real> vect;
  std::function<double(double T)> thermoFunc;
};

struct CVMSolverInfoPacket {};

class CVMLogger;
class CVMOptimizerDataHolder {
 public:  // logging
  void addLogger(CVMLogger *logger) { loggers.push_back(logger); }
  void removeLogger(CVMLogger *logger) { loggers.remove(logger); }
  void log() const;
  CVMInfoPacket cvminfo;
  CVMSolverInfoPacket solverinfo;
  void addInfo(const VectorXd &optcorr, double optfe, double ordfe,
	       double disordfe, double T);

 private:
  std::list<CVMLogger *> loggers;

 public:
  CVMOptimizerDataHolder(const string &maxclusfilename = "maxclus.in",
			 const string &latfilename = "lat.in",
			 const string &ecifilename = "eci.out",
			 const string &strfilename = "str.in");

  static void find_subclusters(LinkedList<MultiCluster> *pcluslist,
			       const MultiCluster &maxclus,
			       const SpaceGroup &spacegroup);
  static void generate_functions_on_clusters(
      LinkedList<MultiCluster> *pmulticluslist,
      const LinkedList<MultiCluster> &cluslist, const SpaceGroup &spacegroup);
  static void generate_config_on_clusters(
      LinkedList<LinkedList<MultiCluster>> *pconfiglistlist,
      LinkedList<LinkedList<double>> *pmultlistlist,
      const LinkedList<MultiCluster> &cluslist, const SpaceGroup &spacegroup);
  static void calc_v_matrix(
      LinkedList<Array2d<double>> *pvmatlist,
      const LinkedList<LinkedList<MultiCluster>> &configlistlist,
      const LinkedList<LinkedList<double>> &configmultlistlist,
      const LinkedList<MultiCluster> &multicluslist,
      const SpaceGroup &spacegroup, const CorrFuncTable &corrfunc);
  static void calc_kikuchi_barker(Array<double> *pkbcoef,
				  const LinkedList<MultiCluster> &cluslist,
				  const SpaceGroup &spacegroup);
  static void get_structure_from_file(Structure &str, const string &strfilename,
				      const Structure &lattice,
				      const Array<AutoString> &label,
				      const Array<Arrayint> &labellookup);

  VectorXd GetDisorderedCorrelation() const;
  VectorXd GetSampleCorrelation() const;
  double GetConvFactor(const UnitType type) const;

  void saveClusterInformation(const string &clusterfname = "clusters.out",
			      const string &clustermultfname = "clusmult.out",
			      const string &configfname = "config.out",
			      const string &configmultfname = "configmult.out",
			      const string &vmatrixfname = "vmat.out",
			      const string &kbfname = "configkb.out");

  void PrintOptimizationVectorsAndMatrices();
  void PrettyPrintRho(const VectorXd &corrs, ostream &out = std::cout,
		      const int sigdig = 3);
  int get_num_clusters() const { return ecilist.get_size(); }
  int get_num_configs() const;
  int get_num_point_clusters() const;
  VectorXd get_mults_eci() const;
  VectorXd get_multconfig_kb() const;
  MatrixXd get_vmatrix() const;
  int get_num_lattice_atoms() const { return lat.atom_pos.get_size(); }

  LinkedList<MultiCluster> GetClusterList() const { return multicluslist; }
  Array<double> GetMultList() const { return multlist; }
  LinkedList<LinkedList<MultiCluster>> GetConfigList() const {
    return configlistlist;
  }
  LinkedList<LinkedList<double>> GetConfigMultList() const {
    return multlistlist;
  }
  LinkedList<Array2d<double>> GetVmatrixList() const { return vmatlist; }
  LinkedList<Array<MultiCluster>> GetEquivClusterList() const {
    return equivcluslist;
  }
  LinkedList<double> GeetECIList() const { return ecilist; }

 private:
  LinkedList<MultiCluster> multicluslist;
  Array<double> multlist;
  LinkedList<LinkedList<MultiCluster>> configlistlist;
  LinkedList<LinkedList<double>> multlistlist;
  LinkedList<Array2d<double>> vmatlist;
  LinkedList<Array<MultiCluster>> equivcluslist;
  LinkedList<double> ecilist;
  Array<double> kbcoef;
  SpaceGroup spacegroup;
  Structure lat;
  Structure str;
  rMatrix3d axes;
  CorrFuncTable *pcorrfunc;

  UnitType UNITTYPE;
};
#endif	// __CVMDATAHOLDER_HPP__
