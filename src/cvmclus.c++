#include <ifopt/ipopt_solver.h>
#include <ifopt/problem.h>

#include <algorithm>
#include <fstream>
#include <memory>
#include <ostream>
#include <random>

#include "CVMLogger.h"
#include "CVMModel.h"
#include "curvefit.h"
#include "getvalue.h"
#include "parse.h"
#include "thermofunctions.h"
#include "version.h"

extern const char *helpstring;
using namespace ifopt;
typedef Problem CVMModel;
typedef std::vector<std::pair<VectorXd, double>> IterMapType;

int main(int argc, char *argv[]) {
  // parsing command line. See getvalue.hh for details;
  const char *latfilename = "lat.in";
  const char *maxclusfilename = "maxclus.in";
  int sigdig = 5;
  const char *corrfunc_label = "trigo";
  int dohelp = 0;
  int dummy = 0;
  double Tmin = .0;
  double Tmax = 2000.0;
  double Tinc = 500.0;
  int nlocal = 10;
  int verbosity = 5;

  AskStruct options[] = {
      {"",
       "Cluster Variation Method CLUSter generator " MAPS_VERSION
       ", by Axel van de Walle",
       TITLEVAL, NULL},
      {"-h", "Display more help", BOOLVAL, &dohelp},
      {"-l", "Input file defining the lattice (Default: lat.in)", STRINGVAL,
       &latfilename},
      {"-m", "Maximal cluster (Default: maxclus.in)", STRINGVAL,
       &maxclusfilename},
      {"-sig", "Number of significant digits printed (Default: 5)", INTVAL,
       &sigdig},
      {"-crf", "Select correlation functions (default: trigo)", STRINGVAL,
       &corrfunc_label},
      {"-Tmin", "Starting Temperature (default 300 K)", INTVAL, &Tmin},
      {"-Tmax", "Starting Temperature (default 300 K)", INTVAL, &Tmax},
      {"-Tinc", "Starting Temperature (default 300 K)", INTVAL, &Tinc},
      {"-nlocal", "No. of local searches for global minima search (default 10)",
       INTVAL, &nlocal},
      {"-v", "verbosity level of IPOPT solver 1--20 (default 5)", INTVAL,
       &verbosity},
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

  auto cvmdata =
      std::make_shared<CVMOptimizerDataHolder>(maxclusfilename, latfilename);
  cvmdata->saveClusterInformation();

  VectorXd disordered_correlation = cvmdata->GetDisorderedCorrelation();
  int num_clusters = cvmdata->get_num_clusters();
  int num_config = cvmdata->get_num_configs();
  int num_point_clusters = cvmdata->get_num_point_clusters();

  auto rhoConstraintSet = std::make_shared<CVMRhoConstraints>(
      "constraint-rho", num_config, num_point_clusters + 1,
      cvmdata->get_vmatrix());
  auto EnergyCost =
      std::make_shared<CVMEnergy>("energy", cvmdata->get_mults_eci());
  auto FreeEnergyCost = std::make_shared<CVMFreeEnergy>(
      "free-energy", cvmdata->get_mults_eci(), cvmdata->get_multconfig_kb(),
      cvmdata->get_vmatrix());
  FreeEnergyCost->setT(Tmin);
  auto correlationSet = std::make_shared<CVMCorrelations>(
      "cvm-correlations", num_clusters, num_point_clusters);
  correlationSet->SetVariables(disordered_correlation);

  // Solve for ordered state
  CVMModel OrderedStateModel;
  OrderedStateModel.AddVariableSet(correlationSet);
  OrderedStateModel.AddConstraintSet(rhoConstraintSet);
  OrderedStateModel.AddCostSet(EnergyCost);

  IpoptSolver ipopt_ordered, ipopt_opt;

  ipopt_ordered.SetOption("linear_solver", "mumps");
  ipopt_ordered.SetOption("jacobian_approximation", "exact");
  ipopt_ordered.SetOption("print_level", verbosity);
  ipopt_ordered.SetOption("derivative_test", "first-order");
  ipopt_ordered.Solve(OrderedStateModel);

  VectorXd orderedcorr = OrderedStateModel.GetOptVariables()->GetValues();

  // Solve for optimised state
  CVMModel OptimisedModel;
  OptimisedModel.AddVariableSet(correlationSet);
  OptimisedModel.AddConstraintSet(rhoConstraintSet);
  OptimisedModel.AddCostSet(FreeEnergyCost);

  auto normConstraintSet = std::make_shared<CVMNormConstraints>(
      "constraint-norm", disordered_correlation, orderedcorr);
  OptimisedModel.AddConstraintSet(normConstraintSet);

  ipopt_opt.SetOption("linear_solver", "mumps");
  ipopt_opt.SetOption("jacobian_approximation", "exact");
  ipopt_opt.SetOption("derivative_test", "first-order");
  ipopt_opt.SetOption("print_level", verbosity);

  cvmdata->cvminfo.ordered_correlation = orderedcorr;
  cvmdata->cvminfo.disordered_correlation = disordered_correlation;
  CVMResultsLogger energyresults("cvmresult.csv", sigdig, cvmdata);
  CVMCorrelationsLogger corrresults("cvmsteps.out", sigdig, cvmdata);

  for (int T = Tmin; T <= Tmax; T += Tinc) {
    FreeEnergyCost->setT(static_cast<double>(T));
    IterMapType mymap;

    for (int i = 0; i < nlocal; i++) {
      correlationSet->SetVariables(cvmdata->GetSampleCorrelation());
      ipopt_opt.Solve(OptimisedModel);

      VectorXd opttemp = OptimisedModel.GetOptVariables()->GetValues();
      double fetemp = FreeEnergyCost->GetCost(opttemp);
      mymap.push_back(std::make_pair(std::move(opttemp), fetemp));
    }

    auto minresult = std::min_element(
	mymap.begin(), mymap.end(), [](const auto &first, const auto &second) {
	  return first.second < second.second;
	});
    cvmdata->addInfo(minresult->first, minresult->second,
		     FreeEnergyCost->GetCost(orderedcorr),
		     FreeEnergyCost->GetCost(disordered_correlation), T);
    cvmdata->log();
  }

  Array<double> initvalues({0.5, -0.5, 0.1});
  vector<double> correction;
  std::transform(cvmdata->cvminfo.opt_fe.begin(), cvmdata->cvminfo.opt_fe.end(),
		 cvmdata->cvminfo.disordered_fe.begin(),
		 std::back_inserter(correction), std::minus<double>());
  Array<double> ys = correction;
  Array<double> xs = cvmdata->cvminfo.temperature;

  Array<double> sroparams =
      curve_fit(sroCorrectionFunction, initvalues, xs, ys);
  std::ofstream fsroparams("sro_params.out");
  fsroparams << sroparams;
}
