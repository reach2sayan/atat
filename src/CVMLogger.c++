#include "CVMLogger.h"

#include "CVMDataHolder.h"

CVMCorrelationsLogger::CVMCorrelationsLogger(
    const std::string &outfilename, const int sigdig,
    std::shared_ptr<CVMOptimizerDataHolder> data)
    : data_(data) {
  sigdig_ = sigdig;
  outfile.open(outfilename);
  outfile << "H: Disordered Correlations : "
	  << data_->cvminfo.disordered_correlation.transpose() << endl;
  outfile << "H: Ordered Correlations : "
	  << data_->cvminfo.ordered_correlation.transpose() << endl;
  data_->addLogger(this);
}

CVMSolverLogger::CVMSolverLogger(const std::string &outfilename,
				 std::shared_ptr<CVMOptimizerDataHolder> data)
    : data_(data) {
  solveroutfile.open(outfilename);
  data_->addLogger(this);
}

void CVMSolverLogger::log() {}

void CVMCorrelationsLogger::log() {
  outfile.precision(sigdig_);
  outfile << "T: T = " << data_->cvminfo.temperature.back()
	  << " Optimised Correlation : "
	  << data_->cvminfo.opt_correlations.back().transpose() << endl;
  outfile << "T: T = " << data_->cvminfo.temperature.back()
	  << " Rhos : " << endl;
  data_->PrettyPrintRho(data_->cvminfo.opt_correlations.back(), outfile,
			sigdig_);
  outfile << "================================================================="
	     "==============="
	  << endl
	  << endl;
}

CVMResultsLogger::CVMResultsLogger(const std::string &outfilename,
				   const int sigdig,
				   std::shared_ptr<CVMOptimizerDataHolder> data)
    : data_(data) {
  sigdig_ = sigdig;
  outfile.open(outfilename);
  outfile << "T,Ordered FE,Disordered FE,Optimised FE" << endl;
  data_->addLogger(this);
}

void CVMResultsLogger::log() {
  outfile << data_->cvminfo.temperature.back() << ","
	  << data_->cvminfo.ordered_fe.back() << ","
	  << data_->cvminfo.disordered_fe.back() << ","
	  << data_->cvminfo.opt_fe.back() << endl;
}
