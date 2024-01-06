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
  outfile << "Disordered Correlations : "
	  << data_->cvminfo.disordered_correlation.transpose() << endl;
  outfile << "Ordered Correlations : "
	  << data_->cvminfo.ordered_correlation.transpose() << endl;
  outfile << "T | " << setw(20) << std::right << " Ordered FE | " << std::right
	  << setw(20) << " Disordered FE | " << std::right << setw(20)
	  << " Optimised FE " << endl;
  outfile << "----------------------------------------------------------------"
	  << endl;
  data_->addLogger(this);
}

void CVMResultsLogger::log() {
  outfile << std::right << data_->cvminfo.temperature.back() << setw(15)
	  << data_->cvminfo.ordered_fe.back() << setw(15)
	  << data_->cvminfo.disordered_fe.back() << setw(15)
	  << data_->cvminfo.opt_fe.back() << endl;
}
