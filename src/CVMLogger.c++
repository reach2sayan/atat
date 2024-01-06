#include "CVMLogger.h"
#include "CVMDataHolder.h"

CVMFileLogger::CVMFileLogger(const std::string &outfilename,
                             std::shared_ptr<CVMOptimizerDataHolder> data)
    : data_(data) {
  outfile.open(outfilename);
  outfile << "Disordered Correlations : "
          << data_->cvminfo.disordered_correlation.transpose() << endl;
  outfile << "Ordered Correlations : "
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

void CVMFileLogger::log() {
  outfile << "T | " << setw(20) << std::right << " Ordered FE | " << setw(20)
          << std::right << " Disordered FE | " << setw(20) << std::right
          << " Optimised FE " << endl;
  outfile << "-----------------------------------------------------------------"
             "---------------"
          << endl;
  outfile << "T = " << data_->cvminfo.temperature.back() << setw(20)
          << std::right << data_->cvminfo.ordered_fe.back() << setw(20)
          << std::right << data_->cvminfo.disordered_fe.back() << setw(20)
          << std::right << data_->cvminfo.opt_fe.back() << endl;
  outfile << "T = " << data_->cvminfo.temperature.back() << setw(40)
          << "Optimised Correlation : "
          << data_->cvminfo.opt_correlations.back().transpose() << endl;
  outfile << "T = " << data_->cvminfo.temperature.back() << setw(20)
          << "Rhos : " << endl;
  data_->PrettyPrintRho(data_->cvminfo.opt_correlations.back(), outfile);
  outfile << "================================================================="
             "==============="
          << endl
          << endl;
}

CVMConsoleLogger::CVMConsoleLogger(std::shared_ptr<CVMOptimizerDataHolder> data)
    : data_(data) {
  cout << "Disordered Correlations : "
       << data_->cvminfo.disordered_correlation.transpose() << endl;
  cout << "Ordered Correlations : "
       << data_->cvminfo.ordered_correlation.transpose() << endl;
  cout << "T | " << setw(20) << std::right << " Ordered FE | " << std::right
       << setw(20) << " Disordered FE | " << std::right << setw(20)
       << " Optimised FE " << endl;
  cout << "--------------------------------------------------------------------"
          "------------"
       << endl;
  data_->addLogger(this);
}

void CVMConsoleLogger::log() {
  cout << data_->cvminfo.temperature.back() << setw(15) << std::right
       << data_->cvminfo.ordered_fe.back() << setw(15) << std::right
       << data_->cvminfo.disordered_fe.back() << setw(15) << std::right
       << data_->cvminfo.opt_fe.back() << setw(15) << endl;
}
