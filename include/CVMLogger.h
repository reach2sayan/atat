#include <ifopt/ipopt_solver.h>

#include <fstream>
#include <iostream>
#include <memory>
#include <string>

class CVMOptimizerDataHolder;
class CVMLogger {
 public:
  virtual ~CVMLogger() {}
  virtual void log() = 0;
};

class CVMCorrelationsLogger : public CVMLogger {
 public:
  CVMCorrelationsLogger() : CVMCorrelationsLogger("cvmsteps.out", 3, nullptr) {}
  CVMCorrelationsLogger(const std::string &outfilename, const int sigdig,
			std::shared_ptr<CVMOptimizerDataHolder> data);
  void log() override;

 private:
  std::shared_ptr<CVMOptimizerDataHolder> data_;
  int sigdig_;
  std::ofstream outfile;
};

class CVMResultsLogger : public CVMLogger {
 public:
  CVMResultsLogger() : CVMResultsLogger("cvmresult.out", 3, nullptr) {}
  CVMResultsLogger(const std::string &outfilename, const int sigdig,
		   std::shared_ptr<CVMOptimizerDataHolder> data);
  void log() override;

 private:
  std::shared_ptr<CVMOptimizerDataHolder> data_;
  int sigdig_;
  std::ofstream outfile;
};

class CVMSolverLogger : public CVMLogger {
 public:
  CVMSolverLogger() : CVMSolverLogger("ipopt.log", nullptr) {}
  CVMSolverLogger(const std::string &outfilename,
		  std::shared_ptr<CVMOptimizerDataHolder> data);
  void log() override;

 private:
  std::shared_ptr<CVMOptimizerDataHolder> data_;
  std::ofstream solveroutfile;
};
