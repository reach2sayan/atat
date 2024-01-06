#include <fstream>
#include <ifopt/ipopt_solver.h>
#include <iostream>
#include <memory>
#include <string>

class CVMOptimizerDataHolder;
class CVMLogger {

public:
  virtual ~CVMLogger() {}
  virtual void log() = 0;
};

class CVMFileLogger : public CVMLogger {

public:
  CVMFileLogger() : CVMFileLogger("cvm.out", nullptr) {}
  CVMFileLogger(const std::string &outfilename,
                std::shared_ptr<CVMOptimizerDataHolder> data);
  void log() override;

private:
  std::shared_ptr<CVMOptimizerDataHolder> data_;
  std::ofstream outfile;
};

class CVMConsoleLogger : public CVMLogger {

public:
  CVMConsoleLogger() : CVMConsoleLogger(nullptr) {}
  CVMConsoleLogger(std::shared_ptr<CVMOptimizerDataHolder> data);
  void log() override;

private:
  std::shared_ptr<CVMOptimizerDataHolder> data_;
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
