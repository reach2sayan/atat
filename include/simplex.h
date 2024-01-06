#include "CVMOptimizer.h"

class Simplex {

private:
  CVMOptimizerDataHolder *data;
  double maximum;
  bool checkOptimality(const VectorXd &costs, const VectorXd &basic_vars);
  int findPivotColumn(const VectorXd costs, const VectorXd &basic_vars);
  int findPivotRow(const MatrixXd A, const VectorXd &b, int column);

  MatrixXd get_inverse_base_matrix(const MatrixXd &A,
                                   const VectorXd basic_vars);
  VectorXd get_simplex_mults(const MatrixXd &inverse_base,
                             const VectorXd &costs, const VectorXd &basic_vars);

  bool simplexAlgorithmCalculation();

public:
  Simplex(CVMOptimizerDataHolder *data);
  void calculate();
};
