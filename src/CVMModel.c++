#include "CVMModel.h"
#include <limits>
using namespace ifopt;

CVMCorrelations::CVMCorrelations(const std::string &name,
                                 const int num_clusters,
                                 const int num_point_clusters)
    : VariableSet(num_clusters, name) {

  correlations = VectorXd::Zero(num_clusters);
  correlations[0] = 1.;
  num_point_clusters_ = num_point_clusters;
}

CVMCorrelations::VecBound CVMCorrelations::GetBounds() const {

  const int num_rows = GetRows();
  VecBound CorrelationBounds(num_rows);
  for (int i = 0; i < num_rows; i++) {
    if (i <= num_point_clusters_)
      CorrelationBounds[i] = Bounds(correlations[i], correlations[i]);
    else
      CorrelationBounds[i] = Bounds(-1., +1.);
  }

  return CorrelationBounds;
}

CVMRhoConstraints::CVMRhoConstraints(const std::string &name,
                                     const int num_config,
                                     const int num_point_configs,
                                     const MatrixXd &vmatrix)
    : ConstraintSet(num_config, name) {

  vmatrix_ = vmatrix;
  num_point_configs_ = num_point_configs;
}

CVMRhoConstraints::VecBound CVMRhoConstraints::GetBounds() const {

  const int num_rows = GetRows();
  VecBound RhoBounds(num_rows);
  RhoBounds[0] = Bounds(1., 1.);
  for (int i = 1; i < num_rows; i++) {
    RhoBounds[i] = Bounds(0., 1);
  }
  return RhoBounds;
}

CVMNormConstraints::VectorXd CVMNormConstraints::GetValues() const {

  VectorXd g(GetRows());
  VectorXd correlations =
      GetVariables()->GetComponent(DEFAULT_CORRELATIONS_SET_NAME)->GetValues();
  g(0) = (disordered_corr_ - ordered_corr_).squaredNorm() / 2 -
         (correlations - disordered_corr_).squaredNorm();
  return g;
}

CVMNormConstraints::CVMNormConstraints(const std::string &name,
                                       const VectorXd &disordered_corr,
                                       const VectorXd &ordered_corr)
    : ConstraintSet(1, name) {

  disordered_corr_ = disordered_corr;
  ordered_corr_ = ordered_corr;
}

CVMNormConstraints::VecBound CVMNormConstraints::GetBounds() const {

  const int num_rows = GetRows();
  VecBound NormBounds(num_rows);
  NormBounds.at(0) = Bounds(0., 0.);

  return NormBounds;
}

void CVMRhoConstraints::FillJacobianBlock(std::string var_set,
                                          Jacobian &jac_block) const {

  if (var_set == DEFAULT_CORRELATIONS_SET_NAME)
    jac_block = vmatrix_.sparseView();
}

VectorXd CVMRhoConstraints::GetValues() const {

  VectorXd g(GetRows());
  VectorXd correlations =
      GetVariables()->GetComponent(DEFAULT_CORRELATIONS_SET_NAME)->GetValues();
  return vmatrix_ * correlations;
}

double CVMEnergy::GetCost() const {
  VectorXd corrs =
      GetVariables()->GetComponent(DEFAULT_CORRELATIONS_SET_NAME)->GetValues();
  return GetCost(corrs);
  ;
};

double CVMEnergy::GetCost(const VectorXd &corrs) const {
  return mult_eci_.dot(corrs);
};

double CVMFreeEnergy::GetCost() const {
  VectorXd corrs =
      GetVariables()->GetComponent(DEFAULT_CORRELATIONS_SET_NAME)->GetValues();
  return GetCost(corrs);
}

double CVMFreeEnergy::GetCost(const VectorXd &corrs) const {
  auto rhologrho = [](double x) {
    return x * std::log(std::abs(x) + std::numeric_limits<double>::epsilon());
  };
  return mult_eci_.dot(corrs) +
         T * kB * multconfig_kb_.dot((vmatrix_ * corrs).unaryExpr(rhologrho));
};

void CVMEnergy::FillJacobianBlock(std::string var_set,
                                  Jacobian &jac_block) const {

  if (var_set == DEFAULT_CORRELATIONS_SET_NAME) {
    for (int i = 0; i < mult_eci_.size(); i++)
      jac_block.coeffRef(0, i) = mult_eci_[i];
  }
}

void CVMFreeEnergy::FillJacobianBlock(std::string var_set,
                                      Jacobian &jac_block) const {

  VectorXd corrs =
      GetVariables()->GetComponent(DEFAULT_CORRELATIONS_SET_NAME)->GetValues();
  auto onepluslog = [](double x) {
    return 1 + std::log(std::abs(x) + std::numeric_limits<double>::epsilon());
  };
  auto mult_T_kB = [T_ = getT()](double x) { return T_ * kB * x; };
  auto cwiseMult = [](double x, double y) { return x * y; };

  VectorXd jac =
      mult_eci_ + (vmatrix_.transpose() *
                   multconfig_kb_.binaryExpr(
                       (vmatrix_ * corrs).unaryExpr(onepluslog), cwiseMult))
                      .unaryExpr(mult_T_kB);
  if (var_set == DEFAULT_CORRELATIONS_SET_NAME) {
    for (int i = 0; i < mult_eci_.size(); i++)
      jac_block.coeffRef(0, i) = jac[i];
  }
}
