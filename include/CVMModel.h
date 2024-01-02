#ifndef __CVM_PROBLEM_DEFINITION_HPP__
#define __CVM_PROBLEM_DEFINITION_HPP__

#include <ifopt/constraint_set.h>
#include <ifopt/cost_term.h>
#include <ifopt/variable_set.h>
#include <vector>
#include "CVMDataHolder.h"
#include "linklist.h"
#include "xtalutil.h"

namespace ifopt {

using Eigen::VectorXd;
using Eigen::MatrixXd;

const std::string DEFAULT_CORRELATIONS_SET_NAME{"cvm-correlations"};
const std::string DEFAULT_RHO_CONST_NAME{"constraint-rho"};
const std::string DEFAULT_NORM_CONST_NAME{"constraint-norm"};
const std::string DEFAULT_ENERGY_NAME{"energy"};
const std::string DEFAULT_FREEENERGY_NAME{"free-energy"};

class CVMCorrelations : public VariableSet {

	public:
		//CVMCorrelations() : CVMCorrelations(DEFAULT_CORRELATIONS_SET_NAME, 0, 0) {};
		CVMCorrelations(const std::string& name, const int num_clusters, const int num_point_clusters);

		void SetVariables(const VectorXd& x) override { correlations = x ; }
		VectorXd GetValues() const override { return correlations; }
		VecBound GetBounds() const override;
		EIGEN_MAKE_ALIGNED_OPERATOR_NEW

	private:
		VectorXd correlations;
		int num_point_clusters_;
};

class CVMNormConstraints : public ConstraintSet {

	//CVMNormConstraints() : CVMNormConstraints(DEFAULT_NORM_CONST_NAME, 0, 0, VectorXd::Zero) {}
	public:
		CVMNormConstraints(const std::string& name, const VectorXd& disordered_corr, const VectorXd& ordered_corr);
		VectorXd GetValues() const override;
		VecBound GetBounds() const override;
		void FillJacobianBlock(std::string var_set, Jacobian& jac_block) const override {};
		EIGEN_MAKE_ALIGNED_OPERATOR_NEW

	private:
		VectorXd disordered_corr_;
		VectorXd ordered_corr_;
};

class CVMRhoConstraints : public ConstraintSet {

	public:
		//CVMRhoConstraints() : CVMRhoConstraints(DEFAULT_RHO_CONST_NAME, 0, 0, MatrixXd::Zero(0,0)) {}
		CVMRhoConstraints(const std::string& name, const int num_config, const int num_point_configs, const MatrixXd& vmatrix);

		VectorXd GetValues() const override;
		VecBound GetBounds() const override;
		void FillJacobianBlock(std::string var_set, Jacobian& jac_block) const override;
		EIGEN_MAKE_ALIGNED_OPERATOR_NEW

	private:
		MatrixXd vmatrix_;
		int num_point_configs_;
};

class CVMEnergy : public CostTerm {

	public: 
		CVMEnergy() : CVMEnergy(DEFAULT_ENERGY_NAME, VectorXd::Zero(0)) {}
		CVMEnergy(const std::string& name, const VectorXd& mult_eci) : CostTerm(name) , mult_eci_(mult_eci) {}

		double GetCost() const override;
		double GetCost(const VectorXd& corrs) const;
		void FillJacobianBlock(std::string var_set, Jacobian& jac_block) const override;
		EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	
	private:
		VectorXd mult_eci_;
};

class CVMFreeEnergy : public CostTerm {

	public:
		CVMFreeEnergy() : CVMFreeEnergy(DEFAULT_FREEENERGY_NAME, VectorXd::Zero(0), VectorXd::Zero(0), MatrixXd::Zero(0,0)) {}
		CVMFreeEnergy(const std::string& name, const VectorXd& mult_eci, const VectorXd& multconfig_kb, const MatrixXd& vmatrix) : CostTerm(name), mult_eci_(mult_eci), multconfig_kb_(multconfig_kb), vmatrix_(vmatrix) { T = 0; }

		double GetCost() const override;
		double GetCost(const VectorXd& corrs) const;
		void FillJacobianBlock(std::string var_set, Jacobian& jac_block) const override;

		double getT() const { return T; }
		void setT(const double T_) { T = T_; }
		EIGEN_MAKE_ALIGNED_OPERATOR_NEW

	private:
		VectorXd mult_eci_;
		VectorXd multconfig_kb_;
		MatrixXd vmatrix_;
		static constexpr double kB = 8.61733326e-05;
		double T;
};

}

#endif // __CVM_PROBLEM__DEFINITION_HPP__
