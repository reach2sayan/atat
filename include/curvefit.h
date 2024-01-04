#ifndef __ATAT_CURVE_FIT_HPP__
#define __ATAT_CURVE_FIT_HPP__

#include <gsl/gsl_vector.h>
#include <gsl/gsl_multifit_nlinear.h>
#include <functional>
#include <tuple>
#include <cassert>
#include "array.h"

template<class R, class... ARGS>
struct function_ripper {
	static constexpr size_t n_args = sizeof...(ARGS);
};

template <class R, class... ARGS>
auto constexpr n_params(R (ARGS...) ){
    return function_ripper<R, ARGS...>();
}

template <typename F, size_t... Is>
auto gen_tuple_impl(F func, std::index_sequence<Is...> ) {
	return std::make_tuple(func(Is)...);
}

template <size_t N, typename F>
auto gen_tuple(F func) {
	return gen_tuple_impl(func, std::make_index_sequence<N>{} );
}

template<typename C>
struct fit_data{
	const Array<double>& t;
	const Array<double>& y;
	C f;
};

template<typename FitData, int n_params>
int internal_f(const gsl_vector* x, void* params, gsl_vector *f) {
	auto* d  = static_cast<FitData*>(params);
	// Convert the parameter values from gsl_vector (in x) into std::tuple
	auto init_args = [x](int index){ return gsl_vector_get(x, index); };
	auto parameters = gen_tuple<n_params>(init_args);

	// Calculate the error for each...
	for (size_t i = 0; i < d->t.get_size(); ++i){
		double ti = d->t[i];
		double yi = d->y[i];
		auto func = [ti, &d](auto ...xs) { return d->f(ti, xs...); };
		auto y = apply(func, parameters);
		gsl_vector_set(f, i, yi - y);
	}
	return GSL_SUCCESS;
}

using func_f_type   = int (*) (const gsl_vector*, void*, gsl_vector*);
using func_df_type  = int (*) (const gsl_vector*, void*, gsl_matrix*);
using func_fvv_type = int (*) (const gsl_vector*, const gsl_vector *, void *, gsl_vector *);

gsl_vector* internal_make_gsl_vector_ptr(const Array<double>& vec);
Array<double> internal_solve_system(gsl_vector* initial_params, gsl_multifit_nlinear_fdf *fdf, gsl_multifit_nlinear_parameters *params);

template<typename C>
Array<double> curve_fit_impl(func_f_type f, func_df_type df, func_fvv_type fvv, gsl_vector* initial_params, fit_data<C>& fd) {

	assert(fd.t.get_size() == fd.y.get_size());

	auto fdf = gsl_multifit_nlinear_fdf();
	auto fdf_params = gsl_multifit_nlinear_default_parameters();

	fdf.f   = f;
	fdf.df  = df;
	fdf.fvv = fvv;
	fdf.n   = fd.t.get_size();
	fdf.p   = initial_params->size;
	fdf.params = &fd;

	// "This selects the Levenberg-Marquardt algorithm with geodesic acceleration."
	fdf_params.trs = gsl_multifit_nlinear_trs_lmaccel;
	return internal_solve_system(initial_params, &fdf, &fdf_params);
}

/**
 * Performs a non-linear least-squares fit.
 * 
 * @param f a function of type double (double x, double c1, double c2, ..., double cn)
 * where  c1, ..., cn are the coefficients to be fitted.
 * @param initial_params intial guess for the parameters. The size of the array must to 
 * be equal to the number of coefficients to be fitted.
 * @param x the idependent data.
 * @param y the dependent data, must to have the same size as x.
 * @return Array<double> with the computed coefficients
 */
template<typename Callable>
Array<double> curve_fit(Callable f, const Array<double>& initial_params, const Array<double>& x, const Array<double>& y) {
	// We can't pass lambdas without convert to std::function.
	constexpr auto n = decltype(n_params(f))::n_args - 1;
	assert(initial_params.get_size() == n);

	auto params = internal_make_gsl_vector_ptr(initial_params);
	auto fd = fit_data<Callable>{x, y, f};
	return curve_fit_impl(internal_f<decltype(fd), n>, nullptr, nullptr, params,  fd);
}
#endif //__ATAT_CURVE_FIT_HPP__
