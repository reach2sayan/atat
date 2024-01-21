#include "normal.h"

/*
 * Lower tail quantile for standard normal distribution function.
 *
 * This function returns an approximation of the inverse cumulative
 * standard normal distribution function.  I.e., given P, it returns
 * an approximation to the X satisfying P = Pr{Z <= X} where Z is a
 * random variable from the standard normal distribution.
 *
 * The algorithm uses a minimax approximation by rational functions
 * and the result has a relative error whose absolute value is less
 * than 1.15e-9.
 *
 * Author:      Peter John Acklam
 * Time-stamp:  2002-06-09 18:45:44 +0200
 * E-mail:      jacklam@math.uio.no
 * WWW URL:     http://www.math.uio.no/~jacklam
 *
 * C implementation adapted from Peter's Perl version
 * downloaded from:
 * http://home.online.no/~pjacklam/notes/invnorm/impl/sprouse/
 */

#include <math.h>
#include <errno.h>

/* Coefficients in rational approximations. */
static const Real invnormalcdf_a[] =
  {
    -3.969683028665376e+01,
    2.209460984245205e+02,
    -2.759285104469687e+02,
    1.383577518672690e+02,
    -3.066479806614716e+01,
     2.506628277459239e+00
  };

static const Real invnormalcdf_b[] =
  {
    -5.447609879822406e+01,
    1.615858368580409e+02,
    -1.556989798598866e+02,
    6.680131188771972e+01,
    -1.328068155288572e+01
  };

static const Real invnormalcdf_c[] =
  {
    -7.784894002430293e-03,
    -3.223964580411365e-01,
    -2.400758277161838e+00,
    -2.549732539343734e+00,
    4.374664141464968e+00,
     2.938163982698783e+00
  };

static const Real invnormalcdf_d[] =
  {
    7.784695709041462e-03,
    3.224671290700398e-01,
    2.445134137142996e+00,
    3.754408661907416e+00
  };

#define invnormalcdf_LOW 0.02425
#define invnormalcdf_HIGH 0.97575

Real invnormalcdf(Real p)
{
  Real q, r;

  errno = 0;

  if (p < 0 || p > 1)
    {
      errno = EDOM;
      return 0.0;
    }
  else if (p == 0)
    {
      errno = ERANGE;
      return -HUGE_VAL /* minus "infinity" */;
    }
  else if (p == 1)
    {
      errno = ERANGE;
      return HUGE_VAL /* "infinity" */;
    }
  else if (p < invnormalcdf_LOW)
    {
      /* Rational approximation for lower region */
      q = sqrt(-2*log(p));
      return (((((invnormalcdf_c[0]*q+invnormalcdf_c[1])*q+invnormalcdf_c[2])*q+invnormalcdf_c[3])*q+invnormalcdf_c[4])*q+invnormalcdf_c[5]) /
	((((invnormalcdf_d[0]*q+invnormalcdf_d[1])*q+invnormalcdf_d[2])*q+invnormalcdf_d[3])*q+1);
    }
  else if (p > invnormalcdf_HIGH)
    {
      /* Rational approximation for upper region */
      q  = sqrt(-2*log(1-p));
      return -(((((invnormalcdf_c[0]*q+invnormalcdf_c[1])*q+invnormalcdf_c[2])*q+invnormalcdf_c[3])*q+invnormalcdf_c[4])*q+invnormalcdf_c[5]) /
	((((invnormalcdf_d[0]*q+invnormalcdf_d[1])*q+invnormalcdf_d[2])*q+invnormalcdf_d[3])*q+1);
    }
  else
    {
      /* Rational approximation for central region */
      q = p - 0.5;
      r = q*q;
      return (((((invnormalcdf_a[0]*r+invnormalcdf_a[1])*r+invnormalcdf_a[2])*r+invnormalcdf_a[3])*r+invnormalcdf_a[4])*r+invnormalcdf_a[5])*q /
	(((((invnormalcdf_b[0]*r+invnormalcdf_b[1])*r+invnormalcdf_b[2])*r+invnormalcdf_b[3])*r+invnormalcdf_b[4])*r+1);
    }
}

Real normal01(void) {
  return invnormalcdf(uniform01());
}
