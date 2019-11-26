#ifndef QUADRATURE_H
#define QUADRATURE_H

#include <limits>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>

namespace quadrature
{
  // Wrapper for GSL's adaptive Gaussian quadrature on the semi infinite interval.
  template <class F>
  double QAGIUIntegrate(F func, double a, double epsabs, double epsrel)
  {
    gsl_function f;
    f.function = [](double x, void* p) { return (*static_cast<F*>(p))(x); };
    f.params = &func;

    gsl_integration_workspace* w
      = gsl_integration_workspace_alloc(1000);

    double result, error;
    gsl_integration_qagiu(&f, a, epsabs, epsrel, 1000, w, &result, &error);

    gsl_integration_workspace_free(w);
    return result;
  }

  // Wrapper for GSL's Gauss-Legendre quadrature.
  template <class F>
  double GaussLaguerre(F func, double a, double b, double alpha, std::size_t nodes)
  {
    gsl_function f;
    f.function = [](double x, void* p) { return (*static_cast<F*>(p))(x); };
    f.params = &func;

    const gsl_integration_fixed_type* T = gsl_integration_fixed_laguerre;

    gsl_integration_fixed_workspace* w
      = gsl_integration_fixed_alloc(T, nodes, a, b, alpha, 0);

    double result;
    gsl_integration_fixed(&f, &result, w);

    gsl_integration_fixed_free(w);

    return result;
  }

  // Wrapper for GSL's Gauss-Laguerre quadrature.
  template <class F>
  double GaussLegendre(F func, double a, double b, double alpha, double beta, std::size_t nodes)
  {
    gsl_function f;
    f.function = [](double x, void* p) { return (*static_cast<F*>(p))(x); };
    f.params = &func;

    const gsl_integration_fixed_type* T = gsl_integration_fixed_legendre;

    gsl_integration_fixed_workspace* w
      = gsl_integration_fixed_alloc(T, nodes, a, b, alpha, beta);

    double result;
    gsl_integration_fixed(&f, &result, w);

    gsl_integration_fixed_free(w);

    return result;
  }

  // Wrapper for GSL's cubic spline integration.
  template <class T>
  double CubicSplineIntegrate(const std::vector<T> x, const std::vector<T> y)
  {
    int m = x.size();

    gsl_interp_accel *acc = gsl_interp_accel_alloc();
    gsl_spline *spline = gsl_spline_alloc(gsl_interp_cspline, m);

    gsl_spline_init(spline, x.data(), y.data(), m);

    double integral = gsl_spline_eval_integ(spline, x.front(), x.back(), acc);

    gsl_spline_free (spline);
    gsl_interp_accel_free (acc);

    return integral;
  }
}

#endif
