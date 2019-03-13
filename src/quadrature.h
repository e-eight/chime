#ifndef QUADRATURE_H
#define QUADRATURE_H

#include <limits>
#include <gsl/gsl_integration.h>

namespace quadrature
{
  template <class F>
  double GaussLaguerre(F func, double a, double b, double alpha, int nodes)
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

  template <class F>
  double GaussLegendre(F func, double a, double b, double alpha, double beta, int nodes)
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
}

#endif
