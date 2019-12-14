#include <cmath>
#include "basis/lsjt_scheme.h"
#include "rme_extras.h"
#include "constants.h"
#include "utility"
#include "integrals.h"
#include "threedho.h"
#include "gt.h"

namespace chiral
{
  // Leading order matrix element.
  double GTOperator::LOMatrixElement(const basis::RelativeStateLSJT& bra,
                                     const basis::RelativeStateLSJT& ket,
                                     const ho::OscillatorParameter& b,
                                     const bool& regularize,
                                     const double& regulator,
                                     const std::size_t& T0,
                                     const std::size_t& Abody)
  {
    if (T0 != 1)
      return 0;
    if (Abody != 1)
      return 0;
    auto result = constants::sqrt2 * LO1Body(bra, ket);
    return result;
  }

  double GTOperator::LOMatrixElement(const basis::RelativeCMStateLSJT& bra,
                                     const basis::RelativeCMStateLSJT& ket,
                                     const ho::OscillatorParameter& b,
                                     const bool& regularize,
                                     const double& regulator,
                                     const std::size_t& T0,
                                     const std::size_t& Abody)
  {
    return 0;
  }


  double GTOperator::NLOMatrixElement(const basis::RelativeStateLSJT& bra,
                                      const basis::RelativeStateLSJT& ket,
                                      const ho::OscillatorParameter& b,
                                      const bool& regularize,
                                      const double& regulator,
                                      const std::size_t& T0,
                                      const std::size_t& Abody)
  {
    return 0;
  }
  double GTOperator::NLOMatrixElement(const basis::RelativeCMStateLSJT& bra,
                                      const basis::RelativeCMStateLSJT& ket,
                                      const ho::OscillatorParameter& b,
                                      const bool& regularize,
                                      const double& regulator,
                                      const std::size_t& T0,
                                      const std::size_t& Abody)
  {
    return 0;
  }

  // Next to next to leading order matrix element
  double GTOperator::N2LOMatrixElement(const basis::RelativeStateLSJT& bra,
                                       const basis::RelativeStateLSJT& ket,
                                       const ho::OscillatorParameter& b,
                                       const bool& regularize,
                                       const double& regulator,
                                       const std::size_t& T0,
                                       const std::size_t& Abody)
  {
    if (T0 != 1)
      return 0;
    if (Abody != 2)
      return 0;
    auto result = (constants::sqrt2 * (c3Term(bra, ket, b, regularize, regulator)
                                       + c4Term(bra, ket, b, regularize, regulator)
                                       + DTerm(bra, ket, b, regularize, regulator)));
    return result;
  }

  double GTOperator::N2LOMatrixElement(const basis::RelativeCMStateLSJT& bra,
                                       const basis::RelativeCMStateLSJT& ket,
                                       const ho::OscillatorParameter& b,
                                       const bool& regularize,
                                       const double& regulator,
                                       const std::size_t& T0,
                                       const std::size_t& Abody)
  {
    return 0;
  }

  double GTOperator::N3LOMatrixElement(const basis::RelativeStateLSJT& bra,
                                       const basis::RelativeStateLSJT& ket,
                                       const ho::OscillatorParameter& b,
                                       const bool& regularize,
                                       const double& regulator,
                                       const std::size_t& T0,
                                       const std::size_t& Abody)
  {
    return 0;
  }
  double GTOperator::N3LOMatrixElement(const basis::RelativeCMStateLSJT& bra,
                                       const basis::RelativeCMStateLSJT& ket,
                                       const ho::OscillatorParameter& b,
                                       const bool& regularize,
                                       const double& regulator,
                                       const std::size_t& T0,
                                       const std::size_t& Abody)
  {
    return 0;
  }

  double GTOperator::N4LOMatrixElement(const basis::RelativeStateLSJT& bra,
                                       const basis::RelativeStateLSJT& ket,
                                       const ho::OscillatorParameter& b,
                                       const bool& regularize,
                                       const double& regulator,
                                       const std::size_t& T0,
                                       const std::size_t& Abody)
  {
    return 0;
  }
  double GTOperator::N4LOMatrixElement(const basis::RelativeCMStateLSJT& bra,
                                       const basis::RelativeCMStateLSJT& ket,
                                       const ho::OscillatorParameter& b,
                                       const bool& regularize,
                                       const double& regulator,
                                       const std::size_t& T0,
                                       const std::size_t& Abody)
  {
    return 0;
  }

  double LO1Body(const basis::RelativeStateLSJT& bra,
                 const basis::RelativeStateLSJT& ket)
  {
    std::size_t nr = ket.n(), nrp = bra.n();
    std::size_t L = ket.L(), Lp = bra.L();
    std::size_t S = ket.S(), Sp = bra.S();
    std::size_t J = ket.J(), Jp = bra.J();
    std::size_t T = ket.T(), Tp = bra.T();

    if (am::AllowedTriangle(Jp, J, 1) && nr == nrp && L == Lp)
      {
        if (Tp == T && Sp == S)
          {
            double result = ParitySign(S + Jp + L + 1) * HatProduct(S, J);
            result *= am::Wigner6J(Sp, S, 1, J, Jp, L);
            result *= (std::sqrt(S * (S + 1) * T * (T + 1)));
            return result;
          }

        if (Tp != T && Sp != S)
          {
            double result = ParitySign(S + Jp + L + 1) * HatProduct(Sp, J, S, T);
            result *= am::Wigner6J(S, S, 1, J, Jp, L);
            result *= ((Sp - S) * (Tp - T));
            return result;
          }
      }
    return 0;
  }

  double c3Term(const basis::RelativeStateLSJT& bra,
                const basis::RelativeStateLSJT& ket,
                const ho::OscillatorParameter& b,
                const bool& regularize,
                const double& regulator)
  {
    std::size_t nr = ket.n(), nrp = bra.n();
    std::size_t L = ket.L(), Lp = bra.L();
    std::size_t S = ket.S(), Sp = bra.S();
    std::size_t J = ket.J(), Jp = bra.J();
    std::size_t T = ket.T(), Tp = bra.T();

    double mpi = constants::pion_mass_fm;
    double gpi = (constants::gA / square(constants::pion_decay_constant_fm));
    double c3 =  constants::c3_fm;
    double brel = b.relative();

    if (am::AllowedTriangle(Jp, J, 1))
      {
        if (Tp == T && Sp == S)
          {
            double result = std::sqrt(10) * HatProduct(Lp, 1);
            result *= am::Wigner9J(L, S, J, 2, 1, 1, Lp, Sp, Jp);
            result *= am::SphericalHarmonicCRME(Lp, L, 2);
            result *= quadrature::IntegralWPiYPiR(nrp, Lp, nr, L, brel, mpi, regularize, regulator);

            if (Lp == L)
              {
                double s1_term = ParitySign(S + L + Jp + 1);
                s1_term *= am::Wigner6J(Sp, S, 1, J, Jp, L);
                s1_term *= quadrature::IntegralYPiR(nrp, Lp, nr, L, brel, mpi, regularize, regulator);
                result += s1_term;
              }

            result *= HatProduct(Sp, J);
            result *= std::sqrt(S * (S + 1) * T * (T + 1));

            result *= (gpi * c3 * cube(mpi) / 6);
            return result;
          }

        if (Tp != T && Sp != S)
          {
            double result = std::sqrt(10) * HatProduct(Lp, 1);
            result *= am::Wigner9J(L, S, J, 2, 1, 1, Lp, Sp, Jp);
            result *= am::SphericalHarmonicCRME(Lp, L, 2);
            result *= quadrature::IntegralWPiYPiR(nrp, Lp, nr, L, brel, mpi, regularize, regulator);

            if (Lp == L)
              {
                double s1_term = ParitySign(S + L + Jp + 1);
                s1_term *= am::Wigner6J(Sp, S, 1, J, Jp, L);
                s1_term *= quadrature::IntegralYPiR(nrp, Lp, nr, L, brel, mpi, regularize, regulator);
                result += s1_term;
              }

            result *= HatProduct(Sp, J, S, T);
            result *= (Sp - S) * (Tp - T);

            result *= (gpi * c3 * cube(mpi) / 6);
            return result;
          }
        return 0;
      }
    return 0;
  }

  double c4Term(const basis::RelativeStateLSJT& bra,
                const basis::RelativeStateLSJT& ket,
                const ho::OscillatorParameter& b,
                const bool& regularize,
                const double& regulator)
  {
    std::size_t nr = ket.n(), nrp = bra.n();
    std::size_t L = ket.L(), Lp = bra.L();
    std::size_t S = ket.S(), Sp = bra.S();
    std::size_t J = ket.J(), Jp = bra.J();
    std::size_t T = ket.T(), Tp = bra.T();

    double mpi = constants::pion_mass_fm;
    double gpi = (constants::gA / square(constants::pion_decay_constant_fm));
    double c4 =  constants::c4_fm;
    double brel = b.relative();

    if (am::AllowedTriangle(Jp, J, 1) && Tp != T && Sp != S)
      {
        double result = std::sqrt(10) * HatProduct(Lp, 1);
        result *= am::Wigner9J(L, S, J, 2, 1, 1, Lp, Sp, Jp);
        result *= am::SphericalHarmonicCRME(Lp, L, 2);
        result *= quadrature::IntegralWPiYPiR(nrp, Lp, nr, L, brel, mpi, regularize, regulator);

        if (Lp == L)
          {
            double s1_term = 2 * ParitySign(S + L + Jp + 1);
            s1_term *= am::Wigner6J(Sp, S, 1, J, Jp, L);
            s1_term *= quadrature::IntegralYPiR(nrp, Lp, nr, L, brel, mpi, regularize, regulator);
            result += s1_term;
          }

        result *= HatProduct(Sp, J, S, T);

        result *= (gpi * c4 * cube(mpi) / 6);
        return result;
      }
    return 0;
  }

  double DTerm(const basis::RelativeStateLSJT& bra,
               const basis::RelativeStateLSJT& ket,
               const ho::OscillatorParameter& b,
               const bool& regularize,
               const double& regulator)
  {
    std::size_t nr = ket.n(), nrp = bra.n();
    std::size_t L = ket.L(), Lp = bra.L();
    std::size_t S = ket.S(), Sp = bra.S();
    std::size_t J = ket.J(), Jp = bra.J();
    std::size_t T = ket.T(), Tp = bra.T();

    double brel = b.relative();
    double D = constants::D_fm;

    if (am::AllowedTriangle(J, Jp, 1) && L == 0 && Lp == 0)
      {
        if (Tp == T && Sp == S)
          {
            double result = (quadrature::IntegralRegularizedDelta(nr, brel, regulator)
                             * quadrature::IntegralRegularizedDelta(nrp, brel, regulator));
            result *= std::sqrt(S * (S + 1) * T * (T + 1));
            result *= (-0.5 * D);
            return result;
          }
        if (Tp != T && Sp != S)
          {
            double result = (quadrature::IntegralRegularizedDelta(nr, brel, regulator)
                             * quadrature::IntegralRegularizedDelta(nrp, brel, regulator));
            result *= HatProduct(S, T) * (Sp - S) * (Tp - T);
            result *= (-0.5 * D);
            return result;
          }
        return 0;
      }
    return 0;
  }

}
