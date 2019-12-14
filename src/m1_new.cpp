#include <cmath>
#include "basis/lsjt_scheme.h"
#include "rme_extras.h"
#include "constants.h"
#include "utility.h"
#include "integrals.h"
#include "threedho.h"
#include "m1.h"

namespace chiral
{
  // Leading order matrix element.
  // Under the LENPIC power counting, there is no contribution to M1 at LO.
  double M1Operator::LOMatrixElement(const basis::RelativeStateLSJT& bra,
                                     const basis::RelativeStateLSJT& ket,
                                     const ho::OscillatorParameter& b,
                                     const bool& regularize,
                                     const double& regulator,
                                     const std::size_t& T0,
                                     const std::size_t& Abody)
  {
    return 0;
  }

  double M1Operator::LOMatrixElement(const basis::RelativeCMStateLSJT& bra,
                                     const basis::RelativeCMStateLSJT& ket,
                                     const ho::OscillatorParameter& b,
                                     const bool& regularize,
                                     const double& regulator,
                                     const std::size_t& T0,
                                     const std::size_t& Abody)
  {
    return 0;
  }

  double M1Operator::NLOMatrixElement(const basis::RelativeStateLSJT& bra,
                                      const basis::RelativeStateLSJT& ket,
                                      const ho::OscillatorParameter& b,
                                      const bool& regularize,
                                      const double& regulator,
                                      const std::size_t& T0,
                                      const std::size_t& Abody)
  {
    if (Abody == 1)
      return NLO1Body(bra, ket, T0);
    else if (Abody == 2)
      return NLO2Body(bra, ket, b, regularize, regulator, T0);
    else
      return 0;
  }

  double M1Operator::NLOMatrixElement(const basis::RelativeCMStateLSJT& bra,
                                      const basis::RelativeCMStateLSJT& ket,
                                      const ho::OscillatorParameter& b,
                                      const bool& regularize,
                                      const double& regulator,
                                      const std::size_t& T0,
                                      const std::size_t& Abody)
  {
    if (Abody == 1)
      return NLO1Body(bra, ket, T0);
    else if (Abody == 2)
      return NLO2Body(bra, ket, b, regularize, regulator, T0);
    else
      return 0;
  }

  // Next to next to leading order. There are no chiral eft correction at N2LO.
  double M1Operator::N2LOMatrixElement(const basis::RelativeStateLSJT& bra,
                                       const basis::RelativeStateLSJT& ket,
                                       const ho::OscillatorParameter& b,
                                       const bool& regularize,
                                       const double& regulator,
                                       const std::size_t& T0,
                                       const std::size_t& Abody)
  {
    return 0;
  }

  double M1Operator::N2LOMatrixElement(const basis::RelativeCMStateLSJT& bra,
                                       const basis::RelativeCMStateLSJT& ket,
                                       const ho::OscillatorParameter& b,
                                       const bool& regularize,
                                       const double& regulator,
                                       const std::size_t& T0,
                                       const std::size_t& Abody)
  {
    return 0;
  }

    double M1Operator::N3LOMatrixElement(const basis::RelativeStateLSJT& bra,
                                         const basis::RelativeStateLSJT& ket,
                                         const ho::OscillatorParameter& b,
                                         const bool& regularize,
                                         const double& regulator,
                                         const std::size_t& T0,
                                         const std::size_t& Abody)
  {
    auto result = N3LO2BodyIsoscalar(bra, ket, b, regularize, regulator, T0);
    return result;
  }

  double M1Operator::N3LOMatrixElement(const basis::RelativeCMStateLSJT& bra,
                                       const basis::RelativeCMStateLSJT& ket,
                                       const ho::OscillatorParameter& b,
                                       const bool& regularize,
                                       const double& regulator,
                                       const std::size_t& T0,
                                       const std::size_t& Abody)
  {
    if (Abody == 1)
      return 0;
    else if (Abody == 2)
      return N3LO2BodyIsoscalar(bra, ket, b, regularize, regulator, T0);
    else
      return 0;
  }

  // Next to next to next to next to leading order.
  // At present there are no results for N4LO.
  double M1Operator::N4LOMatrixElement(const basis::RelativeStateLSJT& bra,
                                       const basis::RelativeStateLSJT& ket,
                                       const ho::OscillatorParameter& b,
                                       const bool& regularize,
                                       const double& regulator,
                                       const std::size_t& T0,
                                       const std::size_t& Abody)
  {
    return 0;
  }

  double M1Operator::N4LOMatrixElement(const basis::RelativeCMStateLSJT& bra,
                                       const basis::RelativeCMStateLSJT& ket,
                                       const ho::OscillatorParameter& b,
                                       const bool& regularize,
                                       const double& regulator,
                                       const std::size_t& T0,
                                       const std::size_t& Abody)
  {
    return 0;
  }

  // Next to leading order reduced matrix element.
  // There are both one body and two body corrections. One body correction is
  // the same as impulse approximation. Two body correction is isovector in
  // nature, and has both center of mass and relative components. The one body
  // component is not regularized.

  // Relative NLO reduced matrix element.

  double NLO1Body(const basis::RelativeStateLSJT& bra,
                  const basis::RelativeStateLSJT& ket,
                  const std::size_t& T0)
  {
    int nr = ket.n(), nrp = bra.n();
    int L = ket.L(), Lp = bra.L();
    int S = ket.S(), Sp = bra.S();
    int J = ket.J(), Jp = bra.J();
    int T = ket.T(), Tp = bra.T();

    if (am::AllowedTriangle(Jp, J, 1) && nrp == nr && Lp == L)
      {
        // Isoscalar term, \frac{1}{2}\vec{L}_r + μS \vec{S}.
        if (T0 == 0 && Sp == S && Tp == T)
          {
            double oam_term = 0.5 * Hat(L) * ParitySign(L + J + 1 + S);
            oam_term *= (am::Wigner6J(L, L, 1, Jp, J, S) * std::sqrt(L * (L + 1)));

            double spin_term = (constants::isoscalar_nucleon_magnetic_moment
                                    * Hat(S) * ParitySign(S + Jp + L + 1));
            spin_term *= (am::Wigner6J(S, S, 1, J, Jp, L) * std::sqrt(S * (S + 1)));

            double result = Hat(J) * (oam_term + spin_term);
            return result;
          }
        // Isovector term, symmetric. \frac{1}{2}Tz\vec{L}_r + μV Tz\vec{S}.
        if (T0 == 1 && Sp == S && Tp == T)
          {
            double oam_term = 0.5 * Hat(L) * ParitySign(L + S + J + 1);
            oam_term *= (am::Wigner6J(L, L, 1, Jp, J, S)
                         * std::sqrt(L * (L + 1)));

            double spin_term = (constants::isovector_nucleon_magnetic_moment
                                * Hat(S) * ParitySign(L + S + Jp + 1));
            spin_term *= (am::Wigner6J(S, S, 1, J, Jp, L)
                          * std::sqrt(S * (S + 1)));

            double result = Hat(J) * std::sqrt(T * (T + 1)) * (oam_term + spin_term);
            return result;
          }
        // Isovector term, antisymmetric. μV \frac{τ1z-τ2z}{2}\frac{\vec{σ}1z-\vec{σ}2z}{2}.
        if (T0 == 1 && Sp != S && Tp != T)
          {
            double spin_term = (constants::isovector_nucleon_magnetic_moment
                                * Hat(Sp) * ParitySign(L + S + Jp + 1));
            spin_term *= (am::Wigner6J(Sp, S, 1, J, Jp, L)
                          * Hat(S) * (Sp - S));

            double isospin_term = Hat(T) * (Tp - T);

            double result = Hat(J) * spin_term * isospin_term;
            return result;
          }
        return 0;
      }
    return 0;
  }

  double NLO2Body(const basis::RelativeStateLSJT& bra,
                  const basis::RelativeStateLSJT& ket,
                  const ho::OscillatorParameter& b,
                  const bool& regularize,
                  const double& regulator,
                  const std::size_t& T0)
  {
    if (T0 != 1)
      return 0;

    int nr = ket.n(), nrp = bra.n();
    int L = ket.L(), Lp = bra.L();
    int S = ket.S(), Sp = bra.S();
    int J = ket.J(), Jp = bra.J();
    int T = ket.T(), Tp = bra.T();

    double mpi = constants::pion_mass_fm;
    double mN = constants::nucleon_mass_fm;
    double gpi = constants::gA / (2 * constants::pion_decay_constant_fm);
    double brel = b.relative();

    if (am::AllowedTriangle(J, 1, Jp) && Sp != S && Tp != T)
      {
        if (Sp < S && Lp == Jp)
          {
            double result = ParitySign(L + Jp) * std::sqrt(30) * Hat(J);
            result *= am::Wigner6J(2, 1, 1, J, L, Jp) * am::SphericalHarmonicCRME(Lp, L, 2);
            result *= quadrature::IntegralZPiYPiR(nrp, Lp, nr, L, brel, mpi, regularize, regulator);

            if (Lp == L)
              {
                double s1_term = ParitySign(J + Jp + 1) * Hat(J) / Hat(Jp);
                s1_term *= quadrature::IntegralTPiYPiR(nrp, Lp, nr, L, brel, mpi, regularize, regulator);
                result += s1_term;
              }

            result *= Hat(T);
            result *= -mpi * mN * square(gpi) / (3 * constants::pi);
            return result;
          }
        if (Sp > S && L == J)
          {
            double result = ParitySign(L + Jp + 1) * std::sqrt(30) * Hat(Lp);
            result *= am::Wigner6J(1, 1, 2, J, Lp, Jp) * am::SphericalHarmonicCRME(Lp, L, 2);
            result *= quadrature::IntegralZPiYPiR(nrp, Lp, nr, L, brel, mpi, regularize, regulator);

            if (Lp == L)
              {
                result += quadrature::IntegralTPiYPiR(nrp, Lp, nr, L, brel, mpi, regularize, regulator);
              }

            result *= Hat(T);
            result *= -mpi * mN * square(gpi) / (3 * constants::pi);
            return result;
          }
        return 0;
      }
    return 0;
  }

  // \hat{l'}⟨n'l'||∇||nl⟩.
  double GradRME(const int& np,
                 const int& lp,
                 const int& n,
                 const int& l,
                 const double& b)
  {
    if (lp == l + 1)
      {
        double result = 0;
        if (np + 1 == n)
          result += std::sqrt(np + 1);
        if (np == n)
          result += std::sqrt(n + l + 1.5);
        result *= (-std::sqrt(l + 1) / b);
        return result;
      }
    if (lp + 1 == l)
      {
        return GradRME(n, l, np, lp, b);
      }
    return 0;
  }

  // \hat{l'}⟨n'l'||r||nl⟩.
  double VecRRME(const int& np,
                 const int& lp,
                 const int& n,
                 const int& l,
                 const double& b)
  {
    if (lp == l + 1)
      {
        double result = 0;
        if (np == n)
          result += std::sqrt(n + l + 1.5);
        if (np + 1 == n)
          result -= std::sqrt(np + 1);
        result *= std::sqrt(l + 1) * b;
        return result;
      }
    if (lp + 1 == l)
      {
        return -VecRRME(n, l, np, lp, b);
      }
    return 0;
  }

  double NLO1Body(const basis::RelativeCMStateLSJT& bra,
                  const basis::RelativeCMStateLSJT& ket,
                  const std::size_t& T0)
  {
    int nr = ket.Nr(), nrp = bra.Nr();
    int lr = ket.lr(), lrp = bra.lr();
    int nc = ket.Nc(), ncp = bra.Nc();
    int lc = ket.lc(), lcp = bra.lc();
    int L = ket.L(), Lp = bra.L();
    int S = ket.S(), Sp = bra.S();
    int J = ket.J(), Jp = bra.J();
    int T = ket.T(), Tp = bra.T();

    if (am::AllowedTriangle(J, Jp, 1))
      {
        if (T0 == 0)
          {
            if (nrp == nr && lrp == lr && ncp == nc && lcp == lc && Sp == S && Tp == T)
              {
                double result = (am::Wigner9J(lr, lc, L, 1, 0, 1, lrp, lcp, Lp) * std::sqrt(lr * (lr + 1))
                                 + am::Wigner9J(lr, lc, L, 0, 1, 1, lrp, lcp, Lp) * std::sqrt(lc * (lc + 1)));
                result *= (Hat(1) * am::Wigner9J(L, S, J, 1, 0, 1, Lp, Sp, Jp) / 2.0);

                if (Lp == L)
                  {
                    result += (constants::isoscalar_nucleon_magnetic_moment
                               * am::Wigner9J(L, S, J, 0, 1, 1, Lp, Sp, Jp)
                               * std::sqrt(S * (S + 1)) / HatProduct(lr, lc, L));
                  }

                result *= HatProduct(Lp, Sp, J, 1, lrp, lcp, L);

                return result;
              }
            return 0;
          }

        if (T0 == 1)
          {
            if (Tp == T)
              {
                if(nrp == nr && lrp == lr && ncp == nc && lcp == lc && Sp == S)
                  {
                    double result = (am::Wigner9J(lr, lc, L, 1, 0, 1, lrp, lcp, Lp) * std::sqrt(lr * (lr + 1))
                                     + am::Wigner9J(lr, lc, L, 0, 1, 1, lrp, lcp, Lp) * std::sqrt(lc * (lc + 1)));
                    result *= (Hat(1) * am::Wigner9J(L, S, J, 1, 0, 1, Lp, Sp, Jp) / 2.0);

                    if (Lp == L)
                      {
                        result += (constants::isovector_nucleon_magnetic_moment
                                   * am::Wigner9J(L, S, J, 0, 1, 1, Lp, Sp, Jp)
                                   * std::sqrt(S * (S + 1)) / HatProduct(lr, lc, L));
                      }

                    result *= HatProduct(Lp, Sp, J, 1, lrp, lcp, L);

                    result *= std::sqrt(T * (T + 1));

                    return result;
                  }
                return 0;
              }

            if (Tp != T)
              {
                double result = Hat(T) * (Tp - T);

                if (Sp == S)
                  {
                    double hat_product = 3 * std::sqrt(2) * HatProduct(Lp, Sp, J, L);
                    double wigner_product = (am::Wigner9J(L, S, J, 1, 0, 1, Lp, Sp, Jp)
                                             * am::Wigner9J(lr, lc, L, 1, 1, 1, lrp, lcp, Lp));
                    double rme = (GradRME(nrp, lrp, nr, lr, 1) * VecRRME(ncp, lcp, nc, lc, 1)
                                  - VecRRME(nrp, lrp, nr, lr, 1) * GradRME(ncp, lcp, nc, lc, 1));
                    result *= (0.5 * hat_product * wigner_product * rme);
                  }

                if (Sp != S)
                  {
                    if (nrp != nr || lrp != lr || ncp != nc || lcp != lc || Lp != L)
                      return 0;

                    double hat_product = HatProduct(Lp, Sp, J, 1, S);
                    double wigner_product = am::Wigner9J(L, S, J, 0, 1, 1, Lp, Sp, Jp);
                    double rme = (constants::isovector_nucleon_magnetic_moment * (Sp - S));
                    result *= (hat_product * wigner_product * rme);
                  }
                return result;
              }
          }
      }
    return 0;
  }


  double U1aRME(const int& lrp,
                const int& lcp,
                const int& Lp,
                const int& Sp,
                const int& Jp,
                const int& lr,
                const int& lc,
                const int& L,
                const int& S,
                const int& J)
  {
    if (Sp == S)
      {
        double result = 3 * ParitySign(Lp + J) * HatProduct(Lp, J, lrp, lcp, L);
        result *= (am::Wigner6J(L, Lp, 1, Jp, J, S)
                   * am::Wigner9J(lr, lc, L, 1, 1, 1, lrp, lcp, Lp));
        result *= (am::SphericalHarmonicCRME(lrp, lr, 1)
                   * am::SphericalHarmonicCRME(lcp, lc, 1));
        result *= (Hat(1 - S) / Hat(S));
        return result;
      }
    return 0;
  }

  double U1bcdeRME(const int& lrp,
                   const int& lcp,
                   const int& Lp,
                   const int& Sp,
                   const int& Jp,
                   const int& lr,
                   const int& lc,
                   const int& L,
                   const int& S,
                   const int& J)
  {
    if (Sp == 1 && S == 1)
      {
        double result_b = (3 * am::Wigner9J(L, 1, J, 1, 2, 1, Lp, 1, Jp)
                           * am::Wigner9J(lr, lc, L, 1, 1, 1, lrp, lcp, Lp)
                           * am::SphericalHarmonicCRME(lrp, lr, 1));

        double result_c = (3 * std::sqrt(5)
                           * am::Wigner9J(L, 1, J, 2, 2, 1, Lp, 1, Jp)
                           * am::Wigner9J(lr, lc, L, 1, 1, 1, lrp, lcp, Lp)
                           * am::SphericalHarmonicCRME(lrp, lr, 1));

        double result_d = (std::sqrt(70)
                           * am::Wigner9J(L, 1, J, 2, 2, 1, Lp, 1, Jp)
                           * am::Wigner9J(lr, lc, L, 3, 1, 2, lrp, lcp, Lp)
                           * am::SphericalHarmonicCRME(lrp, lr, 3));

        double result_e = (7 * am::Wigner9J(L, 1, J, 3, 2, 1, Lp, 1, Jp)
                           * am::Wigner9J(lr, lc, L, 3, 1, 3, lrp, lcp, Lp)
                           * am::SphericalHarmonicCRME(lrp, lr, 3));

        double result = result_b + result_c + result_d + result_e;
        result *= 2 * std::sqrt(3) * HatProduct(Lp, J, lrp, lcp, L);
        result *= am::SphericalHarmonicCRME(lcp, lc, 1);
        return result;
      }
    return 0;
  }

  double U1fS1RME(const int& lrp,
                  const int& lcp,
                  const int& Lp,
                  const int& Sp,
                  const int& Jp,
                  const int& lr,
                  const int& lc,
                  const int& L,
                  const int& S,
                  const int& J)
  {
    double result = std::abs(Sp - S);
    result *= HatProduct(Lp, Sp, J, 1, lrp, lcp, L, 2, S);
    result *= (am::Wigner9J(L, S, J, 2, 1, 1, Lp, Sp, Jp)
               * am::Wigner9J(lr, lc, L, 2, 0, 2, lrp, lcp, Lp));
    result *= am::SphericalHarmonicCRME(lrp, lr, 2);
    result *= 2 * std::sqrt(5);
    return result;
  }

  double S1RME(const int& lrp,
               const int& lcp,
               const int& Lp,
               const int& Sp,
               const int& Jp,
               const int& lr,
               const int& lc,
               const int& L,
               const int& S,
               const int& J)
  {
    if (Sp != S && lrp == lr && lcp == lc && Lp == L)
      {
        if (Sp < S)
          return std::sqrt(2) * ParitySign(J + Jp + 1) * Hat(J) / Hat(Jp);

        if (Sp > S)
          return std::sqrt(2);
      }
    return 0;
  }

  double NLO2Body(const basis::RelativeCMStateLSJT& bra,
                  const basis::RelativeCMStateLSJT& ket,
                  const ho::OscillatorParameter& b,
                  const bool& regularize,
                  const double& regulator,
                  const std::size_t& T0)
  {
    if (T0 != 1)
      return 0;

    int nr = ket.Nr(), nrp = bra.Nr();
    int lr = ket.lr(), lrp = bra.lr();
    int nc = ket.Nc(), ncp = bra.Nc();
    int lc = ket.lc(), lcp = bra.lc();
    int L = ket.L(), Lp = bra.L();
    int S = ket.S(), Sp = bra.S();
    int J = ket.J(), Jp = bra.J();
    int T = ket.T(), Tp = bra.T();

    if (am::AllowedTriangle(Jp, 1, J) && Tp != T)
      {
        double mpi = constants::pion_mass_fm;
        double mN = constants::nucleon_mass_fm;
        double gpi = constants::gA / (2 * constants::pion_decay_constant_fm);
        double brel = b.relative();
        double bcm = b.cm();

        if (Sp == S)
          {
            double u1a_rme = U1aRME(lrp, lcp, Lp, Sp, Jp, lr, lc, S, S, J);
            double mpi_ypi_integral
              = quadrature::IntegralMPiRYPiR(nrp, lrp, nr, lr, brel, mpi, regularize, regulator);

            double result = mpi_ypi_integral * u1a_rme;

            if (Sp == 1)
              {
                double u1bcde_rme = U1bcdeRME(lrp, lcp, Lp, Sp, Jp, lr, lc, L, S, J);
                double mpi_wpi_ypi_integral
                  = quadrature::IntegralMPiRWPiRYPiR(nrp, lrp, nr, lr, brel, mpi, regularize, regulator);

                result += mpi_wpi_ypi_integral * u1bcde_rme;
              }

            double cm_integral = quadrature::IntegralMPiR(ncp, lcp, nc, lc, bcm, mpi);
            result *= cm_integral;

            result *= std::sqrt(2) * Hat(T);
            result *= -mpi * mN * square(gpi) / (6 * constants::pi);

            return result;
          }

        if (Sp != S && ncp == nc && lcp == lc)
          {
            double u1f_s1_rme = U1fS1RME(lrp, lcp, Lp, Sp, Jp, lr, lc, L, S, J);
            double zpi_ypi_integral
              = quadrature::IntegralZPiYPiR(nrp, lrp, nr, lr, brel, mpi, regularize, regulator);

            double result = u1f_s1_rme * zpi_ypi_integral;

            if (lrp == lr && Lp == L)
              {
                double s1_rme = S1RME(lrp, lcp, Lp, Sp, Jp, lr, lc, L, S, J);
                double tpi_ypi_integral
                  = quadrature::IntegralTPiYPiR(nrp, lrp, nr, lr, brel, mpi, regularize, regulator);

                result += tpi_ypi_integral * s1_rme;
              }

            result *= std::sqrt(2) * Hat(T);
            result *= -mpi * mN * square(gpi) / (6 * constants::pi);

            return result;
          }
        return 0;
      }
    return 0;
  }

  // Next to next to next to leading order.
  // There are both isoscalar and isovector two body chiral eft corrections at
  // N3LO. Currently only the isoscalar part has been implemented, as that is
  // important for the deuteron.

  double N3LO2BodyIsoscalar(const basis::RelativeStateLSJT& bra,
                            const basis::RelativeStateLSJT& ket,
                            const ho::OscillatorParameter& b,
                            const bool& regularize,
                            const double& regulator,
                            const std::size_t& T0)
  {
    // // Relative quantum numbers.
    // int nr = ket.n(), nrp = bra.n();
    // int L = ket.L(), Lp = bra.L();
    // int S = ket.S(), Sp = bra.S();
    // int J = ket.J(), Jp = bra.J();
    // int T = ket.T(), Tp = bra.T();

    // // Spin rmes.
    // auto S_rme = am::RelativeSpinSymmetricRME(Lp, L, Sp, S, Jp, J, 0, 1);

    // // GSL parameters for radial integrals.
    // auto brel = b.relative();
    // auto scaled_regulator_rel = regulator / brel;
    // auto scaled_pion_mass_rel = constants::pion_mass_fm * brel;
    // quadrature::gsl_params_2n prel{nrp, Lp, nr, L, regularize, scaled_pion_mass_rel, scaled_regulator_rel};

    // // d9 term.
    // // d9 isospin rme.
    // auto T0_rme = am::PauliProductRME(Tp, T, 0);
    // // d9 radial integrals.
    // auto norm_product = (ho::CoordinateSpaceNorm(nr, L, 1)
    //                      * ho::CoordinateSpaceNorm(nrp, Lp, 1));
    // auto ypi_integral = norm_product * quadrature::IntegralYPiR(prel);
    // auto wpi_integral = norm_product * quadrature::IntegralWPiRYPiR(prel);
    // // d9  angular momentum rmes.
    // auto A6S_rme = (std::sqrt(10) * am::RelativeSpinSymmetricRME(Lp, L, Sp, S, Jp, J, 2, 1));
    // // d9 prefactor.
    // auto d9_prefactor = (constants::gA * constants::d9_fm
    //                      * cube(constants::pion_mass_fm));
    // d9_prefactor /= (std::sqrt(3) * constants::pi
    //                  * square(constants::pion_decay_constant_fm));
    // auto d9_term = (d9_prefactor * T0_rme
    //                 * (wpi_integral * A6S_rme - ypi_integral * S_rme));

    // // L2 term.
    // double L2_term = 0;
    // if (L == 0 && L == 0 && Tp == T)
    //   {
    //     auto delta_integral = quadrature::IntegralRegularizedDelta(prel) / cube(brel);
    //     L2_term += (2 * constants::L2_fm * S_rme * delta_integral);
    //   }

    // // Overall result.
    // auto result = d9_term + L2_term;
    // result *= 2 * constants::nucleon_mass_fm;
    // if (isnan(result))
    //   result = 0;
    // return result;
    return 0;
  }

  double N3LO2BodyIsoscalar(const basis::RelativeCMStateLSJT& bra,
                            const basis::RelativeCMStateLSJT& ket,
                            const ho::OscillatorParameter& b,
                            const bool& regularize,
                            const double& regulator,
                            const std::size_t& T0)
  {
    return 0;
  }

}
