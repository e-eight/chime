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
                                     const util::OscillatorParameter& b,
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

  double LO1Body(const basis::RelativeStateLSJT& bra,
                 const basis::RelativeStateLSJT& ket)
  {
    std::size_t nr = ket.n(), nrp = bra.n();
    std::size_t L = ket.L(), Lp = bra.L();
    std::size_t S = ket.S(), Sp = bra.S();
    std::size_t J = ket.J(), Jp = bra.J();
    std::size_t T = ket.T(), Tp = bra.T();

    if (nr != nrp || L != Lp)
      return 0;

    // Angular momentum & isospin RMEs.
    auto symm_term = am::RelativeSpinSymmetricRME(Lp, L, Sp, S, Jp, J, 0, 1);
    symm_term *= am::SpinSymmetricRME(Tp, T);
    auto asymm_term = am::RelativeSpinAntisymmetricRME(Lp, L, Sp, S, Jp, J, 0, 1);
    asymm_term *= am::SpinAntisymmetricRME(Tp, T);

    // Result.
    auto result = (-constants::gA * (symm_term + asymm_term));
    return result;
  }

  double GTOperator::LOMatrixElement(const basis::RelativeCMStateLSJT& bra,
                                     const basis::RelativeCMStateLSJT& ket,
                                     const util::OscillatorParameter& b,
                                     const bool& regularize,
                                     const double& regulator,
                                     const std::size_t& T0,
                                     const std::size_t& Abody)
  {
    return 0;
  }


  double GTOperator::NLOMatrixElement(const basis::RelativeStateLSJT& bra,
                                      const basis::RelativeStateLSJT& ket,
                                      const util::OscillatorParameter& b,
                                      const bool& regularize,
                                      const double& regulator,
                                      const std::size_t& T0,
                                      const std::size_t& Abody)
  {
    return 0;
  }
  double GTOperator::NLOMatrixElement(const basis::RelativeCMStateLSJT& bra,
                                      const basis::RelativeCMStateLSJT& ket,
                                      const util::OscillatorParameter& b,
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
                                       const util::OscillatorParameter& b,
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

  double c3Term(const basis::RelativeStateLSJT& bra,
                const basis::RelativeStateLSJT& ket,
                const util::OscillatorParameter& b,
                const bool& regularize,
                const double& regulator)
  {
    std::size_t nr = ket.n(), nrp = bra.n();
    std::size_t L = ket.L(), Lp = bra.L();
    std::size_t S = ket.S(), Sp = bra.S();
    std::size_t J = ket.J(), Jp = bra.J();
    std::size_t T = ket.T(), Tp = bra.T();

    // Isospin RMEs.
    auto symm_rme_isospin = am::SpinSymmetricRME(Tp, T);
    auto asymm_rme_isospin = am::SpinAntisymmetricRME(Tp, T);

    // Angular momentum RMEs.
    auto symm_rme_spin = am::RelativeSpinSymmetricRME(Lp, L, Sp, S, Jp, J, 0, 1);
    auto symm_rme_spin_A6 = (std::sqrt(10)
                             * am::RelativeSpinSymmetricRME(Lp, L, Sp, S, Jp, J, 2, 1));
    auto asymm_rme_spin = am::RelativeSpinAntisymmetricRME(Lp, L, Sp, S, Jp, J, 0, 1);
    auto asymm_rme_spin_A6 = (std::sqrt(10)
                              * am::RelativeSpinAntisymmetricRME(Lp, L, Sp, S, Jp, J, 2, 1));

    // Radial integrals.
    auto brel = b.relative();
    auto scaled_regulator = regulator / brel;
    auto scaled_pion_mass = constants::pion_mass_fm * brel;
    quadrature::gsl_params_2n p {nrp, Lp, nr, L, regularize, scaled_regulator, scaled_pion_mass};
    auto wpi_ypi_integral = quadrature::IntegralWPiRYPiR(p);
    auto ypi_integral = quadrature::IntegralYPiR(p);

    // Result.
    auto prefactor = (2 * constants::gA * constants::c3_fm
                      * cube(constants::pion_mass_fm));
    prefactor /= (12 * square(constants::pion_decay_constant_fm)
                  * constants::pi);
    symm_rme_spin *= symm_rme_isospin;
    symm_rme_spin_A6 *= symm_rme_isospin;
    asymm_rme_spin *= asymm_rme_isospin;
    asymm_rme_spin_A6 *= asymm_rme_isospin;
    auto result = ((symm_rme_spin_A6 + asymm_rme_spin_A6) * wpi_ypi_integral
                   - (symm_rme_spin + asymm_rme_spin) * ypi_integral);
    result *= (constants::sqrt2 * prefactor);
    return result;
  }

  double c4Term(const basis::RelativeStateLSJT& bra,
                const basis::RelativeStateLSJT& ket,
                const util::OscillatorParameter& b,
                const bool& regularize,
                const double& regulator)
  {
    std::size_t nr = ket.n(), nrp = bra.n();
    std::size_t L = ket.L(), Lp = bra.L();
    std::size_t S = ket.S(), Sp = bra.S();
    std::size_t J = ket.J(), Jp = bra.J();
    std::size_t T = ket.T(), Tp = bra.T();

    // Isospin RME.
    auto pp_rme_isospin = am::PauliProductRME(Tp, T, 1);

    // Angular momentum RMEs.
    auto pp_rme_spin = am::RelativePauliProductRME(Lp, L, Sp, S, Jp, J, 0, 1, 1);
    auto pp_rme_spin_A6 = (std::sqrt(10)
                           * am::RelativePauliProductRME(Lp, L, Sp, S, Jp, J, 2, 1, 1));

    // Radial integrals.
    auto brel = b.relative();
    auto scaled_regulator = regulator / brel;
    auto scaled_pion_mass = constants::pion_mass_fm * brel;
    quadrature::gsl_params_2n p {nrp, Lp, nr, L, regularize, scaled_regulator, scaled_pion_mass};
    auto wpi_ypi_integral = quadrature::IntegralWPiRYPiR(p);
    auto ypi_integral = quadrature::IntegralYPiR(p);

    // Result.
    auto prefactor = (constants::gA * constants::c4_fm
                      * cube(constants::pion_mass_fm));
    prefactor /= (12 * square(constants::pion_decay_constant_fm)
                  * constants::pi);
    auto result = (pp_rme_isospin * (pp_rme_spin_A6 * wpi_ypi_integral
                                     + 2 * pp_rme_spin * ypi_integral));
    result *= (constants::sqrt2 * prefactor);
    return result;
  }

  double DTerm(const basis::RelativeStateLSJT& bra,
               const basis::RelativeStateLSJT& ket,
               const util::OscillatorParameter& b,
               const bool& regularize,
               const double& regulator)
  {
    std::size_t nr = ket.n(), nrp = bra.n();
    std::size_t L = ket.L(), Lp = bra.L();
    std::size_t S = ket.S(), Sp = bra.S();
    std::size_t J = ket.J(), Jp = bra.J();
    std::size_t T = ket.T(), Tp = bra.T();

    if (Lp != 0 || L != 0)
      return 0;

    // Angular momentum and isospin RMEs.
    auto symm_term = am::RelativeSpinSymmetricRME(Lp, L, Sp, S, Jp, J, 0, 1);
    symm_term *= am::SpinSymmetricRME(Tp, T);
    auto asymm_term = am::RelativeSpinAntisymmetricRME(Lp, L, Sp, S, Jp, J, 0, 1);
    asymm_term *= am::SpinAntisymmetricRME(Tp, T);

    // Radial integral.
    auto brel = b.relative();
    auto scaled_regulator = regulator / brel;
    auto scaled_pion_mass = constants::pion_mass_fm * brel;
    quadrature::gsl_params_2n p {nrp, Lp, nr, L, regularize, scaled_regulator, scaled_pion_mass};
    auto delta_integral = quadrature::IntegralRegularizedDelta(p) / cube(brel);

    // Result.
    auto result = (symm_term + asymm_term) * delta_integral;
    result *= (-0.5 * constants::D_fm);
    return result;
  }

  double GTOperator::N2LOMatrixElement(const basis::RelativeCMStateLSJT& bra,
                                       const basis::RelativeCMStateLSJT& ket,
                                       const util::OscillatorParameter& b,
                                       const bool& regularize,
                                       const double& regulator,
                                       const std::size_t& T0,
                                       const std::size_t& Abody)
  {
    return 0;
  }

  double GTOperator::N3LOMatrixElement(const basis::RelativeStateLSJT& bra,
                                       const basis::RelativeStateLSJT& ket,
                                       const util::OscillatorParameter& b,
                                       const bool& regularize,
                                       const double& regulator,
                                       const std::size_t& T0,
                                       const std::size_t& Abody)
  {
    return 0;
  }
  double GTOperator::N3LOMatrixElement(const basis::RelativeCMStateLSJT& bra,
                                       const basis::RelativeCMStateLSJT& ket,
                                       const util::OscillatorParameter& b,
                                       const bool& regularize,
                                       const double& regulator,
                                       const std::size_t& T0,
                                       const std::size_t& Abody)
  {
    return 0;
  }

  double GTOperator::N4LOMatrixElement(const basis::RelativeStateLSJT& bra,
                                       const basis::RelativeStateLSJT& ket,
                                       const util::OscillatorParameter& b,
                                       const bool& regularize,
                                       const double& regulator,
                                       const std::size_t& T0,
                                       const std::size_t& Abody)
  {
    return 0;
  }
  double GTOperator::N4LOMatrixElement(const basis::RelativeCMStateLSJT& bra,
                                       const basis::RelativeCMStateLSJT& ket,
                                       const util::OscillatorParameter& b,
                                       const bool& regularize,
                                       const double& regulator,
                                       const std::size_t& T0,
                                       const std::size_t& Abody)
  {
    return 0;
  }

}
