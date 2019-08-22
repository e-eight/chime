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
                                     const double& regulator)
  {
    std::size_t nr = ket.n(), nrp = bra.n();
    std::size_t lr = ket.L(), lrp = bra.L();
    std::size_t s = ket.S(), sp = bra.S();
    std::size_t j = ket.J(), jp = bra.J();
    std::size_t t = ket.T(), tp = bra.T();

    bool kronecker = (nr == nrp && lr == lrp);
    if (!kronecker)
      return 0;

    auto symm_term = am::RelativeSpinSymmetricRME(lrp, lr, sp, s, jp, j, 0, 1);
    symm_term *= am::SpinSymmetricRME(tp, t);

    auto asymm_term = am::RelativeSpinAntisymmetricRME(lrp, lr, sp, s, jp, j, 0, 1);
    asymm_term *= am::SpinAntisymmetricRME(tp, t);

    auto result = ((-constants::gA / constants::sqrt2pi)
                   * (symm_term + asymm_term));
    if (isnan(result))
      result = 0;
    return result;
  }

  double GTOperator::LOMatrixElement(const basis::RelativeCMStateLSJT& bra,
                                     const basis::RelativeCMStateLSJT& ket,
                                     const util::OscillatorParameter& b,
                                     const bool& regularize,
                                     const double& regulator)
  {
    return 0;
  }


  double GTOperator::NLOMatrixElement(const basis::RelativeStateLSJT& bra,
                                      const basis::RelativeStateLSJT& ket,
                                      const util::OscillatorParameter& b,
                                      const bool& regularize,
                                      const double& regulator)
  {
    return 0;
  }
  double GTOperator::NLOMatrixElement(const basis::RelativeCMStateLSJT& bra,
                                      const basis::RelativeCMStateLSJT& ket,
                                      const util::OscillatorParameter& b,
                                      const bool& regularize,
                                      const double& regulator)
  {
    return 0;
  }

  // Next to next to leading order matrix element
  double GTOperator::N2LOMatrixElement(const basis::RelativeStateLSJT& bra,
                                       const basis::RelativeStateLSJT& ket,
                                       const util::OscillatorParameter& b,
                                       const bool& regularize,
                                       const double& regulator)
  {
    std::size_t nr = ket.n(), nrp = bra.n();
    std::size_t lr = ket.L(), lrp = bra.L();
    std::size_t s = ket.S(), sp = bra.S();
    std::size_t j = ket.J(), jp = bra.J();
    std::size_t t = ket.T(), tp = bra.T();

    //OscillatorParameter and scaling.
    auto brel = b.relative();
    auto scaled_regulator = regulator / brel;
    auto scaled_pion_mass = constants::pion_mass_fm * brel;

    // Parameters for integration routines.
    quadrature::gsl_params_2n p {nrp, lrp, nr, lr, regularize, scaled_regulator, scaled_pion_mass};

    // Common symmetric and asymmetric spin and isospin RMEs.
    auto symm_rme_spin = am::RelativeSpinSymmetricRME(lrp, lr, sp, s, jp, j, 0, 1);
    auto symm_rme_isospin = am::SpinSymmetricRME(tp, t);
    auto asymm_rme_spin = am::RelativeSpinAntisymmetricRME(lrp, lr, sp, s, jp, j, 0, 1);
    auto asymm_rme_isospin = am::SpinAntisymmetricRME(tp, t);

    // Radial integrals.
    auto norm_product = (ho::CoordinateSpaceNorm(nr, lr, 1)
                         * ho::CoordinateSpaceNorm(nrp, lrp, 1));
    auto ypi_integral = norm_product * quadrature::IntegralYPiR(p);
    auto wpi_integral = norm_product * quadrature::IntegralWPiRYPiR(p);
    auto delta_integral = norm_product * quadrature::IntegralRegularizedDelta(p);
    delta_integral /= util::cube(brel);

    // Common part of pion exchange term prefactors.
    auto pion_prefactor = (constants::gA
                           * util::cube(constants::pion_mass_fm));
    pion_prefactor /= (6 * constants::pi
                       * util::square(constants::pion_decay_constant_fm));

    // C3 term.
    // Symmetric term.
    auto symm_term_c3 = symm_rme_spin * ypi_integral;
    symm_term_c3 += (am::RelativeSpinSymmetricRME(lrp, lr, sp, s, jp, j, 2, 1)
                     * wpi_integral);
    symm_term_c3 *= symm_rme_isospin;
    // Asymmetric term.
    auto asymm_term_c3 = asymm_rme_spin * ypi_integral;
    asymm_term_c3 += (am::RelativeSpinAntisymmetricRME(lrp, lr, sp, s, jp, j, 2, 1)
                      * wpi_integral);
    asymm_term_c3 *= asymm_rme_isospin;
    // C3 prefactor.
    auto c3_prefactor = -pion_prefactor * constants::c3_fm;
    // C3 final result.
    auto c3_result = (c3_prefactor * (symm_term_c3 + asymm_term_c3));
    if (isnan(c3_result))
      c3_result = 0;

    // C4 term.
    auto spin_rme = (2 * ypi_integral
                     * am::RelativePauliProductRME(lrp, lr, sp, s, jp, j, 0, 1, 1));
    spin_rme += (std::sqrt(10) * wpi_integral
                 * am::RelativePauliProductRME(lrp, lr, sp, s, jp, j, 2, 2, 1));
    auto isospin_rme = 0.5 * am::PauliProductRME(tp, t, 1);
    auto c4_prefactor = pion_prefactor * constants::c4_fm;
    auto c4_result = (c4_prefactor * isospin_rme * spin_rme);
    if (isnan(c4_result))
      c4_result = 0;

    // D term.
    double D_result = 0;
    if (lrp == 0 && lr == 0 && nrp == nr)
      {
        norm_product = (ho::CoordinateSpaceNorm(nrp, 0, 1)
                         * ho::CoordinateSpaceNorm(nr, 0, 1));
        auto wf_product = sf::Laguerre(nrp, 0.5, 0) * sf::Laguerre(nr, 0.5, 0);
        D_result = ((symm_rme_spin * symm_rme_isospin)
                    + (asymm_rme_spin * asymm_rme_isospin));
        D_result *= (-0.5 * constants::D_fm * norm_product * wf_product);
        D_result /= util::cube(brel);
        if (isnan(D_result))
          D_result = 0;
      }

    // Overall result
    auto result = ((c3_result + c4_result + D_result) / (constants::sqrt2pi));
    return result;
  }

  double GTOperator::N2LOMatrixElement(const basis::RelativeCMStateLSJT& bra,
                                       const basis::RelativeCMStateLSJT& ket,
                                       const util::OscillatorParameter& b,
                                       const bool& regularize,
                                       const double& regulator)
  {
    return 0;
  }

  double GTOperator::N3LOMatrixElement(const basis::RelativeStateLSJT& bra,
                                       const basis::RelativeStateLSJT& ket,
                                       const util::OscillatorParameter& b,
                                       const bool& regularize,
                                       const double& regulator)
  {
    return 0;
  }
  double GTOperator::N3LOMatrixElement(const basis::RelativeCMStateLSJT& bra,
                                       const basis::RelativeCMStateLSJT& ket,
                                       const util::OscillatorParameter& b,
                                       const bool& regularize,
                                       const double& regulator)
  {
    return 0;
  }

  double GTOperator::N4LOMatrixElement(const basis::RelativeStateLSJT& bra,
                                       const basis::RelativeStateLSJT& ket,
                                       const util::OscillatorParameter& b,
                                       const bool& regularize,
                                       const double& regulator)
  {
    return 0;
  }
  double GTOperator::N4LOMatrixElement(const basis::RelativeCMStateLSJT& bra,
                                       const basis::RelativeCMStateLSJT& ket,
                                       const util::OscillatorParameter& b,
                                       const bool& regularize,
                                       const double& regulator)
  {
    return 0;
  }

}
