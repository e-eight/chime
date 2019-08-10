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
    symm_term *= 2 * am::SpinSymmetricRME(tp, t);

    auto asymm_term = am::RelativeSpinAsymmetricRME(lrp, lr, sp, s, jp, j, 0, 1);
    asymm_term *= 2 * am::SpinAsymmetricRME(tp, t);

    auto result = (std::sqrt(2) * constants::gA * (symm_term + asymm_term));
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
    quadrature::gsl_params_pion p {nrp, lrp, nr, lr, regularize, scaled_regulator, scaled_pion_mass};
    quadrature::gsl_params_contact c {nr, lr, 2.0 / scaled_regulator};
    quadrature::gsl_params_contact cp {nrp, lrp, 2.0 / scaled_regulator};

    // Symmetric and asymmetric RMEs.
    auto symm_rme_spin = am::SpinSymmetricRME(sp, s);
    auto symm_rme_isospin = 2 * am::SpinSymmetricRME(tp, t);
    auto asymm_rme_spin = am::SpinAsymmetricRME(sp, s);
    auto asymm_rme_isospin = 2 * am::SpinAsymmetricRME(tp, t);

    // Radial integrals.
    auto norm_product = (ho::CoordinateSpaceNorm(nr, lr, 1)
                         * ho::CoordinateSpaceNorm(nrp, lrp, 1));
    auto ypi_integral = norm_product * quadrature::IntegralYPiR(p);
    auto wpi_integral = norm_product * quadrature::IntegralWPiRYPiR(p);
    auto contact_integral = (quadrature::IntegralMomentumGaussian(c)
                             * quadrature::IntegralMomentumGaussian(cp));
    contact_integral *= (norm_product / util::cube(brel));

    // Common part of pion exchange term prefactors.
    auto pion_prefactor = (constants::gA
                           * util::cube(constants::pion_mass_fm));
    pion_prefactor /= (12 * constants::pi
                       * util::square(constants::pion_decay_constant_fm));

    // C3 term.
    // Symmetric term.
    auto symm_term_c3 = symm_rme_spin * ypi_integral;
    symm_term_c3 += (am::RelativeSpinSymmetricRME(lrp, lr, sp, s, jp, j, 2, 1)
                     * wpi_integral);
    symm_term_c3 *= symm_rme_isospin;
    // Asymmetric term.
    auto asymm_term_c3 = asymm_rme_spin * ypi_integral;
    asymm_term_c3 += (am::RelativeSpinAsymmetricRME(lrp, lr, sp, s, jp, j, 2, 1)
                      * wpi_integral);
    asymm_term_c3 *= asymm_rme_isospin;
    // C3 prefactor.
    auto c3_prefactor = pion_prefactor * constants::c3_fm;
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
    auto c4_prefactor = -pion_prefactor * constants::c4_fm;
    auto c4_result = (c4_prefactor * isospin_rme * spin_rme);
    if (isnan(c4_result))
      c4_result = 0;

    // D term.
    // Evaluating the lec D.
    double lec_cD = 0;
    double regulator_tol = 10e-10;
    if (abs(regulator - 0.9) < regulator_tol)
      lec_cD = 1.7;
    if (abs(regulator - 1) < regulator_tol)
      lec_cD = 7.2;
    auto lec_D = (lec_cD / (util::square(constants::pion_decay_constant_fm)
                            * constants::lambda_chiral_fm));
    // D final result.
    double D_result = 0;
    if (lrp == 0 && lr == 0)
      {
        D_result = ((symm_rme_spin * symm_rme_isospin)
                    + (asymm_rme_spin * asymm_rme_isospin));
        D_result *= (0.25 * lec_D * contact_integral);
        if (isnan(D_result))
          D_result = 0;
      }

    // Overall result
    auto result = std::sqrt(2) * (c3_result + c4_result + D_result);
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
