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
                                     const util::OscillatorParameter& b,
                                     const bool& regularize,
                                     const double& regulator)
  {
    return 0;
  }

  double M1Operator::LOMatrixElement(const basis::RelativeCMStateLSJT& bra,
                                     const basis::RelativeCMStateLSJT& ket,
                                     const util::OscillatorParameter& b,
                                     const bool& regularize,
                                     const double& regulator)
  {
    return 0;
  }

  // Next to leading order matrix element.
  // There are both one body and two body corrections. One body correction is
  // the same as impulse approximation. Two body correction is isovector in
  // nature, and has both center of mass and relative components. The one body
  // component is not regularized.
  double NLO1Body(const basis::RelativeStateLSJT& bra,
                  const basis::RelativeStateLSJT& ket,
                  const util::OscillatorParameter& b);

  double NLO1Body(const basis::RelativeCMStateLSJT& bra,
                  const basis::RelativeCMStateLSJT& ket,
                  const util::OscillatorParameter& b);

  double NLO2Body(const basis::RelativeStateLSJT& bra,
                  const basis::RelativeStateLSJT& ket,
                  const util::OscillatorParameter& b,
                  const bool& regularize,
                  const double& regulator);

  double NLO2Body(const basis::RelativeCMStateLSJT& bra,
                  const basis::RelativeCMStateLSJT& ket,
                  const util::OscillatorParameter& b,
                  const bool& regularize,
                  const double& regulator);

  double M1Operator::NLOMatrixElement(const basis::RelativeStateLSJT& bra,
                                      const basis::RelativeStateLSJT& ket,
                                      const util::OscillatorParameter& b,
                                      const bool& regularize,
                                      const double& regulator)
  {
    return NLO1Body(bra, ket, b) + NLO2Body(bra, ket, b, regularize, regulator);
  }

  double M1Operator::NLOMatrixElement(const basis::RelativeCMStateLSJT& bra,
                                      const basis::RelativeCMStateLSJT &ket,
                                      const util::OscillatorParameter &b,
                                      const bool &regularize,
                                      const double &regulator)
  {
    return NLO1Body(bra, ket, b) + NLO2Body(bra, ket, b, regularize, regulator);
  }

  // The one body rme is same for both relative and relative-cm. However since
  // the NLO1Body function accepts either RelativeStateLSJT or RelativeCMStateLSJT
  // structures for bra and ket, we will first write OneBodyRME to do the actual
  // computation and use NLO1Body as a wrapper around it.
  double OneBodyRME(const std::size_t nr, const std::size_t nrp,
                    const std::size_t lr, const std::size_t lrp,
                    const std::size_t s, const std::size_t sp,
                    const std::size_t j, const std::size_t jp,
                    const std::size_t t, const std::size_t tp);

  double NLO1Body(const basis::RelativeStateLSJT& bra,
                  const basis::RelativeStateLSJT& ket,
                  const util::OscillatorParameter& b)
  {
    std::size_t nr = ket.n(), nrp = bra.n();
    std::size_t lr = ket.L(), lrp = bra.L();
    std::size_t s = ket.S(), sp = bra.S();
    std::size_t j = ket.J(), jp = bra.J();
    std::size_t t = ket.T(), tp = bra.T();

    return OneBodyRME(nr, nrp, lr, lrp, s, sp, j, jp, t, tp);
  }

  double NLO1Body(const basis::RelativeCMStateLSJT& bra,
                  const basis::RelativeCMStateLSJT& ket,
                  const util::OscillatorParameter& b)
  {
    std::size_t nr = ket.Nr(), nrp = bra.Nr();
    std::size_t nc = ket.Nc(), ncp = bra.Nc();
    std::size_t lr = ket.lr(), lrp = bra.lc();
    std::size_t lc = ket.lc(), lcp = bra.lc();
    std::size_t L = ket.L(), Lp = bra.L();
    std::size_t S = ket.S(), Sp = bra.S();
    std::size_t J = ket.J(), Jp = bra.J();
    std::size_t T = ket.T(), Tp = bra.T();

    bool kronecker = (nc == ncp && lc == lcp);
    if (!kronecker)
      return 0;
    return OneBodyRME(nr, nrp, lr, lrp, S, Sp, J, Jp, T, Tp);
  }

  double OneBodyRME(const std::size_t nr, const std::size_t nrp,
                    const std::size_t lr, const std::size_t lrp,
                    const std::size_t s, const std::size_t sp,
                    const std::size_t j, const std::size_t jp,
                    const std::size_t t, const std::size_t tp)
  {
    bool kronecker = (nr == nrp && lr == lrp);
    if (!kronecker)
      return 0;

    auto symm_term_angular = ((1 + 0.5 * std::sqrt(t * (t + 1)))
                              * ParitySign(1 + lr + s + j)
                              * Hat(lr) * std::sqrt(lr * (lr + 1))
                              * am::Wigner6J(lr, j, s, jp, lrp, 1));

    auto symm_term_spin = ((constants::isoscalar_nucleon_magnetic_moment
                            + (std::sqrt(t * (t + 1))
                               * constants::isovector_nucleon_magnetic_moment))
                               * ParitySign(lr + jp)
                               * Hat(s) * std::sqrt(s * (s + 1))
                               * am::Wigner6J(1, j, lr, jp, 1, 1));
    auto symm_term = ((symm_term_angular + symm_term_spin)
                      * (t == tp) * (s == sp));

    auto asymm_term = (0.75 * (ParitySign(t) - ParitySign(tp))
                       * (ParitySign(s) - ParitySign(sp))
                       * ParitySign(lr + s + jp + 1) * Hat(sp)
                       * am::Wigner6J(s, j, lr, jp, sp, 1)
                       * (t != tp) * (s != sp));

    auto result = Hat(j) * (symm_term + asymm_term);
    if (isnan(result))
      result = 0;
    return result;
  }

  // Relative two body part.
  double NLO2Body(const basis::RelativeStateLSJT& bra,
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

    // OscillatorParameter and scaling.
    auto brel = b.relative();
    auto scaled_regulator = regulator / brel;
    auto scaled_pion_mass = constants::pion_mass_fm / brel;

    // Parameters for integration routines.
    quadrature::gsl_params p{nr, lr, nrp, lrp, regularize, scaled_regulator, scaled_pion_mass};

    // Radial integrals.
    auto norm_product = (ho::CoordinateSpaceNorm(nr, lr, 1)
                         * ho::CoordinateSpaceNorm(nrp, lrp, 1));
    auto zpi_integral = norm_product * quadrature::IntegralZPiYPiR(p);
    auto tpi_integral = norm_product * quadrature::IntegralTPiYPiR(p);

    // Angular momentum rmes.
    auto A6S1_rme = std::sqrt(10) * am::RelativePauliProductRME(lrp, lr, sp, s, jp, j, 2, 1, 1);
    auto S1_rme = am::RelativePauliProductRME(lrp, lr, sp, s, jp, j, 0, 1, 1);

    // Isospin rme
    auto T1_rme = am::PauliProductRME(tp, t, 1);

    // LEC prefactor.
    auto num = constants::gA * util::cube(constants::pion_mass_fm) * constants::d18_fm;
    auto denom = 12 * constants::pi * util::square(constants::pion_decay_constant_fm);
    auto lec_prefactor = num / denom;
    lec_prefactor /= constants::nuclear_magneton_fm;

    // Final result.
    auto result = (lec_prefactor * T1_rme *
                   (zpi_integral * A6S1_rme + tpi_integral * S1_rme));
    if (isnan(result))
      result = 0;
    return result;
  }

  // Relative-cm two body component.
  double NLO2Body(const basis::RelativeCMStateLSJT& bra,
                  const basis::RelativeCMStateLSJT& ket,
                  const util::OscillatorParameter& b,
                  const bool& regularize,
                  const double& regulator)
  {
    std::size_t nr = ket.Nr(), nrp = bra.Nr();
    std::size_t nc = ket.Nc(), ncp = bra.Nc();
    std::size_t lr = ket.lr(), lrp = bra.lc();
    std::size_t lc = ket.lc(), lcp = bra.lc();
    std::size_t L = ket.L(), Lp = bra.L();
    std::size_t S = ket.S(), Sp = bra.S();
    std::size_t J = ket.J(), Jp = bra.J();
    std::size_t T = ket.T(), Tp = bra.T();

    // CM oscillator parameter and scaling.
    auto bcm = b.cm();
    auto scaled_regulator_cm = regulator / bcm;
    auto scaled_pion_mass_cm = constants::pion_mass_fm / bcm;

    // Relative oscillator parameter and scaling.
    auto brel = b.relative();
    auto scaled_regulator_rel = regulator / brel;
    auto scaled_pion_mass_rel = constants::pion_mass_fm / brel;

    // Parameters for integration routines.
    quadrature::gsl_params pcm{nc, lc, ncp, lcp, regularize, scaled_regulator_cm, scaled_pion_mass_cm};
    quadrature::gsl_params prel{nr, lr, nrp, lrp, regularize, scaled_regulator_rel, scaled_pion_mass_rel};

    // Radial integrals.
    // CM integral.
    auto norm_product_cm = (ho::CoordinateSpaceNorm(nc, lc, 1)
                            * ho::CoordinateSpaceNorm(ncp, lcp, 1));
    auto mpir_integral = norm_product_cm * quadrature::IntegralMPiR(pcm);
    // Relative integrals.
    auto norm_product_rel = (ho::CoordinateSpaceNorm(nr, lr, 1)
                             * ho::CoordinateSpaceNorm(nrp, lrp, 1));
    auto mpir_wpi_integral = norm_product_rel * quadrature::IntegralMPiRWPiRYPiR(prel);
    auto zpi_integral = norm_product_rel * quadrature::IntegralZPiYPiR(prel);
    auto tpi_integral = norm_product_rel * quadrature::IntegralTPiYPiR(prel);

    // Angular momentum rmes.
    auto A1_rme = (-std::sqrt(3) * am::RelativeCMPauliProductRME(lrp, lr, lcp, lc, Lp, L, Sp, S, Jp, J, 1, 1, 1, 0, 1));
    auto A2_rme = (std::sqrt(3.0 / 5.0) * am::RelativeCMPauliProductRME(lrp, lr, lcp, lc, Lp, L, Sp, S, Jp, J, 1, 1, 1, 2, 1));
    auto A3_rme = (std::sqrt(9.0 / 5.0) * am::RelativeCMPauliProductRME(lrp, lr, lcp, lc, Lp, L, Sp, S, Jp, J, 1, 1, 2, 2, 1));
    auto A4_rme = (std::sqrt(14.0 / 5.0) * am::RelativeCMPauliProductRME(lrp, lr, lcp, lc, Lp, L, Sp, S, Jp, J, 3, 1, 2, 2, 1));
    auto A5_rme = (std::sqrt(28.0 / 5.0) * am::RelativeCMPauliProductRME(lrp, lr, lcp, lc, Lp, L, Sp, S, Jp, J, 3, 1, 3, 2, 1));
    auto A6S1_rme = (std::sqrt(10) * am::RelativeCMPauliProductRME(lrp, lr, lcp, lc, Lp, L, Sp, S, Jp, J, 0, 2, 2, 1, 1));
    auto S1_rme = am::RelativeCMPauliProductRME(lrp, lr, lcp, lc, Lp, L, Sp, S, Jp, J, 0, 0, 0, 1, 1);

    // Isospin rme.
    auto T1_rme = am::PauliProductRME(Tp, T, 1);

    // LEC prefactor. (g_A m_π^3 \bar{d}_18 / 12 π F_π^2 μ_N)
    auto num = constants::gA * util::cube(constants::pion_mass_fm) * constants::d18_fm;
    auto denom = 12 * constants::pi * util::square(constants::pion_decay_constant_fm);
    auto lec_prefactor = num / denom;
    lec_prefactor /= constants::nuclear_magneton_fm;

    // Final result.
    auto api_r = A1_rme + mpir_wpi_integral * (A2_rme + A3_rme + A4_rme + A5_rme);
    auto relative_cm = mpir_integral * api_r;
    auto relative = zpi_integral * A6S1_rme + tpi_integral * S1_rme;
    auto result = lec_prefactor * T1_rme * (relative_cm + relative);
    if (isnan(result))
      result = 0;
    return result;
  }

  // Next to next to leading order. There are no chiral eft correction at N2LO.
  double M1Operator::N2LOMatrixElement(const basis::RelativeStateLSJT& bra,
                                       const basis::RelativeStateLSJT& ket,
                                       const util::OscillatorParameter& b,
                                       const bool& regularize,
                                       const double& regulator)
  {
    return 0;
  }

  double M1Operator::N2LOMatrixElement(const basis::RelativeCMStateLSJT& bra,
                                       const basis::RelativeCMStateLSJT& ket,
                                       const util::OscillatorParameter& b,
                                       const bool& regularize,
                                       const double& regulator)
  {
    return 0;
  }

  // Next to next to next to leading order.
  // There are both isoscalar and isovector two body chiral eft corrections at
  // N3LO. Currently only the isoscalar part has been implemented, as that is
  // important for the deuteron.
  double M1Operator::N3LOMatrixElement(const basis::RelativeStateLSJT& bra,
                                       const basis::RelativeStateLSJT& ket,
                                       const util::OscillatorParameter& b,
                                       const bool& regularize,
                                       const double& regulator)
  {
    return 0;
  }

  double M1Operator::N3LOMatrixElement(const basis::RelativeCMStateLSJT& bra,
                                       const basis::RelativeCMStateLSJT& ket,
                                       const util::OscillatorParameter& b,
                                       const bool& regularize,
                                       const double& regulator)
  {
    return 0;
  }

  // Next to next to next to next to leading order.
  // At present there are no results for N4LO.
  double M1Operator::N4LOMatrixElement(const basis::RelativeStateLSJT& bra,
                                       const basis::RelativeStateLSJT& ket,
                                       const util::OscillatorParameter& b,
                                       const bool& regularize,
                                       const double& regulator)
  {
    return 0;
  }

  double M1Operator::N4LOMatrixElement(const basis::RelativeCMStateLSJT& bra,
                                       const basis::RelativeCMStateLSJT& ket,
                                       const util::OscillatorParameter& b,
                                       const bool& regularize,
                                       const double& regulator)
  {
    return 0;
  }
}
