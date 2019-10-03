
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

  // Next to leading order reduced matrix element.
  // There are both one body and two body corrections. One body correction is
  // the same as impulse approximation. Two body correction is isovector in
  // nature, and has both center of mass and relative components. The one body
  // component is not regularized.

  // Relative NLO reduced matrix element.

  double NLO1Body(const basis::RelativeStateLSJT& bra,
                  const basis::RelativeStateLSJT& ket)
  {
    std::size_t nr = ket.n(), nrp = bra.n();
    std::size_t L = ket.L(), Lp = bra.L();
    std::size_t S = ket.S(), Sp = bra.S();
    std::size_t J = ket.J(), Jp = bra.J();
    std::size_t T = ket.T(), Tp = bra.T();

    if (nr != nrp || L != Lp)
      return 0;

    // Spin and isospin RMEs.
    auto symm_rme_spin = am::RelativeSpinSymmetricRME(Lp, L, Sp, S, Jp, J, 0, 1);
    auto symm_rme_isospin = am::SpinSymmetricRME(Tp, T);
    auto asymm_rme_spin = am::RelativeSpinAntisymmetricRME(Lp, L, Sp, S, Jp, J, 0, 1);
    auto asymm_rme_isospin = am::SpinAntisymmetricRME(Tp, T);

    // Purely spin terms.
    auto spin_symm_term = ((constants::isoscalar_nucleon_magnetic_moment * (Tp == T)
                           + constants::isovector_nucleon_magnetic_moment * symm_rme_isospin) * symm_rme_spin);
    auto spin_asymm_term = constants::isovector_nucleon_magnetic_moment * asymm_rme_isospin * asymm_rme_spin;

    // Purely orbital angular momentum term.
    auto lrel_rme = 0; //am::RelativeLrelRME(Lp, L, Sp, S, Jp, J);
    auto oam_term = (0.5 * ((Tp == T) + symm_rme_isospin) * lrel_rme);

    auto result = oam_term + spin_symm_term + spin_asymm_term;
    return result;
  }

  double NLO2Body(const basis::RelativeStateLSJT& bra,
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

    // Relative oscillator parameter and scaling.
    auto brel = b.relative();
    auto scaled_regulator_rel = regulator / brel;
    auto scaled_pion_mass_rel = constants::pion_mass_fm * brel;

    // Parameters for integration routines.
    quadrature::gsl_params_2n
      prel{nrp, Lp, nr, L, regularize, scaled_regulator_rel, scaled_pion_mass_rel};

    // Radial integrals.
    auto norm_product_rel = (ho::CoordinateSpaceNorm(nr, L, 1)
                             * ho::CoordinateSpaceNorm(nrp, Lp, 1));
    auto zpi_integral = norm_product_rel * quadrature::IntegralZPiYPiR(prel);
    auto tpi_integral = norm_product_rel * quadrature::IntegralTPiYPiR(prel);

    // Angular momentum RMEs.
    auto A6S1_rme = (std::sqrt(10) * am::RelativePauliProductRME(Lp, L, Sp, S, Jp, J, 2, 1, 1));
    auto S1_rme = am::RelativePauliProductRME(Lp, L, Sp, S, Jp, J, 0, 1, 1);

    // Isospin rme.
    auto T1_rme = am::PauliProductRME(Tp, T, 1);

    // LEC prefactor. (g_A m_π^3 \bar{d}_18 / 12 π F_π^2 μ_N)
    auto lecp = (constants::gA * constants::d18_fm
                 * cube(constants::pion_mass_fm));
    lecp /= (12 * constants::pi * constants::nuclear_magneton_fm
             * square(constants::pion_decay_constant_fm));

    // Overall result.
    auto result = A6S1_rme * zpi_integral + S1_rme * tpi_integral;
    result *= lecp * T1_rme;
    if (isnan(result))
      result = 0;
    return result;
  }

  double M1Operator::NLOMatrixElement(const basis::RelativeStateLSJT& bra,
                                      const basis::RelativeStateLSJT& ket,
                                      const util::OscillatorParameter& b,
                                      const bool& regularize,
                                      const double& regulator)
  {
    auto one_body_term = NLO1Body(bra, ket);
    auto two_body_term = NLO2Body(bra, ket, b, regularize, regulator);
    auto result = one_body_term; //+ two_body_term;
    return result;
  }

  // Relative-cm NLO reduced matrix element.
  double M1Operator::NLOMatrixElement(const basis::RelativeCMStateLSJT& bra,
                                      const basis::RelativeCMStateLSJT &ket,
                                      const util::OscillatorParameter &b,
                                      const bool &regularize,
                                      const double &regulator)
  {
    // Relative-cm quantum numbers.
    // std::size_t nr = ket.Nr(), nrp = bra.Nr();
    // std::size_t nc = ket.Nc(), ncp = bra.Nc();
    // std::size_t lr = ket.lr(), lrp = bra.lc();
    // std::size_t lc = ket.lc(), lcp = bra.lc();
    // std::size_t L = ket.L(), Lp = bra.L();
    // std::size_t S = ket.S(), Sp = bra.S();
    // std::size_t J = ket.J(), Jp = bra.J();
    // std::size_t T = ket.T(), Tp = bra.T();

    // auto one_body_term = NLO1Body(nrp, nr, lrp, lr, ncp, nc, lcp, lc,
    //                               lrp, lr, Sp, S, Jp, J, Tp, T);
    // auto two_body_term = NLO2Body(nrp, nr, lrp, lr, ncp, nc, lcp, lc, lrp, lr,
    //                               Sp, S, Jp, J, Tp, T, b, regularize, regulator);
    // auto result = one_body_term + two_body_term;
    // return result;
    return 0;
  }


  // double NLO1Body(const std::size_t& nrp, const std::size_t& nr,
  //                 const std::size_t& lrp, const std::size_t& lr,
  //                 const std::size_t& ncp, const std::size_t& nc,
  //                 const std::size_t& lcp, const std::size_t& lc,
  //                 const std::size_t& Lp, const std::size_t& L,
  //                 const std::size_t& Sp, const std::size_t& S,
  //                 const std::size_t& Jp, const std::size_t& J,
  //                 const std::size_t& Tp, const std::size_t& T)
  // {
  //   // Spin and isospin RMEs.
  //   auto symm_rme_spin = am::RelativeCMSpinSymmetricRME(lrp, lr, lcp, lc, Lp, L, Sp, S, Jp, J, 0, 0, 0, 1);
  //   auto symm_rme_isospin = am::SpinSymmetricRME(Tp, T);
  //   auto asymm_rme_spin = am::RelativeCMSpinAntisymmetricRME(lrp, lr, lcp, lc, Lp, L, Sp, S, Jp, J, 0, 0, 0, 1);
  //   auto asymm_rme_isospin = am::SpinAntisymmetricRME(Tp, T);

  //   // Orbital angular momentum MEs.
  //   auto lsum_me = am::RelativeCMLsumRME(lrp, lr, lcp, lc, Lp, L, Sp, S, Jp, J) * (nrp == nr && ncp == nc);
  //   auto mass_ratio = std::sqrt(constants::reduced_nucleon_mass_fm / (2 * constants::nucleon_mass_fm));
  //   auto rcm_prel_me = mass_ratio * am::GradientME(nrp, nr, lrp, lr) * am::RadiusME(ncp, nc, lcp, lc);
  //   auto rrel_pcm_me = am::RadiusME(nrp, nr, lrp, lr) * am::GradientME(ncp, nc, lcp, lc) / mass_ratio;

  //   // Purely orbital angular momentum terms.
  //   auto oam_diagonal_term = 0.5 * ((Tp == T) + symm_rme_isospin) * lsum_me;
  //   auto oam_cross_term = asymm_rme_isospin * (rcm_prel_me + 0.25 * rrel_pcm_me);

  //   // Purely spin terms.
  //   auto spin_symm_term = ((constants::isoscalar_nucleon_magnetic_moment * (Tp == T)
  //                          + constants::isovector_nucleon_magnetic_moment * symm_rme_isospin) * symm_rme_spin);
  //   auto spin_asymm_term = constants::isovector_nucleon_magnetic_moment * asymm_rme_isospin * asymm_rme_spin;

  //   // Complete result.
  //   auto result = oam_diagonal_term + spin_symm_term + spin_asymm_term + oam_cross_term;
  //   if (isnan(result))
  //     result = 0;
  //   return result;
  // }

  // double NLO2Body(const std::size_t& nrp, const std::size_t& nr,
  //                 const std::size_t& lrp, const std::size_t& lr,
  //                 const std::size_t& ncp, const std::size_t& nc,
  //                 const std::size_t& lcp, const std::size_t& lc,
  //                 const std::size_t& Lp, const std::size_t& L,
  //                 const std::size_t& Sp, const std::size_t& S,
  //                 const std::size_t& Jp, const std::size_t& J,
  //                 const std::size_t& Tp, const std::size_t& T,
  //                 const util::OscillatorParameter& b,
  //                 const bool& regularize, const double regulator)
  // {
  //   // CM oscillator parameter and scaling.
  //   auto bcm = b.cm();
  //   auto scaled_regulator_cm = regulator / bcm;
  //   auto scaled_pion_mass_cm = constants::pion_mass_fm * bcm;

  //   // Relative oscillator parameter and scaling.
  //   auto brel = b.relative();
  //   auto scaled_regulator_rel = regulator / brel;
  //   auto scaled_pion_mass_rel = constants::pion_mass_fm * brel;

  //   // Parameters for integration routines.
  //   quadrature::gsl_params_2n pcm{ncp, lcp, nc, lc, regularize, scaled_regulator_cm, scaled_pion_mass_cm};
  //   quadrature::gsl_params_2n prel{nrp, lrp, nr, lr, regularize, scaled_regulator_rel, scaled_pion_mass_rel};

  //   // Radial integrals.
  //   // CM integral.
  //   auto norm_product_cm = (ho::CoordinateSpaceNorm(nc, lc, 1)
  //                           * ho::CoordinateSpaceNorm(ncp, lcp, 1));
  //   auto mpir_integral = norm_product_cm * quadrature::IntegralMPiR(pcm);
  //   // Relative integrals.
  //   auto norm_product_rel = (ho::CoordinateSpaceNorm(nr, lr, 1)
  //                            * ho::CoordinateSpaceNorm(nrp, lrp, 1));
  //   auto mpir_wpi_integral = norm_product_rel * quadrature::IntegralMPiRWPiRYPiR(prel);
  //   auto zpi_integral = norm_product_rel * quadrature::IntegralZPiYPiR(prel);
  //   auto tpi_integral = norm_product_rel * quadrature::IntegralTPiYPiR(prel);

  //   // Angular momentum rmes.
  //   auto A1_rme = (-std::sqrt(3) * am::RelativeCMPauliProductRME(lrp, lr, lcp, lc, Lp, L, Sp, S, Jp, J, 1, 1, 1, 0, 1));
  //   auto A2_rme = (std::sqrt(3.0 / 5.0) * am::RelativeCMPauliProductRME(lrp, lr, lcp, lc, Lp, L, Sp, S, Jp, J, 1, 1, 1, 2, 1));
  //   auto A3_rme = (std::sqrt(9.0 / 5.0) * am::RelativeCMPauliProductRME(lrp, lr, lcp, lc, Lp, L, Sp, S, Jp, J, 1, 1, 2, 2, 1));
  //   auto A4_rme = (std::sqrt(14.0 / 5.0) * am::RelativeCMPauliProductRME(lrp, lr, lcp, lc, Lp, L, Sp, S, Jp, J, 3, 1, 2, 2, 1));
  //   auto A5_rme = (std::sqrt(28.0 / 5.0) * am::RelativeCMPauliProductRME(lrp, lr, lcp, lc, Lp, L, Sp, S, Jp, J, 3, 1, 3, 2, 1));
  //   auto A6S1_rme = (std::sqrt(10) * am::RelativeCMPauliProductRME(lrp, lr, lcp, lc, Lp, L, Sp, S, Jp, J, 0, 2, 2, 1, 1));
  //   auto S1_rme = am::RelativeCMPauliProductRME(lrp, lr, lcp, lc, Lp, L, Sp, S, Jp, J, 0, 0, 0, 1, 1);

  //   // Isospin rme.
  //   auto T1_rme = am::PauliProductRME(Tp, T, 1);

  //   // LEC prefactor. (g_A m_π^3 \bar{d}_18 / 12 π F_π^2 μ_N)
  //   auto num = constants::gA * util::cube(constants::pion_mass_fm) * constants::d18_fm;
  //   auto denom = 12 * constants::pi * util::square(constants::pion_decay_constant_fm);
  //   auto lec_prefactor = num / denom;
  //   lec_prefactor /= constants::nuclear_magneton_fm;

  //   // Final result.
  //   auto api_r = A1_rme + mpir_wpi_integral * (A2_rme + A3_rme + A4_rme + A5_rme);
  //   auto relative_cm = mpir_integral * api_r;
  //   auto relative = zpi_integral * A6S1_rme + tpi_integral * S1_rme;
  //   auto result = lec_prefactor * T1_rme * (relative_cm + relative);
  //   if (isnan(result))
  //     result = 0;
  //   return result;
  // }


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

  double N3LO2BodyIsoscalar(const basis::RelativeStateLSJT& bra,
                            const basis::RelativeStateLSJT& ket,
                            const util::OscillatorParameter& b,
                            const bool& regularize,
                            const double regulator)
  {
    // Relative quantum numbers.
    std::size_t nr = ket.n(), nrp = bra.n();
    std::size_t L = ket.L(), Lp = bra.L();
    std::size_t S = ket.S(), Sp = bra.S();
    std::size_t J = ket.J(), Jp = bra.J();
    std::size_t T = ket.T(), Tp = bra.T();

    // Spin rmes.
    auto S_rme = am::RelativeSpinSymmetricRME(Lp, L, Sp, S, Jp, J, 0, 1);

    // GSL parameters for radial integrals.
    auto brel = b.relative();
    auto scaled_regulator_rel = regulator / brel;
    auto scaled_pion_mass_rel = constants::pion_mass_fm * brel;
    quadrature::gsl_params_2n prel{nrp, Lp, nr, L, regularize, scaled_pion_mass_rel, scaled_regulator_rel};

    // d9 term.
    // d9 isospin rme.
    auto T0_rme = am::PauliProductRME(Tp, T, 0);
    // d9 radial integrals.
    auto norm_product = (ho::CoordinateSpaceNorm(nr, L, 1)
                         * ho::CoordinateSpaceNorm(nrp, Lp, 1));
    auto ypi_integral = norm_product * quadrature::IntegralYPiR(prel);
    auto wpi_integral = norm_product * quadrature::IntegralWPiRYPiR(prel);
    // d9  angular momentum rmes.
    auto A6S_rme = (std::sqrt(10) * am::RelativeSpinSymmetricRME(Lp, L, Sp, S, Jp, J, 2, 1));
    // d9 prefactor.
    auto d9_prefactor = (constants::gA * constants::d9_fm
                         * cube(constants::pion_mass_fm));
    d9_prefactor /= (std::sqrt(3) * constants::pi
                     * square(constants::pion_decay_constant_fm));
    auto d9_term = (d9_prefactor * T0_rme
                    * (wpi_integral * A6S_rme - ypi_integral * S_rme));

    // L2 term.
    double L2_term = 0;
    if (L == 0 && L == 0 && Tp == T)
      {
        auto delta_integral = quadrature::IntegralRegularizedDelta(prel) / cube(brel);
        L2_term += (2 * constants::L2_fm * S_rme * delta_integral);
      }

    // Overall result.
    auto result = d9_term + L2_term;
    result *= 2 * constants::nucleon_mass_fm;
    if (isnan(result))
      result = 0;
    return result;
  }

  double M1Operator::N3LOMatrixElement(const basis::RelativeStateLSJT& bra,
                                       const basis::RelativeStateLSJT& ket,
                                       const util::OscillatorParameter& b,
                                       const bool& regularize,
                                       const double& regulator)
  {
    auto result = N3LO2BodyIsoscalar(bra, ket, b, regularize, regulator);
    return result;
  }

  double M1Operator::N3LOMatrixElement(const basis::RelativeCMStateLSJT& bra,
                                       const basis::RelativeCMStateLSJT& ket,
                                       const util::OscillatorParameter& b,
                                       const bool& regularize,
                                       const double& regulator)
  {
    // Relative-cm quantum numbers.
    // std::size_t nr = ket.Nr(), nrp = bra.Nr();
    // std::size_t nc = ket.Nc(), ncp = bra.Nc();
    // std::size_t lr = ket.lr(), lrp = bra.lc();
    // std::size_t lc = ket.lc(), lcp = bra.lc();
    // std::size_t L = ket.L(), Lp = bra.L();
    // std::size_t S = ket.S(), Sp = bra.S();
    // std::size_t J = ket.J(), Jp = bra.J();
    // std::size_t T = ket.T(), Tp = bra.T();

    // return N3LO2BodyIsoscalar(nrp, nr, lrp, lr, ncp, nc, lcp, lc, Lp, L,
    //                           Sp, S, Jp, J, Tp, T, b, regularize, regulator);
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
