#include <cmath>
#include "fmt/format.h"
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
                                     const double& regulator,
                                     const std::size_t& T0,
                                     const std::size_t& Abody)
  {
    return 0;
  }

  double M1Operator::LOMatrixElement(const basis::RelativeCMStateLSJT& bra,
                                     const basis::RelativeCMStateLSJT& ket,
                                     const util::OscillatorParameter& b,
                                     const bool& regularize,
                                     const double& regulator,
                                     const std::size_t& T0,
                                     const std::size_t& Abody)
  {
    return 0;
  }

  double M1Operator::NLOMatrixElement(const basis::RelativeStateLSJT& bra,
                                      const basis::RelativeStateLSJT& ket,
                                      const util::OscillatorParameter& b,
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
                                      const util::OscillatorParameter& b,
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
                                       const util::OscillatorParameter& b,
                                       const bool& regularize,
                                       const double& regulator,
                                       const std::size_t& T0,
                                       const std::size_t& Abody)
  {
    return 0;
  }

  double M1Operator::N2LOMatrixElement(const basis::RelativeCMStateLSJT& bra,
                                       const basis::RelativeCMStateLSJT& ket,
                                       const util::OscillatorParameter& b,
                                       const bool& regularize,
                                       const double& regulator,
                                       const std::size_t& T0,
                                       const std::size_t& Abody)
  {
    return 0;
  }

    double M1Operator::N3LOMatrixElement(const basis::RelativeStateLSJT& bra,
                                         const basis::RelativeStateLSJT& ket,
                                         const util::OscillatorParameter& b,
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
                                       const util::OscillatorParameter& b,
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
                                       const util::OscillatorParameter& b,
                                       const bool& regularize,
                                       const double& regulator,
                                       const std::size_t& T0,
                                       const std::size_t& Abody)
  {
    return 0;
  }

  double M1Operator::N4LOMatrixElement(const basis::RelativeCMStateLSJT& bra,
                                       const basis::RelativeCMStateLSJT& ket,
                                       const util::OscillatorParameter& b,
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
                  const util::OscillatorParameter& b,
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
    double gpi = constants::gA / constants::pion_decay_constant_fm;
    double brel = b.relative();

    if (am::AllowedTriangle(J, 1, Jp) && Sp != S && Tp != T)
      {
        if (Sp < S && Lp == Jp)
          {
            double result = ParitySign(L + Lp) * am::SphericalHarmonicCRME(Lp, L, 2);
            result *= am::Wigner6J(2, 1, 1, J, L, Jp);
            result *= quadrature::IntegralZPiYPiR(nrp, Lp, nr, L, brel, mpi, regularize, regulator);

            if (Lp == L)
              {
                double s1_term = ParitySign(J + Jp);
                s1_term *= quadrature::IntegralTPiYPiR(nrp, Lp, nr, L, brel, mpi, regularize, regulator);
                result -= s1_term;
              }

            result *= 2 * HatProduct(1, J, T); // * std::abs(T - Tp);
            result *= -mpi * mN * square(gpi) / (24 * constants::pi);
            return result;
          }
        if (Sp > S && L == J)
          {
            double result = ParitySign(J + Jp + 1) * am::SphericalHarmonicCRME(Lp, L, 2);
            result *= am::Wigner6J(1, 1, 2, J, Lp, Jp);
            result *= quadrature::IntegralZPiYPiR(nrp, Lp, nr, L, brel, mpi, regularize, regulator);

            if (Lp == L)
              {
                result += quadrature::IntegralTPiYPiR(nrp, Lp, nr, L, brel, mpi, regularize, regulator);
              }

            result *= 2 * HatProduct(1, Lp, T); // * std::abs(T - Tp);
            result *= -mpi * mN * square(gpi) / (24 * constants::pi);
            return result;
          }
        return 0;
      }
    return 0;
  }

    // // Relative oscillator parameter and scaling.
    // double brel = b.relative();
    // double scaled_regulator_rel = regulator / brel;
    // double scaled_pion_mass_rel = constants::pion_mass_fm * brel;

    // // Parameters for integration routines.
    // quadrature::gsl_params_2n
    //   prel{nrp, Lp, nr, L, regularize, scaled_regulator_rel, scaled_pion_mass_rel};

    // // Product of HO norms for integrals.
    // double norm_product_rel = (ho::CoordinateSpaceNorm(nr, L, 1)
    //                          * ho::CoordinateSpaceNorm(nrp, Lp, 1));

    // // 9j symbol common to both Uf(S1) and S1.
    // double spin_9j = am::Wigner9J(HalfInt(1, 2), HalfInt(1, 2), S, 1, 1, 1, HalfInt(1, 2), HalfInt(1, 2), Sp);

    // // Calculate the Uf(S1) term.
    // double UfS1_term = (18 * std::sqrt(10) * HatProduct(Lp, Sp, J, S));
    // UfS1_term *= (am::Wigner9J(L, S, J, 2, 1, 1, Lp, Sp, Jp) * spin_9j);
    // UfS1_term *= am::SphericalHarmonicCRME(Lp, L, 2);
    // double zpi_integral = norm_product_rel * quadrature::IntegralZPiYPiR(prel);
    // UfS1_term *= zpi_integral;

    // // Calculate the S1 term only if L and Lp are the same.
    // double S1_term = 0;
    // if (L == Lp)
    //   {
    //     S1_term = (ParitySign(S + Jp + Lp + 1) * 6 * std::sqrt(3));
    //     S1_term *= HatProduct(S, Sp, J);
    //     S1_term *= (am::Wigner6J(Sp, S, 1, J, Jp, L) * spin_9j);
    //     double tpi_integral = norm_product_rel * quadrature::IntegralTPiYPiR(prel);
    //     S1_term *= tpi_integral;
    //   }

    // // Isospin rme.
    // double T1_rme = 6 * std::sqrt(3) * Hat(T);
    // T1_rme *= am::Wigner9J(HalfInt(1, 2), HalfInt(1, 2), T, 1, 1, 1, HalfInt(1, 2), HalfInt(1, 2), Tp);

    // // LEC prefactor. (-g_A^2 m_π / 48 π F_π^2 μ_N)
    // double lecp = -(square(constants::gA) * constants::pion_mass_fm);
    // lecp /= (48 * constants::pi * constants::nuclear_magneton_fm
    //          * square(constants::pion_decay_constant_fm));

    // // Overall result.
    // double result = lecp * T1_rme * (UfS1_term + S1_term);
    // return result;

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

    if (!am::AllowedTriangle(J, Jp, 1))
      return 0;
    // Spin and isospin RMEs.
    // auto symm_rme_spin = am::RelativeCMSpinSymmetricRME(lrp, lr, lcp, lc, Lp, L, Sp, S, Jp, J, 0, 0, 0, 1);
    // auto symm_rme_isospin = am::SpinSymmetricRME(Tp, T);
    // auto asymm_rme_spin = am::RelativeCMSpinAntisymmetricRME(lrp, lr, lcp, lc, Lp, L, Sp, S, Jp, J, 0, 0, 0, 1);
    // auto asymm_rme_isospin = am::SpinAntisymmetricRME(Tp, T);
    // auto delta_T = Tp == T;

    // Orbital angular momentum MEs.
    // auto lsum_me = (am::RelativeCMLsumRME(lrp, lr, lcp, lc, Lp, L, Sp, S, Jp, J)
    //                 * (nrp == nr && ncp == nc));
    // auto mass_ratio_sqrt = 0.5; // std::sqrt(constants::reduced_nucleon_mass_fm / (2 * constants::nucleon_mass_fm));
    // auto rcm_prel_me = (mass_ratio_sqrt * am::GradientME(nrp, nr, lrp, lr)
    //                     * am::RadiusME(ncp, nc, lcp, lc));
    // auto rrel_pcm_me = (am::RadiusME(nrp, nr, lrp, lr)
    //                     * am::GradientME(ncp, nc, lcp, lc) / mass_ratio_sqrt);

    // double result = 0;

    // if (T0 == 0)
    //   {
    //     bool kronecker = (nrp == nr && lrp == lr && ncp == nc && lcp == lc && Sp == S && Tp == T);
    //     if (kronecker)
    //       {
    //         double oam_term = 0;
    //         if (am::AllowedTriangle(L, Lp, 1))
    //           {
    //             double rel_oam_term = Hat(lr) * ParitySign(Lp + J + S + lr + L + lc);
    //             rel_oam_term *= am::Wigner6J(lr, lr, 1, Lp, L, lc);
    //             rel_oam_term *= std::sqrt(lr * (lr + 1));

    //             double cm_oam_term = Hat(lc) * ParitySign(lc + lr + J + S);
    //             cm_oam_term *= am::Wigner6J(lc, lc, 1, L, Lp, lr);
    //             cm_oam_term *= std::sqrt(lc * (lc + 1));

    //             oam_term = (0.5 * HatProduct(Lp, L) * am::Wigner6J(L, Lp, 1, Jp, J, S)
    //                         * (rel_oam_term + cm_oam_term));
    //           }

    //         double spin_term = 0;
    //         if (Lp == L)
    //           {
    //             spin_term = constants::isoscalar_nucleon_magnetic_moment * Hat(S);
    //             spin_term *= (ParitySign(S + Jp + Lp + 1) * am::Wigner6J(Sp, S, 1, J, Jp, L));
    //             spin_term *= std::sqrt(S * (S  +1));
    //           }
    //         // // Purely spin term.
    //         // auto spin_term = (constants::isoscalar_nucleon_magnetic_moment
    //         //                    * symm_rme_spin * delta_T);
    //         // // Purely orbital angular momentum term.
    //         // auto oam_term = (0.5 * am::RelativeLrelRME(Lp, L, Sp, S, Jp, J) * delta_T);
    //         result = Hat(J) * (spin_term + oam_term);
    //       }
    //   }
    // else if (T0 == 1)
    //   {
    //     // Purely spin terms.
    //     // auto spin_symm_term = (constants::isovector_nucleon_magnetic_moment
    //     //                        * symm_rme_spin * symm_rme_isospin);
    //     // auto spin_asymm_term = (constants::isovector_nucleon_magnetic_moment
    //     //                         * asymm_rme_spin * asymm_rme_isospin);
    //     // // Purely orbital angular momentum term.
    //     // auto oam_diagonal_term = (0.5 * lsum_me * symm_rme_isospin);
    //     // auto oam_cross_term = (0.5 * (2 * rcm_prel_me + 0.5 * rrel_pcm_me)
    //     //                        * asymm_rme_isospin);
    //     // result = (spin_symm_term + spin_asymm_term
    //     //           + oam_diagonal_term + oam_cross_term);
    //     result = 0;
    //   }
    // else
    //   {
    //     result = 0;
    //   }
    // return result;
    return 0;
  }

  // double NLO2Body(const basis::RelativeCMStateLSJT& bra,
  //                 const basis::RelativeCMStateLSJT& ket,
  //                 const util::OscillatorParameter& b,
  //                 const bool& regularize,
  //                 const double& regulator,
  //                 const std::size_t& T0)
  // {
  //   // if (T0 != 1)
  //     return 0;

    // int nr = ket.Nr(), nrp = bra.Nr();
    // int lr = ket.lr(), lrp = bra.lr();
    // int nc = ket.Nc(), ncp = bra.Nc();
    // int lc = ket.lc(), lcp = bra.lc();
    // int L = ket.L(), Lp = bra.L();
    // int S = ket.S(), Sp = bra.S();
    // int J = ket.J(), Jp = bra.J();
    // int T = ket.T(), Tp = bra.T();

    // // CM oscillator parameter and scaling.
    // auto bcm = b.cm();
    // auto scaled_regulator_cm = regulator / bcm;
    // auto scaled_pion_mass_cm = constants::pion_mass_fm * bcm;

    // // Relative oscillator parameter and scaling.
    // auto brel = b.relative();
    // auto scaled_regulator_rel = regulator / brel;
    // auto scaled_pion_mass_rel = constants::pion_mass_fm * brel;

    // // Parameters for integration routines.
    // // quadrature::gsl_params_2n pcm{ncp, lcp, nc, lc, false, scaled_regulator_cm, scaled_pion_mass_cm};
    // quadrature::gsl_params_2n prel{nrp, lrp, nr, lr, regularize, scaled_regulator_rel, scaled_pion_mass_rel};

    // // Radial integrals.
    // // CM integral.
    // auto mpir_integral = (constants::pion_mass_fm * bcm * quadrature::IntegralMPiR(ncp, nc, lcp, lc));
    // // Relative integrals.
    // auto norm_product_rel = (ho::CoordinateSpaceNorm(nr, lr, 1)
    //                           * ho::CoordinateSpaceNorm(nrp, lrp, 1));
    // auto mpir_wpi_integral = norm_product_rel * quadrature::IntegralMPiRWPiRYPiR(prel);

    // // Angular momentum rmes.
    // auto Ua_rme = (-std::sqrt(3) * am::RelativeCMPauliProductRME(lrp, lr, lcp, lc, Lp, L, Sp, S, Jp, J, 1, 1, 1, 0, 1));
    // auto Ub_rme = (std::sqrt(3.0 / 5.0) * am::RelativeCMPauliProductRME(lrp, lr, lcp, lc, Lp, L, Sp, S, Jp, J, 1, 1, 1, 2, 1));
    // auto Uc_rme = (std::sqrt(9.0 / 5.0) * am::RelativeCMPauliProductRME(lrp, lr, lcp, lc, Lp, L, Sp, S, Jp, J, 1, 1, 2, 2, 1));
    // auto Ud_rme = (std::sqrt(14.0 / 5.0) * am::RelativeCMPauliProductRME(lrp, lr, lcp, lc, Lp, L, Sp, S, Jp, J, 3, 1, 2, 2, 1));
    // auto Ue_rme = (std::sqrt(28.0 / 5.0) * am::RelativeCMPauliProductRME(lrp, lr, lcp, lc, Lp, L, Sp, S, Jp, J, 3, 1, 3, 2, 1));

    // // Isospin rme.
    // auto T1_rme = am::PauliProductRME(Tp, T, 1);

    // // LEC prefactor. (g_A^2 m_π / 48 π F_π^2 μ_N)
    // auto lecp = (square(constants::gA) * constants::pion_mass_fm);
    // lecp /= (48 * constants::pi * constants::nuclear_magneton_fm
    //          * square(constants::pion_decay_constant_fm));

    // // Final result.
    // auto upi_r = Ua_rme + mpir_wpi_integral * (Ub_rme + Uc_rme + Ud_rme + Ue_rme);
    // auto relative_cm = mpir_integral * upi_r;
    // double relative = 0;
    // if (ncp == nc && lcp == lc)
    //   {
    //     auto UfS1_rme = (std::sqrt(10) * am::RelativePauliProductRME(lrp, lr, Sp, S, Jp, J, 2, 1, 1));
    //     auto S1_rme = am::RelativePauliProductRME(lrp, lr, Sp, S, Jp, J, 0, 1, 1);
    //     auto zpi_integral = norm_product_rel * quadrature::IntegralZPiYPiR(prel);
    //     auto tpi_integral = norm_product_rel * quadrature::IntegralTPiYPiR(prel);
    //     relative = (zpi_integral * UfS1_rme + tpi_integral * S1_rme);
    //   }
    // auto result = -lecp * T1_rme * (relative_cm + relative);
    // return result;
  }

  // Next to next to next to leading order.
  // There are both isoscalar and isovector two body chiral eft corrections at
  // N3LO. Currently only the isoscalar part has been implemented, as that is
  // important for the deuteron.

  double N3LO2BodyIsoscalar(const basis::RelativeStateLSJT& bra,
                            const basis::RelativeStateLSJT& ket,
                            const util::OscillatorParameter& b,
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
                            const util::OscillatorParameter& b,
                            const bool& regularize,
                            const double& regulator,
                            const std::size_t& T0)
  {
    return 0;
  }

}
