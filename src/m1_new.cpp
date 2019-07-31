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
                    const std::size_t l, const std::size_t lp,
                    const std::size_t s, const std::size_t sp,
                    const std::size_t j, const std::size_t jp,
                    const std::size_t t, const std::size_t tp);

  double NLO1Body(const basis::RelativeStateLSJT& bra,
                  const basis::RelativeStateLSJT& ket,
                  const util::OscillatorParameter& b)
  {
    std::size_t nr = ket.n(), nrp = bra.n();
    std::size_t l = ket.L(), lp = bra.L();
    std::size_t s = ket.S(), sp = bra.S();
    std::size_t j = ket.J(), jp = bra.J();
    std::size_t t = ket.T(), tp = bra.T();

    return OneBodyRME(nr, nrp, l, lp, s, sp, j, jp, t, tp);
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
                    const std::size_t l, const std::size_t lp,
                    const std::size_t s, const std::size_t sp,
                    const std::size_t j, const std::size_t jp,
                    const std::size_t t, const std::size_t tp)
  {
    bool kronecker = (nr == nrp && l == lp);
    if (!kronecker)
      return 0;

    auto symm_term_angular = ((1 + 0.5 * std::sqrt(t * (t + 1)))
                              * ParitySign(1 + l + s + j)
                              * Hat(l) * std::sqrt(l * (l + 1))
                              * am::Wigner6J(l, j, s, jp, lp, 1));

    auto symm_term_spin = ((constants::isoscalar_nucleon_magnetic_moment
                            + (std::sqrt(t * (t + 1))
                               * constants::isovector_nucleon_magnetic_moment))
                               * ParitySign(l + jp)
                               * Hat(s) * std::sqrt(s * (s + 1))
                               * am::Wigner6J(1, j, l, jp, 1, 1));

    auto symm_term = ((symm_term_angular + symm_term_spin)
                      * (t == tp) * (s == sp));

    auto asymm_term = (0.75 * (ParitySign(t) - ParitySign(tp))
                       * (ParitySign(s) - ParitySign(sp))
                       * ParitySign(l + s + jp + 1) * Hat(sp)
                       * am::Wigner6J(s, j, l, jp, sp, 1)
                       * (t != tp) * (s != sp));

    auto result = Hat(j) * (symm_term + asymm_term);
    return result;
  }

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
    auto zpi_integral = quadrature::IntegralZPiYPiR(p);
    auto tpi_integral = quadrature::IntegralTPiYPiR(p);

    // Angular momentum rme.
    auto A6S1_rme = am::RelativePauliProductRME(lrp, lr, sp, s, jp, j, 2, 1, 1);
    auto S1_rme = am::RelativePauliProductRME(lrp, lr, sp, s, jp, j, 0, 1, 1);

    // Isospin rme
    auto T1_rme = am::PauliProductRME(tp, t, 1);

    // LEC prefactor.
    auto num = constants::gA * util::cube(constants::pion_mass_fm) * constants::d18_fm;
    auto denom = 12 * constants::pi * util::square(constants::pion_decay_constant_fm);
    auto lec_prefactor = num / denom;
    lec_prefactor /= constants::nuclear_magneton_fm;

    auto result = (lec_prefactor * T1_rme *
                   (zpi_integral * A6S1_rme + tpi_integral * S1_rme));
    return result;
  }

}
