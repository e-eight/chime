#ifndef M1_H
#define M1_H

#include "chiral.h"

using namespace factory;
namespace chiral
{
  struct M1Operator : public Operator::Registrar<M1Operator>
    {
      M1Operator() {}

      ~M1Operator() {}

      // You must include a Name function in your operator, and the name
      // returned by that function must match the name given as input to ChiME.
      static std::string Name() { return "m1"; }

      int G0() override { return 0; } // Parity of the operator

      int J0() override { return 1; } // Angular momentum of the operator

      int T0() override { return 1; } // Isospin of the operator

      double LOMatrixElement(const basis::RelativeStateLSJT& bra,
                             const basis::RelativeStateLSJT& ket,
                             const util::OscillatorParameter& b,
                             const bool& regularize,
                             const double& regulator) override;
      double LOMatrixElement(const basis::RelativeCMStateLSJT& bra,
                             const basis::RelativeCMStateLSJT& ket,
                             const util::OscillatorParameter& b,
                             const bool& regularize,
                             const double& regulator) override;

      double NLOMatrixElement(const basis::RelativeStateLSJT& bra,
                              const basis::RelativeStateLSJT& ket,
                              const util::OscillatorParameter& b,
                              const bool& regularize,
                              const double& regulator) override;
      double NLOMatrixElement(const basis::RelativeCMStateLSJT& bra,
                              const basis::RelativeCMStateLSJT& ket,
                              const util::OscillatorParameter& b,
                              const bool& regularize,
                              const double& regulator) override;

      double N2LOMatrixElement(const basis::RelativeStateLSJT& bra,
                               const basis::RelativeStateLSJT& ket,
                               const util::OscillatorParameter& b,
                               const bool& regularize,
                               const double& regulator) override;
      double N2LOMatrixElement(const basis::RelativeCMStateLSJT& bra,
                               const basis::RelativeCMStateLSJT& ket,
                               const util::OscillatorParameter& b,
                               const bool& regularize,
                               const double& regulator) override;

      double N3LOMatrixElement(const basis::RelativeStateLSJT& bra,
                               const basis::RelativeStateLSJT& ket,
                               const util::OscillatorParameter& b,
                               const bool& regularize,
                               const double& regulator) override;
      double N3LOMatrixElement(const basis::RelativeCMStateLSJT& bra,
                               const basis::RelativeCMStateLSJT& ket,
                               const util::OscillatorParameter& b,
                               const bool& regularize,
                               const double& regulator) override;

      double N4LOMatrixElement(const basis::RelativeStateLSJT& bra,
                               const basis::RelativeStateLSJT& ket,
                               const util::OscillatorParameter& b,
                               const bool& regularize,
                               const double& regulator) override;
      double N4LOMatrixElement(const basis::RelativeCMStateLSJT& bra,
                               const basis::RelativeCMStateLSJT& ket,
                               const util::OscillatorParameter& b,
                               const bool& regularize,
                               const double& regulator) override;
    };

  double NLO1Body(const std::size_t& nrp, const std::size_t& nr,
                  const std::size_t& lrp, const std::size_t& lr,
                  const std::size_t& ncp, const std::size_t& nc,
                  const std::size_t& lcp, const std::size_t& lc,
                  const std::size_t& Lp, const std::size_t& L,
                  const std::size_t& Sp, const std::size_t& S,
                  const std::size_t& Jp, const std::size_t& J,
                  const std::size_t& Tp, const std::size_t& T);

  double NLO2Body(const std::size_t& nrp, const std::size_t& nr,
                  const std::size_t& lrp, const std::size_t& lr,
                  const std::size_t& ncp, const std::size_t& nc,
                  const std::size_t& lcp, const std::size_t& lc,
                  const std::size_t& Lp, const std::size_t& L,
                  const std::size_t& Sp, const std::size_t& S,
                  const std::size_t& Jp, const std::size_t& J,
                  const std::size_t& Tp, const std::size_t& T,
                  const util::OscillatorParameter& b,
                  const bool& regularize, const double regulator);

  double N3LO2BodyIsoscalar(const std::size_t& nrp, const std::size_t& nr,
                            const std::size_t& lrp, const std::size_t& lr,
                            const std::size_t& ncp, const std::size_t& nc,
                            const std::size_t& lcp, const std::size_t& lc,
                            const std::size_t& Lp, const std::size_t& L,
                            const std::size_t& Sp, const std::size_t& S,
                            const std::size_t& Jp, const std::size_t& J,
                            const std::size_t& Tp, const std::size_t& T,
                            const util::OscillatorParameter& b,
                            const bool& regularize, const double regulator);
}

#endif
