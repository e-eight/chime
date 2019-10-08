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

      double LOMatrixElement(const basis::RelativeStateLSJT& bra,
                             const basis::RelativeStateLSJT& ket,
                             const util::OscillatorParameter& b,
                             const bool& regularize,
                             const double& regulator,
                             const std::size_t& T0) override;
      double LOMatrixElement(const basis::RelativeCMStateLSJT& bra,
                             const basis::RelativeCMStateLSJT& ket,
                             const util::OscillatorParameter& b,
                             const bool& regularize,
                             const double& regulator,
                             const std::size_t& T0) override;

      double NLOMatrixElement(const basis::RelativeStateLSJT& bra,
                              const basis::RelativeStateLSJT& ket,
                              const util::OscillatorParameter& b,
                              const bool& regularize,
                              const double& regulator,
                              const std::size_t& T0) override;
      double NLOMatrixElement(const basis::RelativeCMStateLSJT& bra,
                              const basis::RelativeCMStateLSJT& ket,
                              const util::OscillatorParameter& b,
                              const bool& regularize,
                              const double& regulator,
                              const std::size_t& T0) override;

      double N2LOMatrixElement(const basis::RelativeStateLSJT& bra,
                               const basis::RelativeStateLSJT& ket,
                               const util::OscillatorParameter& b,
                               const bool& regularize,
                               const double& regulator,
                               const std::size_t& T0) override;
      double N2LOMatrixElement(const basis::RelativeCMStateLSJT& bra,
                               const basis::RelativeCMStateLSJT& ket,
                               const util::OscillatorParameter& b,
                               const bool& regularize,
                               const double& regulator,
                               const std::size_t& T0) override;

      double N3LOMatrixElement(const basis::RelativeStateLSJT& bra,
                               const basis::RelativeStateLSJT& ket,
                               const util::OscillatorParameter& b,
                               const bool& regularize,
                               const double& regulator,
                               const std::size_t& T0) override;
      double N3LOMatrixElement(const basis::RelativeCMStateLSJT& bra,
                               const basis::RelativeCMStateLSJT& ket,
                               const util::OscillatorParameter& b,
                               const bool& regularize,
                               const double& regulator,
                               const std::size_t& T0) override;

      double N4LOMatrixElement(const basis::RelativeStateLSJT& bra,
                               const basis::RelativeStateLSJT& ket,
                               const util::OscillatorParameter& b,
                               const bool& regularize,
                               const double& regulator,
                               const std::size_t& T0) override;
      double N4LOMatrixElement(const basis::RelativeCMStateLSJT& bra,
                               const basis::RelativeCMStateLSJT& ket,
                               const util::OscillatorParameter& b,
                               const bool& regularize,
                               const double& regulator,
                               const std::size_t& T0) override;

      double BareMatrixElement(const basis::RelativeStateLSJT& bra,
                               const basis::RelativeStateLSJT& ket,
                               const util::OscillatorParameter& b,
                               const std::size_t& T0) override;
      double BareMatrixElement(const basis::RelativeCMStateLSJT& bra,
                               const basis::RelativeCMStateLSJT& ket,
                               const util::OscillatorParameter& b,
                               const std::size_t& T0) override;
    };

  double NLO1Body(const basis::RelativeStateLSJT& bra,
                  const basis::RelativeStateLSJT& ket,
                  const std::size_t& T0);
  double NLO1Body(const basis::RelativeCMStateLSJT& bra,
                  const basis::RelativeCMStateLSJT& ket,
                  const std::size_t& T0);
  double NLO2Body(const basis::RelativeStateLSJT& bra,
                  const basis::RelativeStateLSJT& ket,
                  const util::OscillatorParameter& b,
                  const bool& regularize,
                  const double& regulator,
                  const std::size_t& T0);
  double NLO2Body(const basis::RelativeCMStateLSJT& bra,
                  const basis::RelativeCMStateLSJT& ket,
                  const util::OscillatorParameter& b,
                  const bool& regularize,
                  const double& regulator,
                  const std::size_t& T0);
  double N3LO2BodyIsoscalar(const basis::RelativeStateLSJT& bra,
                            const basis::RelativeStateLSJT& ket,
                            const util::OscillatorParameter& b,
                            const bool& regularize,
                            const double& regulator,
                            const std::size_t& T0);
  double N3LO2BodyIsoscalar(const basis::RelativeCMStateLSJT& bra,
                            const basis::RelativeCMStateLSJT& ket,
                            const util::OscillatorParameter& b,
                            const bool& regularize,
                            const double& regulator,
                            const std::size_t& T0);
}

#endif
