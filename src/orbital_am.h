#ifndef ORBITAL_AM_H
#define ORBITAL_AM_H

#include "chiral.h"

namespace chiral
{
  struct OrbitalAMOperator : public Operator::Registrar<OrbitalAMOperator>
  {
    OrbitalAMOperator() {}
    ~OrbitalAMOperator() {}

    static std::string Name() { return "orbital-am"; }

    int G0() override { return 0; }
    int J0() override { return 1; }

    double LOMatrixElement(const basis::RelativeStateLSJT& bra,
                             const basis::RelativeStateLSJT& ket,
                             const ho::OscillatorParameter& b,
                             const bool& regularize,
                             const double& regulator,
                             const std::size_t& T0,
                             const std::size_t& Abody) override;
      double LOMatrixElement(const basis::RelativeCMStateLSJT& bra,
                             const basis::RelativeCMStateLSJT& ket,
                             const ho::OscillatorParameter& b,
                             const bool& regularize,
                             const double& regulator,
                             const std::size_t& T0,
                             const std::size_t& Abody) override;

      double NLOMatrixElement(const basis::RelativeStateLSJT& bra,
                              const basis::RelativeStateLSJT& ket,
                              const ho::OscillatorParameter& b,
                              const bool& regularize,
                              const double& regulator,
                              const std::size_t& T0,
                              const std::size_t& Abody) override;
      double NLOMatrixElement(const basis::RelativeCMStateLSJT& bra,
                              const basis::RelativeCMStateLSJT& ket,
                              const ho::OscillatorParameter& b,
                              const bool& regularize,
                              const double& regulator,
                              const std::size_t& T0,
                              const std::size_t& Abody) override;

      double N2LOMatrixElement(const basis::RelativeStateLSJT& bra,
                               const basis::RelativeStateLSJT& ket,
                               const ho::OscillatorParameter& b,
                               const bool& regularize,
                               const double& regulator,
                               const std::size_t& T0,
                               const std::size_t& Abody) override;
      double N2LOMatrixElement(const basis::RelativeCMStateLSJT& bra,
                               const basis::RelativeCMStateLSJT& ket,
                               const ho::OscillatorParameter& b,
                               const bool& regularize,
                               const double& regulator,
                               const std::size_t& T0,
                               const std::size_t& Abody) override;

      double N3LOMatrixElement(const basis::RelativeStateLSJT& bra,
                               const basis::RelativeStateLSJT& ket,
                               const ho::OscillatorParameter& b,
                               const bool& regularize,
                               const double& regulator,
                               const std::size_t& T0,
                               const std::size_t& Abody) override;
      double N3LOMatrixElement(const basis::RelativeCMStateLSJT& bra,
                               const basis::RelativeCMStateLSJT& ket,
                               const ho::OscillatorParameter& b,
                               const bool& regularize,
                               const double& regulator,
                               const std::size_t& T0,
                               const std::size_t& Abody) override;

      double N4LOMatrixElement(const basis::RelativeStateLSJT& bra,
                               const basis::RelativeStateLSJT& ket,
                               const ho::OscillatorParameter& b,
                               const bool& regularize,
                               const double& regulator,
                               const std::size_t& T0,
                               const std::size_t& Abody) override;
      double N4LOMatrixElement(const basis::RelativeCMStateLSJT& bra,
                               const basis::RelativeCMStateLSJT& ket,
                               const ho::OscillatorParameter& b,
                               const bool& regularize,
                               const double& regulator,
                               const std::size_t& T0,
                               const std::size_t& Abody) override;
  };
}

#endif
