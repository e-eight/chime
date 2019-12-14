#ifndef GT_H
#define GT_H

#include "chiral.h"

using namespace factory;
namespace chiral
{
  struct GTOperator : public Operator::Registrar<GTOperator>
  {
    GTOperator() {}
    ~GTOperator() {}

    // Name of the operator, that matches the name given as input to ChiME.
    static std::string Name() { return "gt"; }

    int G0() override { return 0; } // Parity of the operator
    int J0() override { return 1; } // Angular momentum of the operator

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

  double LO1Body(const basis::RelativeStateLSJT& bra,
                 const basis::RelativeStateLSJT& ket);

  double c3Term(const basis::RelativeStateLSJT& bra,
                const basis::RelativeStateLSJT& ket,
                const ho::OscillatorParameter& b,
                const bool& regularize,
                const double& regulator);

  double c4Term(const basis::RelativeStateLSJT& bra,
                const basis::RelativeStateLSJT& ket,
                const ho::OscillatorParameter& b,
                const bool& regularize,
                const double& regulator);

  double DTerm(const basis::RelativeStateLSJT& bra,
               const basis::RelativeStateLSJT& ket,
               const ho::OscillatorParameter& b,
               const bool& regularize,
               const double& regulator);
}

#endif
