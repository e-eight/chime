#ifndef CHARGE_RADIUS_H
#define CHARGE_RADIUS_H

#include "chiral.h"

using namespace factory;
namespace chiral
{
  struct ChargeRadiusOperator : public Operator::Registrar<ChargeRadiusOperator>
    {
      ChargeRadiusOperator() {}

      ~ChargeRadiusOperator() {}

      // You must include a Name function in your operator, and the name
      // returned by that function must match the name given as input to ChiME.
      static std::string Name() { return "rsq"; }

      int G0() override { return 0; } // Parity of the operator

      int J0() override { return 0; } // Angular momentum of the operator

      int T0() override { return 0; } // Isospin of the operator

      double LOMatrixElement(const basis::RelativeStateLSJT& bra,
                             const basis::RelativeStateLSJT& ket,
                             const double& osc_b) override;

      double NLOMatrixElement(const basis::RelativeStateLSJT& bra,
                              const basis::RelativeStateLSJT& ket,
                              const double& osc_b) override;

      double N2LOMatrixElement(const basis::RelativeStateLSJT& bra,
                               const basis::RelativeStateLSJT& ket,
                               const double& osc_b) override;

      double N3LOMatrixElement(const basis::RelativeStateLSJT& bra,
                               const basis::RelativeStateLSJT& ket,
                               const double& osc_b) override;

      double N4LOMatrixElement(const basis::RelativeStateLSJT& bra,
                               const basis::RelativeStateLSJT& ket,
                               const double& osc_b) override;
    };
}

#endif
