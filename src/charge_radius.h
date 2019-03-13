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

      static std::string name() { return "rsq"; }

      int G0() override { return 0; }

      int J0() override { return 0; }

      int T0() override { return 0; }

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
