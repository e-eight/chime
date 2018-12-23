#include <string>
#include "chiral.h"
#include "operator_rc_sq.h"

namespace chiral
{
    ChiralOperator::ChiralOperator():
        order(Order::full), J0(0), G0(0), T0(0) {}

    ChiralOperator::~ChiralOperator() {}
    
    ChiralOperator::operator_ptr
    ChiralOperator::create_operator(std::string name)
    {
        if (name == "rc_sq")
            return std::move(std::make_unique<ChargeRadiusOperator>());
    }
}
