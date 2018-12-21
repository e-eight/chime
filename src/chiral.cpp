#include <string>
#include "chiral.h"
#include "operator_rc_sq.h"

namespace chiral
{
    ChiralOperator::ChiralOperator():
        order(Order::full), J0(0), G0(0), T0(0) {}

    // ChiralOperator::ChiralOperator(Name name):
    //     name(name), order(Order::full), J0(0), T0(0) {}

    // ChiralOperator::ChiralOperator(Name name, Order order):
    //     name(name), order(order), J0(0), T0(0) {}

    // ChiralOperator::ChiralOperator(Name name, int J0, int T0):
    //     name(name), order(Order::full), J0(J0), T0(T0) {}

    // ChiralOperator::ChiralOperator(Name name, Order order, int J0, int T0):
    //     name(name), order(order), J0(J0), T0(T0) {}

    ChiralOperator::~ChiralOperator() {}

    ChiralOperator* ChiralOperator::create_operator(std::string name)
    {
        if (name == "rc_sq")
            return new ChargeRadiusOperator();
    }
}
