#ifndef CHIRAL_H
#define CHIRAL_H

#include <string>
#include "lib/basis/lsjt_scheme.h"

namespace chiral
{
    enum class ChiralOrder { lo, nlo, n2lo, n3lo, n4lo };

    class ChiralOperator
    {
        std::string name;
        double matrix_element;

    public:
        // Constructors
        ChiralOperator() {};
        ChiralOperator(const std::string s): name{s} {}

        // Accessors
        std::string get_name() const {return name;};
        double get_matrix_element() const {return matrix_element;};
        double set_matrix_element(const ChiralOperator& O,
                                  const ChiralOrder& order,
                                  const basis::RelativeStateLSJT& bra,
                                  const basis::RelativeStateLSJT& ket,
                                  const double oscillator_energy);
    };
}

#endif
