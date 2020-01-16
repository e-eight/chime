#include <cmath>
#include "basis/lsjt_scheme.h"
#include "rme_extras.h"
#include "constants.h"
#include "utility.h"
#include "integrals.h"
#include "threedho.h"
#include "spin_am.h"

namespace chiral
{
  double SpinAMOperator::LOMatrixElement(const basis::RelativeStateLSJT &bra,
                                         const basis::RelativeStateLSJT &ket,
                                         const ho::OscillatorParameter &b,
                                         const bool &regularize,
                                         const double &regulator,
                                         const std::size_t &T0,
                                         const std::size_t &Abody)
  {
    int nr = ket.n(), nrp = bra.n();
    int L = ket.L(), Lp = bra.L();
    int S = ket.S(), Sp = bra.S();
    int J = ket.J(), Jp = bra.J();
    int T = ket.T(), Tp = bra.T();

    if (am::AllowedTriangle(Jp, J, 1) && Lp == L && nrp == nr)
      {
        if (Abody != 2)
          {
            if (T0 == 0)
              {
                if (Sp == S && Tp == T)
                  {
                    double result = am::Wigner9J(L, S, J, 0, 1, 1, Lp, Sp, Jp);
                    result *= HatProduct(Lp, Sp, J, 1);
                    result *= std::sqrt(S * (S  + 1));
                    return result;
                  }
                return 0;
              }
            if (T0 == 1)
              {
                double result = Hat(S) * (Sp - S) * Hat(T) * (Tp - T);
                if (Sp == S && Tp == T)
                  result += std::sqrt(S * (S + 1)) * std::sqrt(T * (T + 1));
                result *= am::Wigner9J(L, S, J, 0, 1, 1, Lp, Sp, Jp);
                result *= HatProduct(Lp, Sp, J, 1);
                return result;
              }
            return 0;
          }
        return 0;
      }
    return 0;
  }

  double SpinAMOperator::LOMatrixElement(const basis::RelativeCMStateLSJT &bra,
                                         const basis::RelativeCMStateLSJT &ket,
                                         const ho::OscillatorParameter &b,
                                         const bool &regularize,
                                         const double &regulator,
                                         const std::size_t &T0,
                                         const std::size_t &Abody)
  {
    return 0;
  }

  double SpinAMOperator::NLOMatrixElement(const basis::RelativeStateLSJT &bra,
                                          const basis::RelativeStateLSJT &ket,
                                          const ho::OscillatorParameter &b,
                                          const bool &regularize,
                                          const double &regulator,
                                          const std::size_t &T0,
                                          const std::size_t &Abody)
  {
    return 0;
  }

  double SpinAMOperator::NLOMatrixElement(const basis::RelativeCMStateLSJT &bra,
                                          const basis::RelativeCMStateLSJT &ket,
                                          const ho::OscillatorParameter &b,
                                          const bool &regularize,
                                          const double &regulator,
                                          const std::size_t &T0,
                                          const std::size_t &Abody)
  {
    return 0;
  }

  double SpinAMOperator::N2LOMatrixElement(const basis::RelativeStateLSJT &bra,
                                           const basis::RelativeStateLSJT &ket,
                                           const ho::OscillatorParameter &b,
                                           const bool &regularize,
                                           const double &regulator,
                                           const std::size_t &T0,
                                           const std::size_t &Abody)
  {
    return 0;
  }

  double SpinAMOperator::N2LOMatrixElement(const basis::RelativeCMStateLSJT &bra,
                                           const basis::RelativeCMStateLSJT &ket,
                                           const ho::OscillatorParameter &b,
                                           const bool &regularize,
                                           const double &regulator,
                                           const std::size_t &T0,
                                           const std::size_t &Abody)
  {
    return 0;
  }

  double SpinAMOperator::N3LOMatrixElement(const basis::RelativeStateLSJT &bra,
                                           const basis::RelativeStateLSJT &ket,
                                           const ho::OscillatorParameter &b,
                                           const bool &regularize,
                                           const double &regulator,
                                           const std::size_t &T0,
                                           const std::size_t &Abody)
  {
    return 0;
  }

  double SpinAMOperator::N3LOMatrixElement(const basis::RelativeCMStateLSJT &bra,
                                           const basis::RelativeCMStateLSJT &ket,
                                           const ho::OscillatorParameter &b,
                                           const bool &regularize,
                                           const double &regulator,
                                           const std::size_t &T0,
                                           const std::size_t &Abody)
  {
    return 0;
  }

  double SpinAMOperator::N4LOMatrixElement(const basis::RelativeStateLSJT &bra,
                                           const basis::RelativeStateLSJT &ket,
                                           const ho::OscillatorParameter &b,
                                           const bool &regularize,
                                           const double &regulator,
                                           const std::size_t &T0,
                                           const std::size_t &Abody)
  {
    return 0;
  }

  double SpinAMOperator::N4LOMatrixElement(const basis::RelativeCMStateLSJT &bra,
                                           const basis::RelativeCMStateLSJT &ket,
                                           const ho::OscillatorParameter &b,
                                           const bool &regularize,
                                           const double &regulator,
                                           const std::size_t &T0,
                                           const std::size_t &Abody)
  {
    return 0;
  }
}
