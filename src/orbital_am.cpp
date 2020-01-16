#include <cmath>
#include "basis/lsjt_scheme.h"
#include "rme_extras.h"
#include "constants.h"
#include "utility.h"
#include "integrals.h"
#include "threedho.h"
#include "orbital_am.h"

namespace chiral
{
  double OrbitalAMOperator::LOMatrixElement(const basis::RelativeStateLSJT &bra,
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

    if (am::AllowedTriangle(Jp, J, 1) && Sp == S && Tp == T && Lp == L && nrp == nr)
      {
        if (T0 != 2 && Abody != 2)
          {
            double result = am::Wigner9J(L, S, J, 1, 0, 1, Lp, Sp, Jp);
            result *= HatProduct(Lp, Sp, J, 1);
            result *= std::sqrt(L * (L  + 1));
            if (T0 == 1)
              result *= 2 * std::sqrt(T * (T + 1));
            return result;
          }
        return 0;
      }
    return 0;
  }

  double OrbitalAMOperator::LOMatrixElement(const basis::RelativeCMStateLSJT &bra,
                                            const basis::RelativeCMStateLSJT &ket,
                                            const ho::OscillatorParameter &b,
                                            const bool &regularize,
                                            const double &regulator,
                                            const std::size_t &T0,
                                            const std::size_t &Abody)
  {
    return 0;
  }

  double OrbitalAMOperator::NLOMatrixElement(const basis::RelativeStateLSJT &bra,
                                             const basis::RelativeStateLSJT &ket,
                                             const ho::OscillatorParameter &b,
                                             const bool &regularize,
                                             const double &regulator,
                                             const std::size_t &T0,
                                             const std::size_t &Abody)
  {
    return 0;
  }

  double OrbitalAMOperator::NLOMatrixElement(const basis::RelativeCMStateLSJT &bra,
                                             const basis::RelativeCMStateLSJT &ket,
                                             const ho::OscillatorParameter &b,
                                             const bool &regularize,
                                             const double &regulator,
                                             const std::size_t &T0,
                                             const std::size_t &Abody)
  {
    return 0;
  }

  double OrbitalAMOperator::N2LOMatrixElement(const basis::RelativeStateLSJT &bra,
                                              const basis::RelativeStateLSJT &ket,
                                              const ho::OscillatorParameter &b,
                                              const bool &regularize,
                                              const double &regulator,
                                              const std::size_t &T0,
                                              const std::size_t &Abody)
  {
    return 0;
  }

  double OrbitalAMOperator::N2LOMatrixElement(const basis::RelativeCMStateLSJT &bra,
                                              const basis::RelativeCMStateLSJT &ket,
                                              const ho::OscillatorParameter &b,
                                              const bool &regularize,
                                              const double &regulator,
                                              const std::size_t &T0,
                                              const std::size_t &Abody)
  {
    return 0;
  }

  double OrbitalAMOperator::N3LOMatrixElement(const basis::RelativeStateLSJT &bra,
                                              const basis::RelativeStateLSJT &ket,
                                              const ho::OscillatorParameter &b,
                                              const bool &regularize,
                                              const double &regulator,
                                              const std::size_t &T0,
                                              const std::size_t &Abody)
  {
    return 0;
  }

  double OrbitalAMOperator::N3LOMatrixElement(const basis::RelativeCMStateLSJT &bra,
                                              const basis::RelativeCMStateLSJT &ket,
                                              const ho::OscillatorParameter &b,
                                              const bool &regularize,
                                              const double &regulator,
                                              const std::size_t &T0,
                                              const std::size_t &Abody)
  {
    return 0;
  }

  double OrbitalAMOperator::N4LOMatrixElement(const basis::RelativeStateLSJT &bra,
                                              const basis::RelativeStateLSJT &ket,
                                              const ho::OscillatorParameter &b,
                                              const bool &regularize,
                                              const double &regulator,
                                              const std::size_t &T0,
                                              const std::size_t &Abody)
  {
    return 0;
  }

  double OrbitalAMOperator::N4LOMatrixElement(const basis::RelativeCMStateLSJT &bra,
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
