#ifndef UTILITY_H
#define UTILITY_H

#include <cmath>
#include <memory>
#include <map>
#include <tuple>
#include "constants.h"

namespace util
{
  // make_unique, because C++11 does not have one.
  // https://herbsutter.com/gotw/_102/ for details.
  template <class T, class... Args>
  std::unique_ptr<T> make_unique(Args&&... args)
  {
    return std::unique_ptr<T>( new T( std::forward<Args>(args)... ));
  }


  // Returns an array of type T, and length size.
  template <class T, std::size_t size>
  static inline std::array<T, size> GetArray()
  {
    std::array<T, size> a;
    return a;
  }


  // Oscillator parameter for relative & relative-cm cases.
  struct OscillatorParameter
  {
    OscillatorParameter()
      : oscillator_energy(0) {}
    OscillatorParameter(double _energy)
      : oscillator_energy(_energy) {}
    ~OscillatorParameter() = default;

    double relative() const
    {
      auto b = constants::hbarc;
      b /= std::sqrt(constants::reduced_nucleon_mass_MeV * oscillator_energy);
      return b;
    }
    double cm() const
    {
      auto b = constants::hbarc;
      b /= std::sqrt(2 * constants::nucleon_mass_MeV * oscillator_energy);
      return b;
    }
  private:
    double oscillator_energy; // in MeV
  };


  // Square function.
  template <class T>
  T square(const T a)
  {
    return a * a;
  }
  // Cube function.
  template <class T>
  T cube(const T a)
  {
    return a * a * a;
  }


  // Lenpic coordinate space semi local regulator.
  static inline double LenpicSemiLocalRegulator(const double r, const double R)
  {
    return std::pow(1 - std::exp(-square(r / R)), 6);
  }
  static inline double ApplyRegulator(const bool regularize, const double r, const double R)
  {
    return regularize ? LenpicSemiLocalRegulator(r, R) : 1;
  }


  // Yukawa function that appears in many operator matrix elements.
  static inline double YPi(const double scaled_r, const double scaled_pion_mass)
  {
    auto rpi = scaled_pion_mass * scaled_r;
    return std::exp(-rpi) / rpi;
  }

  // W_\pi, T_\pi, and  Z_\pi appearing in multiple operator matrix elements.
  static inline double ZPi(const double scaled_r, const double scaled_pion_mass)
  {
    return 1 + scaled_pion_mass * scaled_r;
  }
  static inline double TPi(const double scaled_r, const double scaled_pion_mass)
  {
    return 2 * ZPi(scaled_r, scaled_pion_mass) - 1;
  }
  static inline double WPi(const double scaled_r, const double scaled_pion_mass)
  {
    auto rpi2 = square(scaled_pion_mass * scaled_r);
    return 1 + (ZPi(scaled_r, scaled_pion_mass) / rpi2);
  }


  // General memoizer.
  // https://stackoverflow.com/a/17807129, and
  // http://slackito.com/2011/03/17/automatic-memoization-in-cplusplus0x/
  // for details.
  template <class T, class... Args>
  auto memoize(T (*fn)(Args...))
  {
    std::map<std::tuple<Args...>, T> table;
    return [fn, table](Args... args) mutable -> T
           {
             auto index = std::make_tuple(args...);
             auto found = table.find(index);
             if (found == table.end())
               {
                 auto result = fn(args...);
                 table[index] = fn(args...);
                 return result;
               }
             return found->second;
           };
  }
}
#endif
