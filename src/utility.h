#ifndef UTILITY_H
#define UTILITY_H

#include <cmath>
#include <memory>
#include <map>
#include <tuple>

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

  // Lenpic coordinate space semi local regulator.
  static inline double LenpicSemiLocalRegulator(const double r, const double R)
  {
    return std::pow(1 - std::exp(-r * r / R / R), 6);
  }

  // T_\pi and Z_\pi appearing in multiple operator matrix elements.
  static inline double ZPi(const double y, const double scale)
  {
    double inv_scaled_y = 1.0 / (scale * std::sqrt(y));
    return (inv_scaled_y * (inv_scaled_y + 1.0));
  }
  static inline double TPi(const double y, const double scale)
  {
    return (1 + 3 * ZPi(y, scale));
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
