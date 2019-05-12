#ifndef UTILITY_H
#define UTILITY_H

#include <memory>
#include <map>
#include <tuple>

namespace util
{
  inline int KroneckerDelta(int a, int b) { return a == b; };

  // make_unique, because C++11 does not have one.
  // https://herbsutter.com/gotw/_102/ for details.
  template <class T, class... Args>
  std::unique_ptr<T> make_unique(Args&&... args)
  {
    return std::unique_ptr<T>( new T( std::forward<Args>(args)... ));
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
