#ifndef FACTORY_H
#define FACTORY_H

/*
 * Generic factory class based on CRTP. For details check
 * http://www.nirfriedman.com/2018/04/29/unforgettable-factory/
 */

#include <string>
#include <unordered_map>
#include "utility.h"

namespace factory
{
  template <class Base, class... Args>
  class Factory
  {
  public:
    template <class ... T>
    static std::unique_ptr<Base> make(const std::string &s, T&&... args)
    {
      return data().at(s)(std::forward<T>(args)...);
    }

    template <class T> struct Registrar : Base
    {
      friend T;

      static bool registerT()
      {
        const auto name = T::Name();
        Factory::data()[name] = [](Args... args) -> std::unique_ptr<Base>
          {
           return util::make_unique<T>(std::forward<Args>(args)...);
          };
        return true;
      }
      static bool registered;

    private:
      Registrar() : Base(Key{}) { (void)registered; }
    };

    friend Base;

  private:
    class Key
    {
      Key(){};
      template <class T> friend struct Registrar;
    };

    using FuncType = std::unique_ptr<Base> (*)(Args...);
    Factory() = default;

    using FuncMap = std::unordered_map<std::string, FuncType>;
    static FuncMap &data()
    {
      static std::unordered_map<std::string, FuncType> s;
      return s;
    }
  };

  template <class Base, class... Args>
  template <class T>
  bool Factory<Base, Args...>::Registrar<T>::registered =
    Factory<Base, Args...>::Registrar<T>::registerT();
}
#endif
