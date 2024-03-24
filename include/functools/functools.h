#ifndef __ATAT_FUNCTIONALTOOLS_HPP__
#define __ATAT_FUNCTIONALTOOLS_HPP__

#include <algorithm>
#include <numeric>
#include <utility>
#include <vector>

namespace ATATFunctionalTools {

constexpr auto compose() {
  return
      [](auto&& x) -> decltype(auto) { return std::forward<decltype(x)>(x); };
}

// base type when this is the last function to be executed
template <typename Fn>
constexpr auto compose(Fn fn) {
  return [=](auto&&... x) -> decltype(auto) {
    return fn(std::forward<decltype(x)>(x)...);
  };
}

// more functions to be execute
template <typename Fn, typename... Fargs>
constexpr auto compose(Fn fn, Fargs... args) {
  return [=](auto&&... x) -> decltype(auto) {
    return fn(compose(args...)(std::forward<decltype(x)>(x)...));
  };
}

template <template <typename... /*,typename*/> class OutContainer = std::vector,
	  typename Fn, class Container>
constexpr auto map(Fn mapfn, const Container& inputs) {
  using OutContainerType = typename std::decay_t<decltype(*inputs.begin())>;
  OutContainer<OutContainerType> result;
  std::transform(inputs.begin(), inputs.end(), std::back_inserter(result),
		 mapfn);
  return result;
}

template <typename Fn, class Container>
constexpr auto reduce(
    Fn reducefn, const Container& inputs,
    typename Container::value_type initval = typename Container::value_type()) {
  auto result =
      std::accumulate(inputs.begin(), inputs.end(), initval, reducefn);
  return result;
}

}  // namespace ATATFunctionalTools

#endif	// __ATAT_FUNCTIONALTOOLS_HPP__
