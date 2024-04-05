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

// compose takes in multiple functions, passes the input to the last function
// the output of which is passed to the penultimate function
// the final output is the output of the (last executed) first input function

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

// maps the fn to an input vector to return a new vector
template <template <typename... /*,typename*/> class OutContainer = std::vector,
	  typename Fn, class Container>
constexpr auto map(Fn mapfn, const Container& inputs) {
  using OutContainerType = typename std::decay_t<decltype(*inputs.begin())>;
  OutContainer<OutContainerType> result;
  std::transform(inputs.begin(), inputs.end(), std::back_inserter(result),
		 mapfn);
  return result;
}

// accumulates each element of a vector by some operation
template <typename Fn, class Container>
constexpr auto reduce(
    Fn reducefn, const Container& inputs,
    typename Container::value_type initval = typename Container::value_type()) {
  auto result =
      std::accumulate(inputs.begin(), inputs.end(), initval, reducefn);
  return result;
}

// takes two operations which operate on the input container
// to return a bool. Multiple Predicates are concatenated as as AND
template <typename Predicate>
constexpr auto filter(Predicate pred) {
  return [=](auto&& x) -> decltype(auto) {
    return pred(std::forward<decltype(x)>(x)) ? x : decltype(x){};
  };
}

template <typename Predicate, typename... Predicates>
constexpr auto filter(Predicate pred, Predicates... preds) {
  return [=](auto&& x) -> decltype(auto) {
    return pred(std::forward<decltype(x)>(x))
	       ? filter(preds...)(std::forward<decltype(x)>(x))
	       : decltype(x){};
  };
}
}  // namespace ATATFunctionalTools

#endif	// __ATAT_FUNCTIONALTOOLS_HPP__
