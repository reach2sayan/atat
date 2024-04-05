#ifndef __ATAT_FUNCTIONALTOOLS_HPP__
#define __ATAT_FUNCTIONALTOOLS_HPP__

#include <functional>
#include <memory>
#include <numeric>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

namespace ATATFunctionalTools {

template <typename F>
inline constexpr decltype(auto) compose(F&& f) {
  return std::forward<F>(f);
}

// something went wrong. This gives better error messages
template <typename P1, typename P2, typename F, typename Tail, typename... T>
void compose_impl(P1, P2, F&&, Tail&&, T&&...) {
  constexpr auto unitail = std::is_invocable_v<Tail, T...>;
  constexpr auto multitail = (std::is_invocable_v<Tail, T> && ...);
  if constexpr (unitail) {
    using tail_type = std::invoke_result_t<Tail, T...>;
    static_assert(std::is_invocable_v<F, tail_type>,
		  "Function not callable with result of next function");
  } else if constexpr (multitail) {
    static_assert(std::is_invocable_v<F, std::invoke_result<Tail, T>...>,
		  "Function not callable with results from multiple calls of "
		  "unary function");
  }
  static_assert(unitail || multitail, "function not callable");
  static_assert(sizeof...(T) == 1U || !(unitail && multitail),
		"ambigous composition");
}

// just one more at the tail, do the calculations
template <typename P, typename F, typename Tail, typename... T>
inline constexpr decltype(auto) compose_impl(std::true_type, P, F&& f,
					     Tail&& tail, T&&... objs) {
  return f(tail(std::forward<T>(objs)...));
}

// more than one-function in tail, recurse...
template <typename F, typename Tail, typename... T>
inline constexpr decltype(auto) compose_impl(std::false_type, std::true_type,
					     F&& f, Tail&& tail, T&&... objs) {
  return f(tail(std::forward<T>(objs))...);
}

template <typename F, typename... Fs>
inline constexpr decltype(auto) compose(F&& f, Fs&&... fs) {
  return [f = std::forward<F>(f), tail = compose(std::forward<Fs>(fs)...)](
	     auto&&... objs) -> decltype(auto) {
    using tail_type = decltype(tail);

    constexpr auto unitail =
	typename std::is_invocable<tail_type, decltype(objs)...>::type{};
    constexpr auto multitail =
	(std::is_invocable_v<tail_type, decltype(objs)> && ...);

    return compose_impl(unitail, std::bool_constant<multitail>{}, f, tail,
			std::forward<decltype(objs)>(objs)...);
  };
}

template <typename Predicate>
inline constexpr auto negate(Predicate&& pred) {
  return [f = std::forward<Predicate>(pred)](auto&&... obj) {
    return !f(std::forward<decltype(obj)>(obj)...);
  };
}

template <typename T>
inline constexpr auto equal(T&& t) {
  return [t = std::forward<T>(t)](const auto& obj) -> bool { return obj == t; };
}

template <typename T>
inline constexpr auto not_equal(T&& t) {
  return [t = std::forward<T>(t)](const auto& obj) -> bool { return obj != t; };
}

template <typename T>
inline constexpr auto less_than(T&& t) {
  return [t = std::forward<T>(t)](const auto& obj) -> bool { return obj < t; };
}

template <typename T>
inline constexpr auto less_than_equal(T&& t) {
  return [t = std::forward<T>(t)](const auto& obj) -> bool { return obj <= t; };
}

template <typename T>
inline constexpr auto greater_than(T&& t) {
  return [t = std::forward<T>(t)](const auto& obj) -> bool { return obj > t; };
}

template <typename T>
inline auto constexpr greater_than_equal(T&& t) {
  return [t = std::forward<T>(t)](const auto& obj) -> bool { return obj >= t; };
}

template <typename Fs, std::size_t... I, typename... T>
inline constexpr bool when_all_impl(Fs& fs, std::index_sequence<I...>,
				    const T&... t) {
  return (std::get<I>(fs)(t...) && ...);
}

template <typename... Fs>
inline constexpr auto when_all(Fs&&... fs) {
  return [funcs =
	      std::tuple(std::forward<Fs>(fs)...)](const auto&... obj) -> bool {
    return when_all_impl(funcs, std::index_sequence_for<Fs...>{}, obj...);
  };
}

template <typename Fs, std::size_t... I, typename... T>
inline constexpr bool when_any_impl(Fs& fs, std::index_sequence<I...>,
				    const T&... t) {
  return (std::get<I>(fs)(t...) || ...);
}

template <typename... Fs>
inline constexpr auto when_any(Fs&&... fs) {
  return [funcs =
	      std::tuple(std::forward<Fs>(fs)...)](const auto&... obj) -> bool {
    return when_any_impl(funcs, std::index_sequence_for<Fs...>{}, obj...);
  };
}

template <typename... Fs>
inline constexpr auto when_none(Fs&&... fs) {
  return negate(when_any(std::forward<Fs>(fs)...));
}

template <typename Predicate, typename Action>
inline constexpr auto if_then(Predicate&& predicate, Action&& action) {
  return [predicate = std::forward<Predicate>(predicate),
	  action = std::forward<Action>(action)](auto&&... obj) mutable {
    if (predicate(obj...)) action(std::forward<decltype(obj)>(obj)...);
  };
}

template <typename Predicate, typename TAction, typename FAction>
inline constexpr auto if_then_else(Predicate&& predicate, TAction&& t_action,
				   FAction&& f_action) {
  return [predicate = std::forward<Predicate>(predicate),
	  t_action = std::forward<TAction>(t_action),
	  f_action = std::forward<FAction>(f_action)](auto&&... obj) mutable
	 -> std::common_type_t<
	     decltype(t_action(std::forward<decltype(obj)>(obj)...)),
	     decltype(f_action(std::forward<decltype(obj)>(obj)...))> {
    return predicate(obj...) ? t_action(std::forward<decltype(obj)>(obj)...)
			     : f_action(std::forward<decltype(obj)>(obj)...);
  };
}

template <typename Fs, std::size_t... I, typename... T>
inline constexpr void do_all_impl(Fs& fs, std::index_sequence<I...>,
				  const T&... t) {
  ((void)(std::get<I>(fs)(t...)), ...);
}

template <typename... Fs>
inline constexpr auto do_all(Fs&&... fs) {
  return [funcs =
	      std::tuple(std::forward<Fs>(fs)...)](const auto&... obj) mutable {
    do_all_impl(funcs, std::index_sequence_for<Fs...>{}, obj...);
  };
}

}  // namespace ATATFunctionalTools

#endif	// __ATAT_FUNCTIONALTOOLS_HPP__
