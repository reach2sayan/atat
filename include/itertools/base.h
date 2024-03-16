#ifndef __ITERATORTOOLSBASE_HPP__
#define __ITERATORTOOLSBASE_HPP__

#include <algorithm>
#include <functional>
#include <tuple>
#include <type_traits>
#include <utility>

namespace ATATIteratorTools {
namespace fancy_getters {
template <typename T, std::size_t N>
T* begin_impl(T (&array)[N], int) {
  return array;
}
template <typename T, typename I = decltype(std::declval<T&>().begin())>
I begin_impl(T& r, int) {
  return r.begin();
}

template <typename T>
auto begin(T& t) -> decltype(begin_impl(std::declval<T&>(), 42)) {
  return begin_impl(t, 42);
}

template <typename T, std::size_t N>
T* end_impl(T (&array)[N], int) {
  return array + N;
}
template <typename T, typename I = decltype(std::declval<T&>().end())>
I end_impl(T& r, int) {
  return r.end();
}

template <typename T>
auto end(T& t) -> decltype(end_impl(std::declval<T&>(), 42)) {
  return end_impl(t, 42);
}
}  // namespace fancy_getters

template <typename T>
using make_const_t = std::add_const_t<decltype(std::declval<T&>())>;

// iterator_t<C> and const_iterator_t<C> is the type of C and const C' iterators
// respectively
template <typename T>
using iterator_t = decltype(fancy_getters::begin(std::declval<T&>()));

template <typename Container>
using const_iterator_t = decltype(fancy_getters::begin(
    std::declval<const std::remove_reference_t<Container>&>()));

template <typename Container>
using iterator_deref = decltype(*std::declval<iterator_t<Container>&>());

// const_iterator_deref is the type obtained through dereferencing
// const iterator&  (not a const_iterator)
// ie: the result of Container::iterator::operator*() const
template <typename Container>
using const_iterator_deref =
    decltype(*std::declval<const_iterator_t<Container>&>());

// the type of dereferencing a const_iterator
template <typename Container>
using const_iterator_t_deref =
    decltype(*std::declval<const_iterator_t<Container>&>());

template <typename Container>
using iterator_traits_deref =
    std::remove_reference_t<iterator_deref<Container>>;

// iterator_end_t
template <typename Container>
using iterator_end_t = decltype(fancy_getters::end(std::declval<Container&>()));

template <typename T, typename = void>
struct is_iterable_type : std::false_type {};

// if the type works with begin(), it is an iterable
template <typename T>
struct is_iterable_type<T, std::void_t<iterator_t<T>>> : std::true_type {};

template <typename T>
constexpr bool is_iterable_t = is_iterable_type<T>::value;

template <typename... Ts>
struct are_same : std::true_type {};

template <typename T, typename U, typename... Ts>
struct are_same<T, U, Ts...>
    : std::integral_constant<bool, std::is_same<T, U>::value &&
				       are_same<T, Ts...>::value> {};

template <typename T, typename = void>
struct ArrowOperatorType {
  using type = void;
  void operator()(T&) const noexcept {}
};  // if no operator->()

template <typename T>
struct ArrowOperatorType<T*, void> {
  using type = T*;
  constexpr type operator()(T* t) const noexcept { return t; }
};  // if obj is pointer

template <typename T>
struct ArrowOperatorType<
    T, std::void_t<decltype(std::declval<T&>().operator->())>> {
  using type = decltype(std::declval<T&>().operator->());
  type operator()(T& t) const { return t.operator->(); }
};  // if it has a operator->()

template <typename T>
using arrow_t = typename ArrowOperatorType<T>::type;

// type of C::iterator::operator->
// void if no such operator
template <typename Container>
using iterator_arrow_t = arrow_t<iterator_t<Container>>;

// apply operator-> to an object
// returns the pointer if
template <typename T>
arrow_t<T> apply_arrow(T& t) {
  return ArrowOperatorType<T>{}(t);
}

// This is used a return-type object to hold the value which was
// returned by the iterator operator*()
template <typename T>
class ArrowProxy {
 public:
  template <typename I>
  constexpr explicit ArrowProxy(I&& obj) : obj(std::forward<I>(obj)) {}
  std::remove_reference_t<T>* operator->() { return &obj; }

 private:
  T obj;
};

// the 'evaluators' only serve the purpose to give the current Types in
// declval(expr)
template <typename... Ts>
std::tuple<iterator_deref<Ts>...> iterator_tuple_deref_evaluator(
    const std::tuple<Ts...>&);

template <typename... Ts>
std::tuple<iterator_t<Ts>...> iterator_tuple_type_evaluator(
    const std::tuple<Ts...>&);

template <typename... Ts>
std::tuple<iterator_deref<const std::remove_reference_t<Ts>>...>
const_iterator_tuple_deref_evaluator(const std::tuple<Ts...>&);

template <typename... Ts>
std::tuple<iterator_t<const std::remove_reference_t<Ts>>...>
const_iterator_tuple_type_evaluator(const std::tuple<Ts...>&);

// Given a tuple-type template argument, gives the type of the tuple of
// iterators for the tuples contained types.
template <typename TupleType>
using iterator_tuple_t =
    decltype(iterator_tuple_type_evaluator(std::declval<TupleType>()));

template <typename TupleType>
using const_iterator_tuple_t =
    decltype(const_iterator_tuple_type_evaluator(std::declval<TupleType>()));

// Given a tuple-type template argument, gives the type of the object
// which the iterators pf the contained type dereferemce to
template <typename TupleType>
using iterator_deref_tuple =
    decltype(iterator_tuple_deref_evaluator(std::declval<TupleType>()));

template <typename TupleType>
using const_iterator_deref_tuple =
    decltype(const_iterator_tuple_deref_evaluator(std::declval<TupleType>()));

// function absorbing all arguments passed to it. used when
// applying a function to a parameter pack but not passing the evaluated
// results anywhere
template <typename... Ts>
inline void blackhole(Ts&&...) {}
}  // namespace ATATIteratorTools

#endif	//__ITERATORTOOLSBASE_HPP__
