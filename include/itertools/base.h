#ifndef __ITERATORTOOLSBASE_HPP__
#define __ITERATORTOOLSBASE_HPP__

#include <algorithm>
#include <cassert>
#include <functional>
#include <optional>
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
using make_const_t = decltype(std::as_const(std::declval<T&>()));

// iterator_t<C> and const_iterator_t<C> is the type of C and const C' iterators
// respectively
template <typename Container>
using iterator_t = decltype(fancy_getters::begin(std::declval<Container&>()));

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
struct is_iterable_t : std::false_type {};

// if the type works with begin(), it is an iterable
template <typename T>
struct is_iterable_t<T, std::void_t<iterator_t<T>>> : std::true_type {};

template <typename T>
constexpr bool is_iterable_v = is_iterable_t<T>::value;

template <typename, typename = void>
struct is_iterator_t : std::false_type {};

template <typename T>
struct is_iterator_t<
    T,
    std::void_t<
	decltype(T(std::declval<const T&>())),			  // copyctor
	decltype(std::declval<T&>() = std::declval<const T&>()),  // copy =
	decltype(*std::declval<T&>()),				  // operator*
	decltype(std::declval<T&>().operator->()),		  // operator->
	decltype(++std::declval<T&>()),				  // prefix ++
	decltype(std::declval<T&>()++),				  // postfix ++
	decltype(std::declval<const T&>() != std::declval<const T&>()),	 //  !=
	decltype(std::declval<const T&>() == std::declval<const T&>())	 //  ==
	>> : std::true_type {};

template <typename T>
constexpr bool is_iterator_v = is_iterator_t<T>::value;

template <typename T>
struct is_forward_iterator_type
    : std::integral_constant<bool, is_iterator_v<T> &&
				       std::is_default_constructible_v<T>> {};

template <typename T>
struct is_move_constructible_only
    : std::integral_constant<bool, !std::is_copy_constructible_v<T> &&
				       !std::is_copy_assignable_v<T> &&
				       !std::is_move_assignable_v<T> &&
				       std::is_move_constructible_v<T>> {};

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
  // can't be explicit since returning an initializer list where a object is
  // expected is equivalent copy-assignment which can't use the copy ctor if
  // explicit
  template <typename I>
  constexpr /*explicit*/ ArrowProxy(I&& obj) : obj(std::forward<I>(obj)) {}
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

// DerefDataHolder holds the value gotten from an iterator dereference
// rather hold a pointer : if dereferences to an lvalue references
// a copy otherwise : but no worries, it is moved
// operator*() returns a reference to the held item
// operator->() returns a pointer to the held item
// reset() replaces the currently held item with a new one

template <typename T>
class DerefDataHolder {
 private:
  static_assert(!std::is_lvalue_reference<T>::value,
		"Non-lvalue-ref specialization used for lvalue ref type");
  // it could still be an rvalue reference
  using TType = std::remove_reference_t<T>;

  TType data;  // hold the real-data, could be nullptr so std::optional

 public:
  using reference = TType&;
  using pointer = TType*;

  static constexpr bool stores_value = true;

  DerefDataHolder() = default;

  reference operator*() const {
    assert(data);
    return *data;
  }

  pointer operator->() const {
    assert(data.has_value());
    return &data.value();
  }

  void reset(T&& item) { data(std::move(item)); }
  // explicit operator bool() const { return static_cast<bool>(data); }
};

// Specialization for when T is an lvalue ref
template <typename T>
class DerefDataHolder<T&> {
 public:
  using reference = T&;
  using pointer = T*;

 private:
  pointer data{};  // just store a pointer to data

 public:
  static constexpr bool stores_value = false;

  DerefDataHolder() = default;

  reference operator*() {
    assert(data);
    return *data;
  }

  pointer operator->() {
    assert(data);
    return data;
  }

  void reset(reference _data) { data = &_data; }
  explicit operator bool() const { return data != nullptr; }
};

template <typename T, typename I>
struct has_T_return : std::false_type {};

template <typename Ret, typename... Args, typename T>
struct has_T_return<Ret(Args...), T> : std::is_same<Ret, T> {};

template <typename Ret, typename... Args, typename T>
struct has_T_return<Ret (*)(Args...), T> : std::is_same<Ret, T> {};

// Makes a new function which inverts the result of the passed function (works
// only if bool
template <typename Predicate>
class Not;

template <typename Predicate>
class Not {
 private:
  Predicate predFn;

 public:
  Not(Predicate _predFn) : predFn(std::move(_predFn)) {}

  template <typename T>
  bool operator()(const T& item) {
    return !bool(std::invoke(predFn, item));
  }
  template <typename T>
  bool operator()(const T& item) const {
    return !bool(std::invoke(predFn, item));
  }
};

// Pipeable Callable generator, where ItImpl is templated on the first
// argument to the call.
template <template <typename> class ItImpl>
struct IteratorToolFnAdapter {
  template <typename T, typename... Ts>
  ItImpl<T> operator()(T&& t, Ts... ts) const {
    return {std::forward<T>(t), std::move(ts)...};
  }
};

template <template <typename, typename> class ItImpl, typename DefaultT>
struct IterToolFnOptionalBindSecond {
 private:
  // T is whatever is being held for later use
  template <typename T>
  struct FnPartial {
    mutable T stored_arg;
    constexpr FnPartial(T in_t) : stored_arg(in_t) {}

    template <typename Container>
    auto operator()(Container&& container) const {
      return IterToolFnOptionalBindSecond{}(std::forward<Container>(container),
					    stored_arg);
    }
  };

 public:
  template <typename Container, typename T>
  ItImpl<Container, T> operator()(Container&& container, T func) const {
    return {std::forward<Container>(container), std::move(func)};
  }

  template <typename T, typename = std::enable_if_t<!is_iterable_v<T>>>
  FnPartial<std::decay_t<T>> operator()(T&& func) const {
    return {std::forward<T>(func)};
  }

  template <typename Container,
	    typename = std::enable_if_t<is_iterable_v<Container>>>
  auto operator()(Container&& container) const {
    return (*this)(std::forward<Container>(container), DefaultT{});
  }
};

template <typename F>
struct PipeableAndBindFirst {
 protected:
  template <typename T>
  struct FnPartial {
    mutable T stored_arg;
    constexpr FnPartial(T in_t) : stored_arg(in_t) {}

    template <typename Container>
    auto operator()(Container&& container) const {
      return F{}(stored_arg, std::forward<Container>(container));
    }
  };

 public:
  template <typename T, typename = std::enable_if_t<!is_iterable_v<T>>>
  FnPartial<std::decay_t<T>> operator()(T&& t) const {
    return {std::forward<T>(t)};
  }
};

template <template <typename, typename> class ItImpl, typename DefaultT>
struct IterToolFnOptionalBindFirst
    : PipeableAndBindFirst<IterToolFnOptionalBindFirst<ItImpl, DefaultT>> {
 private:
  using Base =
      PipeableAndBindFirst<IterToolFnOptionalBindFirst<ItImpl, DefaultT>>;

 protected:
  template <typename Container>
  auto operator()(Container&& container, std::false_type) const {
    return static_cast<const Base&>(*this)(std::forward<Container>(container));
  }

  template <typename Container>
  auto operator()(Container&& container, std::true_type) const {
    return (*this)(DefaultT{}, std::forward<Container>(container));
  }

 public:
  template <typename T>
  auto operator()(T&& t) const {
    return (*this)(std::forward<T>(t), is_iterable_t<T>{});
  }

  template <typename T, typename Container,
	    typename = std::enable_if_t<is_iterable_v<Container>>>
  ItImpl<T, Container> operator()(T func, Container&& container) const {
    return {std::move(func), std::forward<Container>(container)};
  }
};

}  // namespace ATATIteratorTools

#endif	//__ITERATORTOOLSBASE_HPP__
