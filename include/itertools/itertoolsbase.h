#ifndef __ITERTOOLS_BASE_HPP__
#define __ITERTOOLS_BASE_HPP__

#include <cassert>
#include <functional>
#include <iterator>
#include <optional>
#include <tuple>
#include <type_traits>
#include <utility>

namespace itertools {

namespace getters {

// begin() for C arrays
template <typename T, std::size_t N>
T* begin_impl(T (&array)[N], int) {
  return array;
}

// Prefer member begin().
template <typename T, typename I = decltype(std::declval<T&>().begin())>
I begin_impl(T& r, int) {
  return r.begin();
}
// Use ADL otherwises.
template <typename T, typename I = decltype(begin(std::declval<T&>()))>
I begin_impl(T& r, long) {
  return std::begin(r);
}

template <typename T>
auto begin(T& t) -> decltype(get_begin_impl(std::declval<T&>(), 42)) {
  return begin_impl(t, 42);
}

// end() for C arrays
template <typename T, std::size_t N>
T* end_impl(T (&array)[N], int) {
  return array + N;
}

// Prefer member end().
template <typename T, typename I = decltype(std::declval<T&>().end())>
I end_impl(T& r, int) {
  return r.end();
}

// Use ADL otherwise.
template <typename T, typename I = decltype(end(std::declval<T&>()))>
I end_impl(T& r, long) {
  return std::end(r);
}

template <typename T>
auto end(T& t) -> decltype(get_end_impl(std::declval<T&>(), 42)) {
  return end_impl(t, 42);
}

}  // namespace getters

template <typename T>
struct type_is {
  using type = T;
};

template <typename T>
using AsConst = decltype(std::as_const(std::declval<T&>()));

template <typename T>
using iterator_type = decltype(getters::begin(std::declval<T&>()));

// iterator_type<C> is the type of C's iterator
template <typename Container>
using const_iterator_type = decltype(getters::begin(
    std::declval<const std::remove_reference_t<Container>&>()));

// const_iterator_deref is the type obtained through dereferencing
// a const iterator& (note: not a const_iterator).  ie: the result
// of Container::iterator::operator*() const
template <typename Container>
using const_iterator_deref =
    decltype(*std::declval<const iterator_type<Container>&>());

template <typename Container>
using iterator_deref = decltype(*std::declval<iterator_type<Container>&>());

// the type of dereferencing a const_iterator
template <typename Container>
using const_iterator_type_deref =
    decltype(*std::declval<const_iterator_type<Container>&>());

template <typename Container>
using iterator_traits_deref =
    std::remove_reference_t<iterator_deref<Container>>;

template <typename T, typename = void>
struct IsIterable : std::false_type {};

// Assuming that if a type works with begin, it is an iterable.
template <typename T>
struct IsIterable<T, std::void_t<iterator_type<T>>> : std::true_type {};

template <typename T>
constexpr bool is_iterable = IsIterable<T>::value;

struct Identity {
  template <typename T>
  const T& operator()(const T& t) const {
    return t;
  }
};

template <typename T, typename = void>
struct ArrowHelper {
  using type = void;
  void operator()(T&) const noexcept {}
};
template <typename T>
struct ArrowHelper<T*, void> {
  using type = T*;
  constexpr type operator()(T* t) const noexcept { return t; }
};

template <typename T>
struct ArrowHelper<T, std::void_t<decltype(std::declval<T&>().operator->())>> {
  using type = decltype(std::declval<T&>().operator->());
  type operator()(T& t) const { return t.operator->(); }
};

template <typename T>
using arrow = typename ArrowHelper<T>::type;

// type of C::iterator::operator->, also works with pointers
// void if the iterator has no operator->
template <typename C>
using iterator_arrow = arrow<iterator_type<C>>;

// applys the -> operator to an object, if the object is a pointer,
// it returns the pointer
template <typename T>
arrow<T> apply_arrow(T& t) {
  return ArrowHelper<T>{}(t);
}

// For iterators that have an operator* which returns a value
// they can return this type from their operator-> instead, which will
// wrap an object and allow it to be used with arrow
template <typename T>
class ArrowProxy {
 private:
  using TPlain = typename std::remove_reference<T>::type;
  T obj;

 public:
  constexpr ArrowProxy(T&& in_obj) : obj(std::forward<T>(in_obj)) {}

  TPlain* operator->() { return &obj; }
};

template <typename, typename = void>
struct is_random_access_iter : std::false_type {};

template <typename T>
struct is_random_access_iter<
    T, std::enable_if_t<
	   std::is_same<typename std::iterator_traits<T>::iterator_category,
			std::random_access_iterator_tag>::value>>
    : std::true_type {};

template <typename T>
using has_random_access_iter = is_random_access_iter<iterator_type<T>>;
// because std::advance assumes a lot and is actually smart, I need a dumb
// version that will work with most things

template <typename InputIt, typename Distance = std::size_t>
void dumb_advance_unsafe(InputIt& iter, Distance distance) {
  for (Distance i(0); i < distance; ++i) {
    ++iter;
  }
}

template <typename Iter, typename EndIter, typename Distance>
void dumb_advance_impl(Iter& iter, const EndIter& end, Distance distance,
		       std::false_type) {
  for (Distance i(0); i < distance && iter != end; ++i) {
    ++iter;
  }
}

template <typename Iter, typename EndIter, typename Distance>
void dumb_advance_impl(Iter& iter, const EndIter& end, Distance distance,
		       std::true_type) {
  if (static_cast<Distance>(end - iter) < distance) {
    iter = end;
  } else {
    iter += distance;
  }
}

// iter will not be incremented past end
template <typename Iter, typename EndIter, typename Distance = std::size_t>
void dumb_advance(Iter& iter, const EndIter& end, Distance distance) {
  dumb_advance_impl(iter, end, distance, is_random_access_iter<Iter>{});
}

template <typename ForwardIt, typename Distance = std::size_t>
ForwardIt dumb_next(ForwardIt it, Distance distance = 1) {
  dumb_advance_unsafe(it, distance);
  return it;
}

template <typename ForwardIt, typename Distance = std::size_t>
ForwardIt dumb_next(ForwardIt it, const ForwardIt& end, Distance distance = 1) {
  dumb_advance(it, end, distance);
  return it;
}

template <typename Container, typename Distance = std::size_t>
Distance dumb_size(Container&& container) {
  Distance d{0};
  auto end_it = getters::end(container);
  for (auto it = getters::begin(container); it != end_it; ++it) {
    ++d;
  }
  return d;
}

template <typename... Ts>
struct are_same : std::true_type {};

template <typename T, typename U, typename... Ts>
struct are_same<T, U, Ts...>
    : std::integral_constant<bool, std::is_same<T, U>::value &&
				       are_same<T, Ts...>::value> {};

// DerefHolder holds the value gotten from an iterator dereference
// if the iterate dereferences to an lvalue references, a pointer to the
// element is stored
// if it does not, a value is stored instead
//
// get() : returns a reference to the held item
// get_ptr() : returns a pointer to the held item
// reset() : replaces the currently held item

template <typename T>
class DerefHolder {
 private:
  static_assert(!std::is_lvalue_reference<T>::value,
		"Non-lvalue-ref specialization used for lvalue ref type");
  // it could still be an rvalue reference
  using TPlain = std::remove_reference_t<T>;

  std::optional<TPlain> pItem;

 public:
  using reference = TPlain&;
  using pointer = TPlain*;
  static constexpr bool stores_value = true;

  DerefHolder() = default;
  reference get() {
    assert(pItem.has_value());
    return *pItem;
  }

  pointer get_ptr() {
    assert(pItem.has_value());
    return &pItem.value();
  }
  void reset(T&& item) { pItem.emplace(std::move(item)); }

  explicit operator bool() const { return static_cast<bool>(pItem); }
};

// Specialization for when T is an lvalue ref
template <typename T>
class DerefHolder<T&> {
 public:
  using reference = T&;
  using pointer = T*;

 private:
  pointer pItem{};

 public:
  static constexpr bool stores_value = false;

  DerefHolder() = default;

  reference get() {
    assert(pItem);
    return *pItem;
  }

  pointer get_ptr() {
    assert(pItem);
    return pItem;
  }

  void reset(reference item) { pItem = &item; }

  explicit operator bool() const { return pItem != nullptr; }
};

}  // namespace itertools

#endif	//__ITERTOOLS_BASE_HPP__
