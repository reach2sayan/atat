#ifndef __ITERATORTOOLS_ENUMERATOR_HPP__
#define __ITERATORTOOLS_ENUMERATOR_HPP__

#include "base.h"

namespace ATATIteratorTools {
// Held by the Enumerable::Iterator.  Has a .index, and a
// .element referencing the value yielded by the subiterator
template <typename Elem>
class EnumeratorDataHolder : public std::pair<std::size_t, Elem> {
  using std::pair<std::size_t, Elem>::pair;

 public:
  typename std::pair<std::size_t, Elem>::first_type index =
      std::pair<std::size_t, Elem>::first;
  typename std::pair<std::size_t, Elem>::second_type element =
      std::pair<std::size_t, Elem>::second;
};

template <typename Container>
class Enumerable;

// We any need a type for the index since we have to specify it for std::pair
// so why not allow the enumerate to choose it's index type, instead of
// defaulting to some unsigned integer type
struct EnumeratorClosureObject {
 private:
  template <typename Container>
  constexpr Enumerable<Container> operator()(Container&& container,
					     std::size_t start_index) const {
    return {std::forward<Container>(container), std::move(start_index)};
  }

 public:
  template <typename Container,
	    typename = std::enable_if_t<is_iterable_v<Container>>>
  constexpr auto operator()(Container&& container) const {
    return (*this)(std::forward<Container>(container), std::size_t{});
  }
};
constexpr EnumeratorClosureObject enumerate{};

}  // namespace ATATIteratorTools

// Partial specialization of std::tuple related templates for data-holder
namespace std {
template <typename T>
struct tuple_size<ATATIteratorTools::EnumeratorDataHolder<T>>
    : public tuple_size<std::pair<std::size_t, T>> {};

template <std::size_t N, typename T>
struct tuple_element<N, ATATIteratorTools::EnumeratorDataHolder<T>>
    : public tuple_element<N, std::pair<std::size_t, T>> {};
}  // namespace std

template <typename Container>
class ATATIteratorTools::Enumerable {
 private:
  Container container_;
  const std::size_t start_;

  friend EnumeratorClosureObject;

 protected:
  // private Value constructor
  constexpr Enumerable(Container&& container, std::size_t start)
      : container_(std::forward<Container>(container)), start_{start} {}

 public:
  constexpr Enumerable(Enumerable&&) = default;

  template <typename CContainer>
  using Data = EnumeratorDataHolder<iterator_deref_t<CContainer>>;

  //  Holds an iterator of the contained type and an Index for the
  // index_.  Each call to ++ increments both of these data members.
  // Each dereference returns an IterYield.
  template <typename ContainerT>
  class Iterator {
   private:
    template <typename>
    friend class Iterator;

    iterator_t<ContainerT> sub_iter_;
    std::size_t index_;

   public:
    using iterator_category = std::input_iterator_tag;
    using value_type = Data<ContainerT>;
    using difference_type = std::ptrdiff_t;
    using pointer = value_type*;
    using reference = value_type&;

    constexpr Iterator(iterator_t<ContainerT>&& sub_iter, std::size_t start)
	: sub_iter_{std::move(sub_iter)}, index_{start} {}

    Data<ContainerT> operator*() { return {index_, *sub_iter_}; }

    constexpr auto operator->() -> ArrowProxy<decltype(**this)> {
      return {**this};
    }

    Iterator& operator++() {
      ++sub_iter_;
      ++index_;
      return *this;
    }

    Iterator operator++(int) {
      auto ret = *this;
      ++*this;
      return ret;
    }

    template <typename T>
    constexpr bool operator!=(const Iterator<T>& other) const {
      return sub_iter_ != other.sub_iter_;
    }

    template <typename T>
    constexpr bool operator==(const Iterator<T>& other) const {
      return !(*this != other);
    }
  };

  constexpr Iterator<Container> begin() {
    return {fancy_getters::begin(container_), start_};
  }
  constexpr Iterator<Container> end() {
    return {fancy_getters::end(container_), start_};
  }

  constexpr Iterator<make_const_t<Container>> begin() const {
    return {fancy_getters::begin(std::as_const(container_)), start_};
  }
  constexpr Iterator<make_const_t<Container>> end() const {
    return {fancy_getters::end(std::as_const(container_)), start_};
  }
};

#endif	//__ITERATORTOOLS_ENUMERATOR_HPP__
