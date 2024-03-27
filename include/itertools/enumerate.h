#ifndef __ITERATORTOOLS_ENUMERATOR_HPP__
#define __ITERATORTOOLS_ENUMERATOR_HPP__

#include "base.h"

namespace ATATIteratorTools {
// Held by the Enumerable::Iterator.  Has a .index, and a
// .element referencing the value yielded by the subiterator
template <typename Index, typename Elem>
class EnumeratorDataHolder : public std::pair<Index, Elem> {
  using std::pair<Index, Elem>::pair;

 public:
  typename std::pair<Index, Elem>::first_type index =
      std::pair<Index, Elem>::first;
  typename std::pair<Index, Elem>::second_type element =
      std::pair<Index, Elem>::second;
};

template <typename Container, typename Index>
class Enumerable;

using EnumerateObject = EnumeratorClosureObject<Enumerable, std::size_t>;
constexpr EnumerateObject enumerate{};

}  // namespace ATATIteratorTools

// Partial specialization of std::tuple related templates for data-holder
namespace std {
template <typename Index, typename Elem>
struct tuple_size<ATATIteratorTools::EnumeratorDataHolder<Index, Elem>>
    : public tuple_size<std::pair<Index, Elem>> {};

template <std::size_t N, typename Index, typename Elem>
struct tuple_element<N, ATATIteratorTools::EnumeratorDataHolder<Index, Elem>>
    : public tuple_element<N, std::pair<Index, Elem>> {};
}  // namespace std

template <typename Container, typename Index>
class ATATIteratorTools::Enumerable {
 private:
  Container container_;
  const Index start_;

  friend EnumeratorClosureObject<Enumerable, std::decay_t<Index>>;

 protected:
  // private Value constructor
  Enumerable(Container&& container, Index start)
      : container_(std::forward<Container>(container)), start_{start} {}

 public:
  Enumerable(Enumerable&&) = default;

  template <typename T>
  using Data = EnumeratorDataHolder<Index, iterator_deref<T>>;

  //  Holds an iterator of the contained type and an Index for the
  // index_.  Each call to ++ increments both of these data members.
  // Each dereference returns an IterYield.
  template <typename ContainerT>
  class Iterator {
   private:
    template <typename>
    friend class Iterator;

    iterator_t<ContainerT> sub_iter_;
    Index index_;

   public:
    using iterator_category = std::input_iterator_tag;
    using value_type = Data<ContainerT>;
    using difference_type = std::ptrdiff_t;
    using pointer = value_type*;
    using reference = value_type&;

    Iterator(iterator_t<ContainerT>&& sub_iter, Index start)
	: sub_iter_{std::move(sub_iter)}, index_{start} {}

    Data<ContainerT> operator*() { return {index_, *sub_iter_}; }

    ArrowProxy<Data<ContainerT>> operator->() { return {**this}; }

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
    bool operator!=(const Iterator<T>& other) const {
      return sub_iter_ != other.sub_iter_;
    }

    template <typename T>
    bool operator==(const Iterator<T>& other) const {
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
