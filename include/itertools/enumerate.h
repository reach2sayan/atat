#ifndef __ITERTOOLS_ENUMERATE_HPP__
#define __ITERTOOLS_ENUMERATE_HPP__

#include "itertoolsbase.h"
#include "itertoolswrapper.h"

namespace itertools {
template <typename Index, typename Elem>
using EnumBasePair = std::pair<Index, Elem>;

// "yielded" by the Enumerable::Iterator.  Has a .index, and a
// .element referencing the value yielded by the subiterator
template <typename Index, typename Elem>
class EnumIteratorYield : public EnumBasePair<Index, Elem> {
  using BasePair = EnumBasePair<Index, Elem>;
  using BasePair::BasePair;

 public:
  typename BasePair::first_type index = BasePair::first;
  typename BasePair::second_type element = BasePair::second;
};

template <typename Container, typename Index>
class Enumerable;

using EnumerateFn = IterToolFnOptionalBindSecond<Enumerable, std::size_t>;
constexpr EnumerateFn enumerate{};
}  // namespace itertools

namespace std {
template <typename Index, typename Elem>
struct tuple_size<itertools::EnumIteratorYield<Index, Elem>>
    : public tuple_size<itertools::EnumBasePair<Index, Elem>> {};

template <std::size_t N, typename Index, typename Elem>
struct tuple_element<N, itertools::EnumIteratorYield<Index, Elem>>
    : public tuple_element<N, itertools::EnumBasePair<Index, Elem>> {};
}  // namespace std

template <typename Container, typename Index>
class itertools::Enumerable {
 private:
  Container _container;
  const Index _start;
  friend EnumerateFn;
  // Value constructor for use only in the enumerate function
  Enumerable(Container&& container, Index start)
      : _container(std::forward<Container>(container)), _start{start} {}

 public:
  Enumerable(Enumerable&&) = default;

  template <typename T>
  using IteratorYield = EnumIteratorYield<Index, iterator_deref<T>>;

  //  Holds an iterator of the contained type and an Index for the
  //  index_.  Each call to ++ increments both of these data members.
  //  Each dereference returns an IterYield.
  template <typename ContainerT>
  class Iterator {
   private:
    template <typename>
    friend class Iterator;
    IteratorWrapper<ContainerT> _sub_iter;
    Index _index;

   public:
    using iterator_category = std::input_iterator_tag;
    using value_type = IteratorYield<ContainerT>;
    using difference_type = std::ptrdiff_t;
    using pointer = value_type*;
    using reference = value_type&;

    Iterator(IteratorWrapper<ContainerT>&& sub_iter, Index start)
	: _sub_iter{std::move(sub_iter)}, _index{start} {}

    IteratorYield<ContainerT> operator*() { return {_index, *_sub_iter}; }

    ArrowProxy<IteratorYield<ContainerT>> operator->() { return {**this}; }

    Iterator& operator++() {
      ++_sub_iter;
      ++_index;
      return *this;
    }

    Iterator operator++(int) {
      auto ret = *this;
      ++*this;
      return ret;
    }

    template <typename T>
    bool operator!=(const Iterator<T>& other) const {
      return _sub_iter != other._sub_iter;
    }

    template <typename T>
    bool operator==(const Iterator<T>& other) const {
      return !(*this != other);
    }
  };

  Iterator<Container> begin() { return {getters::begin(_container), _start}; }

  Iterator<Container> end() { return {getters::end(_container), _start}; }

  Iterator<AsConst<Container>> begin() const {
    return {getters::begin(std::as_const(_container)), _start};
  }

  Iterator<AsConst<Container>> end() const {
    return {getters::end(std::as_const(_container)), _start};
  }
};
#endif	// __ITERTOOLS_ENUMERATE_HPP__
