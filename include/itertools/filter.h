#ifndef __ITERATORTOOLS_FILTER_HPP__
#define __ITERATORTOOLS_FILTER_HPP__

#include "base.h"

namespace ATATIteratorTools {

template <typename Predicate, typename Container>
class Filtered;

struct Boolean {
  template <typename T>
  constexpr bool operator()(const T& reference_item) const {
    return bool(reference_item);
  }
};

using FilterTrueObject = FilterClosureObject<Filtered, Boolean>;
constexpr FilterTrueObject filter{};

}  // namespace ATATIteratorTools

template <typename Predicate, typename Container>
class ATATIteratorTools::Filtered {
 private:
  Container _container;
  mutable Predicate predFn;
  const bool _use_false;

  friend FilterTrueObject;

 protected:
  // private Value constructor
  Filtered(Predicate Fn, Container&& container, bool use_false)
      : _container(std::forward<Container>(container)),
	predFn{Fn},
	_use_false{use_false} {}

 public:
  Filtered(Filtered&&) = default;

  // this template type includes the cv-qualification of type Container
  template <typename ContainerT>
  class Iterator {
   private:
    template <typename>
    friend class Iterator;

    using DataHolder = DereferencedDataHolder<iterator_deref<ContainerT>>;
    mutable iterator_t<ContainerT> input_iterator;
    iterator_t<ContainerT> input_iterator_end;

    mutable DataHolder reference_item;
    Predicate* predFn;
    bool iterator_use_false;

    void increment() const {
      ++input_iterator;
      if (input_iterator != input_iterator_end) {
	reference_item.set_data(*input_iterator);
      }
    }

    // increment until the iterator points to is true on the
    // predicate.  Called by constructor and operator++
    void next() const {
      while ((input_iterator != input_iterator_end) &&
	     (iterator_use_false != !std::invoke(*predFn, *reference_item))) {
	increment();
      }
    }

    void init() const {
      if (!reference_item && input_iterator != input_iterator_end) {
	reference_item.set_data(*input_iterator);
	next();
      }
    }

   public:
    using iterator_category = std::input_iterator_tag;
    using value_type = iterator_traits_deref<ContainerT>;
    using difference_type = std::ptrdiff_t;
    using pointer = value_type*;
    using reference = value_type&;

    Iterator(iterator_t<ContainerT>&& _input_iterator,
	     iterator_t<ContainerT>&& _input_iterator_end, Predicate& Fn,
	     bool use_false = false)
	: input_iterator{std::move(_input_iterator)},
	  input_iterator_end{std::move(_input_iterator_end)},
	  predFn(&Fn),
	  iterator_use_false(use_false) {}

    typename DataHolder::reference operator*() {
      init();
      return reference_item.operator*();
    }

    typename DataHolder::pointer operator->() {
      init();
      return reference_item.operator->();
    }
    Iterator& operator++() {
      init();
      increment();
      next();
      return *this;
    }

    Iterator operator++(int) {
      auto ret = *this;
      ++*this;
      return ret;
    }

    template <typename T>
    bool operator!=(const Iterator<T>& other) const {
      init();
      other.init();
      return input_iterator != other.input_iterator;
    }

    template <typename T>
    bool operator==(const Iterator<T>& other) const {
      return !(*this != other);
    }
  };

  Iterator<Container> begin() {
    return {fancy_getters::begin(_container), fancy_getters::end(_container),
	    predFn, _use_false};
  }

  Iterator<Container> end() {
    return {fancy_getters::end(_container), fancy_getters::end(_container),
	    predFn, _use_false};
  }

  Iterator<make_const_t<Container>> begin() const {
    return {fancy_getters::begin(std::as_const(_container)),
	    fancy_getters::end(std::as_const(_container)), predFn, _use_false};
  }

  Iterator<Container> end() const {
    return {fancy_getters::end(std::as_const(_container)),
	    fancy_getters::end(std::as_const(_container)), predFn, _use_false};
  }
};

#endif	//__ITERATORTOOLS_FILTER_HPP__
