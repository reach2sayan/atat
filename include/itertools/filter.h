#ifndef __ITERATORTOOLS_FILTER_HPP__
#define __ITERATORTOOLS_FILTER_HPP__

#include "base.h"

namespace ATATIteratorTools {

// this is default predicate when no predicate is supplied, this essentially
// runs the bool() operator on the dataype,
struct Boolean {
  template <typename T>
  constexpr bool operator()(const T& reference_item) const {
    return bool(reference_item);
  }
};

template <typename Predicate, typename Container>
class Filtered;

struct FilterClosureObject {
 public:
  // when only the iterable is passed, we perform (!)bool(object)
  template <typename Container,
	    typename = std::enable_if_t<is_iterable_v<Container>>>
  constexpr auto operator()(Container&& container,
			    const bool use_false = false) const {
    return (*this)(Boolean{}, std::forward<Container>(container), use_false);
  }

  // this is when a explicit predicate is passed
  template <typename Predicate, typename Container,
	    typename = std::enable_if_t<is_iterable_v<Container>>>
  constexpr Filtered<Predicate, Container> operator()(
      Predicate predFn, Container&& container,
      const bool use_false = false) const {
    return {std::move(predFn), std::forward<Container>(container), use_false};
  }
};

constexpr FilterClosureObject filter{};

}  // namespace ATATIteratorTools

template <typename Predicate, typename Container>
class ATATIteratorTools::Filtered {
 private:
  Container _container;
  mutable Predicate predFn;
  const bool _use_false;

  friend FilterClosureObject;

 protected:
  // private Value constructor
  constexpr Filtered(Predicate Fn, Container&& container, bool use_false)
      : _container(std::forward<Container>(container)),
	predFn{Fn},
	_use_false{use_false} {}

 public:
  constexpr Filtered(Filtered&&) = default;

  // this template type includes the cv-qualification of type Container
  template <typename ContainerT>
  class Iterator {
   private:
    template <typename>
    friend class Iterator;

    using DataHolder = DereferencedDataHolder<iterator_deref_t<ContainerT>>;
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
	     // the code below is essentially a XOR gate, that flips the sign of
	     // the second bit if the first bit is true
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
    using value_type = iterator_traits_deref_t<ContainerT>;
    using difference_type = std::ptrdiff_t;
    using pointer = value_type*;
    using reference = value_type&;

    constexpr Iterator(iterator_t<ContainerT>&& _input_iterator,
		       iterator_t<ContainerT>&& _input_iterator_end,
		       Predicate& Fn, bool use_false = false)
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
    constexpr bool operator==(const Iterator<T>& other) const {
      return !(*this != other);
    }
  };

  constexpr Iterator<Container> begin() {
    return {fancy_getters::begin(_container), fancy_getters::end(_container),
	    predFn, _use_false};
  }

  constexpr Iterator<Container> end() {
    return {fancy_getters::end(_container), fancy_getters::end(_container),
	    predFn, _use_false};
  }

  constexpr Iterator<make_const_t<Container>> begin() const {
    return {fancy_getters::begin(std::as_const(_container)),
	    fancy_getters::end(std::as_const(_container)), predFn, _use_false};
  }

  constexpr Iterator<Container> end() const {
    return {fancy_getters::end(std::as_const(_container)),
	    fancy_getters::end(std::as_const(_container)), predFn, _use_false};
  }
};

#endif	//__ITERATORTOOLS_FILTER_HPP__
