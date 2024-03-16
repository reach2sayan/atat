#ifndef __ITERATORTOOLS_FILTER_HPP__
#define __ITERATORTOOLS_FILTER_HPP__

#include "base.h"

namespace ATATIteratorTools {

template <typename Predicate, typename Container>
class Filtered;

template <typename Predicate, typename Container>
class FilteredFalse;

struct Boolean {
  template <typename T>
  constexpr bool operator()(const T& item) const {
    return bool(item);
  }
};

template <typename Container, typename Predicate>
constexpr auto filter = [](Container&& container, Predicate&& predFn,
			   bool filter_false = false) {
  if (filter_false == false)
    return Filtered<Container, Predicate>(std::forward<Container>(container),
					  std::forward<Predicate>(predFn));
  else
    return FilteredFalse<Container, Predicate>(
	std::forward<Container>(container), std::forward<Predicate>(predFn));
};

}  // namespace ATATIteratorTools

template <typename Container, typename Predicate>
class ATATIteratorTools::Filtered {
 private:
  Container _container;
  mutable Predicate predFn;

 protected:
  // private Value constructor
  Filtered(Container&& container, Predicate Fn)
      : _container(std::forward<Container>(container)), predFn(Fn) {}

 public:
  Filtered(Filtered&&) = default;

  // this template type includes the cv-qualification of type Container
  template <typename ContainerT>
  class Iterator {
   private:
    template <typename>
    friend class Iterator;

    using DataHolder = DerefDataHolder<iterator_deref<ContainerT>>;
    mutable iterator_t<ContainerT> sub_iter_;
    iterator_t<ContainerT> sub_end_;

    mutable DataHolder item;
    Predicate* predFn_;

    void increment() const {
      ++sub_iter_;
      if (sub_iter_ != sub_end_) {
	item.init(*sub_iter_);
      }
    }

    // increment until the iterator points to is true on the
    // predicate.  Called by constructor and operator++
    void next() const {
      while (sub_iter_ != sub_end_ && !std::invoke(*predFn_, *item)) {
	increment();
      }
    }

    void init() const {
      if (!item && sub_iter_ != sub_end_) {
	item.reset(*sub_iter_);
	next();
      }
    }

   public:
    using iterator_category = std::input_iterator_tag;
    using value_type = iterator_traits_deref<ContainerT>;
    using difference_type = std::ptrdiff_t;
    using pointer = value_type*;
    using reference = value_type&;

    Iterator(iterator_t<ContainerT>&& sub_iter,
	     iterator_t<ContainerT>&& sub_end, Predicate& Fn)
	: sub_iter_{std::move(sub_iter)},
	  sub_end_{std::move(sub_end)},
	  predFn_(&Fn) {}

    typename DataHolder::reference operator*() {
      init();
      return item.operator*();
    }

    typename DataHolder::pointer operator->() {
      init();
      return item.operator->();
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
      return sub_iter_ != other.sub_iter_;
    }

    template <typename T>
    bool operator==(const Iterator<T>& other) const {
      return !(*this != other);
    }
  };

  Iterator<Container> begin() {
    return {fancy_getters::begin(_container), fancy_getters::end(_container),
	    predFn};
  }

  Iterator<Container> end() {
    return {fancy_getters::end(_container), fancy_getters::end(_container),
	    predFn};
  }

  Iterator<make_const_t<Container>> begin() const {
    return {fancy_getters::begin(std::as_const(_container)),
	    fancy_getters::end(std::as_const(_container)), predFn};
  }

  Iterator<Container> end() const {
    return {fancy_getters::end(std::as_const(_container)),
	    fancy_getters::end(std::as_const(_container)), predFn};
  }
};

template <typename Container, typename Predicate>
class ATATIteratorTools::FilteredFalse
    : public Filtered<Container, Not<Predicate>> {
  FilteredFalse(Container&& container, Predicate predFn)
      : Filtered<Container, Not<Predicate>>(std::forward<Container>(container),
					    {predFn}) {}
};

#endif	//__ITERATORTOOLS_FILTER_HPP__
