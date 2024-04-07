#ifndef __ITERATORTOOLS_ZIPPER_HPP__
#define __ITERATORTOOLS_ZIPPER_HPP__

// This algorithm utilises std::index_sequence to extract out each input
// container (interpreted as a collection of N-tuple set), the ith container and
// its associated iterators are obtained through the std::get<index>(tuple)

#include "base.h"

namespace ATATIteratorTools {
template <typename TupleType, std::size_t... Is>
class Zipped;

template <typename TupleType, std::size_t... Is>
Zipped<TupleType, Is...> zip_impl(TupleType&&, std::index_sequence<Is...>);

template <typename TupleType, std::size_t... Is>
Zipped<TupleType, Is...> zip_impl(TupleType&& containers,
				  std::index_sequence<Is...>) {
  return {std::move(containers)};
}

template <typename... Containers>
constexpr auto zip(Containers&&... containers) {
  return zip_impl(
      std::tuple<Containers...>{std::forward<Containers>(containers)...},
      std::make_index_sequence<sizeof...(Containers)>{});
}
template <typename... Containers>
constexpr auto zip(Containers&&... containers);
}  // namespace ATATIteratorTools

template <typename TupleType, std::size_t... Is>
class ATATIteratorTools::Zipped {
 private:
  TupleType containers_;

  friend Zipped ATATIteratorTools::zip_impl<TupleType, Is...>(
      TupleType&&, std::index_sequence<Is...>);
  Zipped(TupleType&& containers) : containers_(std::move(containers)) {}

 public:
  Zipped(Zipped&&) noexcept = default;

  // template templates here because I need to defer evaluation in the const
  // iteration case for types that don't have non-const begin() and end(). If I
  // passed in the actual types of the tuples of iterators and the type for
  // deref they'd need to be known in the function declarations below.

  template <typename TupleT, template <typename> class TupleIteratorT,
	    template <typename> class TupleDerefT>
  class Iterator {
   public:
    TupleIteratorT<TupleT> iters_;

   public:
    using iterator_category = std::input_iterator_tag;
    using value_type = TupleDerefT<TupleT>;
    using difference_type = std::ptrdiff_t;
    using pointer = value_type*;
    using reference = value_type;

    constexpr Iterator(TupleIteratorT<TupleT>&& iters)
	: iters_(std::move(iters)) {}

    Iterator& operator++() {
      // the blackhole function essentially ignores the return value
      // of the operator++(), but it regardless has to expand the parameters
      // nonetheless, this executes operator++() for the individual iterators
      blackhole(++std::get<Is>(iters_)...);
      return *this;
    }

    Iterator operator++(int) {
      auto ret = *this;
      ++*this;
      return ret;
    }

    // funny - usually it's the other way around
    template <typename T, template <typename> class IT,
	      template <typename> class TD>
    constexpr bool operator!=(const Iterator<T, IT, TD>& other) const {
      if constexpr (sizeof...(Is) == 0) {
	return false;  // empty is equal, hence this construction
      } else {
	return (... && (std::get<Is>(iters_) != std::get<Is>(other.iters_)));
      }
    }

    template <typename T, template <typename> class IT,
	      template <typename> class TD>
    constexpr bool operator==(const Iterator<T, IT, TD>& other) const {
      return !(*this != other);
    }

    TupleDerefT<TupleT> operator*() { return {(*std::get<Is>(iters_))...}; }

    constexpr auto operator->() -> ArrowProxy<decltype(**this)> {
      return {**this};
    }
  };

  constexpr Iterator<TupleType, iterator_tuple_t, iterator_deref_tuple_t>
  begin() {
    return {{fancy_getters::begin(std::get<Is>(containers_))...}};
  }

  constexpr Iterator<TupleType, iterator_tuple_t, iterator_deref_tuple_t>
  end() {
    return {{fancy_getters::end(std::get<Is>(containers_))...}};
  }

  constexpr Iterator<make_const_t<TupleType>, const_iterator_tuple_t,
		     const_iterator_deref_tuple_t>
  begin() const {
    return {
	{fancy_getters::begin(std::as_const(std::get<Is>(containers_)))...}};
  }

  constexpr Iterator<make_const_t<TupleType>, const_iterator_tuple_t,
		     const_iterator_deref_tuple_t>
  end() const {
    return {{fancy_getters::end(std::as_const(std::get<Is>(containers_)))...}};
  }
};

#endif	//__ITERATORTOOLS_ZIPPER_HPP__
