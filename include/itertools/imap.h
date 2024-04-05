#ifndef __ITERATORTOOLS_IMAP_HPP__
#define __ITERATORTOOLS_IMAP_HPP__

#include <array>
#include <iterator>

#include "base.h"

namespace ATATIteratorTools {
template <typename Func, typename Container>
class IMapper;

struct IMapClosureObject {
  template <typename MapFunc, typename Container>
  auto operator()(MapFunc&& map_func, Container&& container) const {
    return IMapper<MapFunc, Container>{std::forward<MapFunc&&>(map_func),
				       std::forward<Container>(container)};
  }
};
constexpr IMapClosureObject imap{};

}  // namespace ATATIteratorTools

template <typename Func, typename Container>
class ATATIteratorTools::IMapper {
 private:
  mutable Func func_;
  Container container_;

  using IterDeref =
      decltype(std::apply(func_, std::declval<iterator_deref<Container>>()));
  using IterDerefValue = std::remove_cv_t<std::remove_reference_t<IterDeref>>;

  template <typename FFunc, typename CContainer>
  constexpr IMapper(FFunc&& f, CContainer&& c)
      : func_(std::forward<FFunc>(f)),
	container_(std::forward<CContainer>(c)) {}

  friend IMapClosureObject;

 public:
  template <typename ContainerT>
  class Iterator {
   private:
    template <typename>
    friend class Iterator;
    Func* func_;
    iterator_t<ContainerT> sub_iter_;

   public:
    using iterator_category = std::input_iterator_tag;
    using value_type = IterDerefValue;
    using difference_type = std::ptrdiff_t;
    using pointer = value_type*;
    using reference = IterDeref;

    constexpr Iterator(Func& f, iterator_t<ContainerT>&& sub_iter)
	: func_(&f), sub_iter_(std::move(sub_iter)) {}

    template <typename T>
    constexpr bool operator!=(const Iterator<T>& other) const {
      return sub_iter_ != other.sub_iter_;
    }

    template <typename T>
    constexpr bool operator==(const Iterator<T>& other) const {
      return !(*this != other);
    }

    Iterator& operator++() {
      ++sub_iter_;
      return *this;
    }

    Iterator operator++(int) {
      auto ret = *this;
      ++*this;
      return ret;
    }

    decltype(auto) operator*() { return std::apply(*func_, *sub_iter_); }

    auto operator->() -> ArrowProxy<decltype(**this)> { return {**this}; }
  };

  constexpr Iterator<Container> begin() {
    return {func_, get_begin(container_)};
  }

  constexpr Iterator<Container> end() { return {func_, get_end(container_)}; }

  constexpr Iterator<make_const_t<Container>> begin() const {
    return {func_, fancy_getters::begin(std::as_const(container_))};
  }

  constexpr Iterator<make_const_t<Container>> end() const {
    return {func_, fancy_getters::end(std::as_const(container_))};
  }
};

#endif	// __ITERATORTOOLS_IMAP_HPP__
