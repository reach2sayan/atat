#ifndef __ITERTOOLS_ITERATOR_WRAPPER_HPP__
#define __ITERTOOLS_ITERATOR_WRAPPER_HPP__

#include <cassert>
#include <functional>
#include <variant>

#include "itertoolsbase.h"

namespace itertools {

// iterator_end_type<C> is the type of C's end iterator
template <typename Container>
using iterator_end_type = decltype(getters::end(std::declval<Container&>()));

template <typename SubIter, typename SubEnd>
class IteratorWrapperImpl;

// If begin and end return the same type, type will be
// iterator_type<Container>
// If begin and end return different types, type will be IteratorWrapperImpl
template <typename Container, bool same_types>
struct IteratorWrapperImplType;

template <typename Container>
struct IteratorWrapperImplType<Container, true>
    : type_is<iterator_type<Container>> {};

template <typename Container>
struct IteratorWrapperImplType<Container, false>
    : type_is<IteratorWrapperImpl<iterator_type<Container>,
				  iterator_end_type<Container>>> {};

template <typename Container>
using IteratorWrapper = typename IteratorWrapperImplType<
    Container, std::is_same_v<iterator_type<Container>,
			      iterator_end_type<Container>>>::type;
}  // namespace itertools

template <typename SubIter, typename SubEnd>
class itertools::IteratorWrapperImpl {
 private:
  static_assert(!std::is_same_v<SubIter, SubEnd>);

  SubIter& sub_iter() {
    auto* sub = std::get_if<SubIter>(&sub_iter_or_end);
    assert(sub);
    return *sub;
  }
  const SubIter& sub_iter() const {
    auto* sub = std::get_if<SubIter>(&sub_iter_or_end);
    assert(sub);
    return *sub;
  }

  std::variant<SubIter, SubEnd> sub_iter_or_end;

 public:
  IteratorWrapperImpl() : IteratorWrapperImpl(SubIter{}) {}

  IteratorWrapperImpl(SubIter&& it) : sub_iter_or_end{std::move(it)} {}

  IteratorWrapperImpl(SubEnd&& it) : sub_iter_or_end(std::move(it)) {}
  IteratorWrapperImpl& operator++() {
    ++sub_iter();
    return *this;
  }

  decltype(auto) operator*() { return *sub_iter(); }

  decltype(auto) operator*() const { return *sub_iter(); }

  decltype(auto) operator->() { return apply_arrow(sub_iter()); }

  decltype(auto) operator->() const { return apply_arrow(sub_iter()); }

  bool operator!=(const IteratorWrapperImpl& other) const {
    constexpr static struct : std::not_equal_to<void> {
      // specially compare Ends because rangev3 sentinels are not equality
      // comparable
      bool operator()(const SubEnd&, const SubEnd&) const { return false; }
      using std::not_equal_to<void>::operator();
    } not_equal;
    return std::visit(not_equal, sub_iter_or_end, other.sub_iter_or_end);
  }

  bool operator==(const IteratorWrapperImpl& other) const {
    return !(*this != other);
  }
};

#endif	// __ITERTOOLS_ITERATOR_WRAPPER_HPP__
