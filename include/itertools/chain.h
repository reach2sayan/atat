#ifndef __ITERATORTOOLS_CHAIN_HPP__
#define __ITERATORTOOLS_CHAIN_HPP__

#include "base.h"

namespace ATATIteratorTools {

template <typename TupType, std::size_t... Is>
class Chained;

template <typename>
struct TupleOfConstImpl;

template <typename... Ts>
struct TupleOfConstImpl<std::tuple<Ts...>>
    : type_is<std::tuple<make_const_t<Ts>...>> {};

template <typename T>
using make_tuple_const_t = typename TupleOfConstImpl<T>::type;

class ChainerClosureObject {
 private:
  template <typename TupleType, std::size_t... Is>
  Chained<TupleType, Is...> chain_impl(TupleType&& containers,
				       std::index_sequence<Is...>) const {
    return {std::move(containers)};
  }

 public:
  template <typename... Containers>
  auto operator()(Containers&&... cs) const {
    return chain_impl(
	std::tuple<Containers...>{std::forward<Containers>(cs)...},
	std::index_sequence_for<Containers...>{});
  }
};

constexpr auto chain = ChainerClosureObject{};

}  // namespace ATATIteratorTools

template <typename TupleType, std::size_t... Is>
class ATATIteratorTools::Chained {
 private:
  friend ChainerClosureObject;

  template <typename TupleT>
  class IteratorData {
    static_assert(std::tuple_size<std::decay_t<TupleT>>::value == sizeof...(Is),
		  "tuple size != sizeof Is");

    static_assert(
	are_same<iterator_deref_t<std::tuple_element_t<Is, TupleT>>...>::value,
	"All chained iterables must have iterators that "
	"dereference to the same type, including cv-qualifiers "
	"and references.");

   public:
    using TupleIteratorT = iterator_tuple_t<TupleT>;
    using DerefT = iterator_deref_t<std::tuple_element_t<0, TupleT>>;
    using ArrowT = iterator_arrow_t<std::tuple_element_t<0, TupleT>>;

    IteratorData() = delete;

    template <std::size_t Idx>
    static DerefT get_and_deref(TupleIteratorT& iters) {
      return *std::get<Idx>(iters);
    }

    template <std::size_t Idx>
    static ArrowT get_and_arrow(TupleIteratorT& iters) {
      return apply_arrow(std::get<Idx>(iters));
    }

    template <std::size_t Idx>
    static void get_and_increment(TupleIteratorT& iters) {
      ++std::get<Idx>(iters);
    }

    template <std::size_t Idx>
    static bool get_and_check_not_equal(const TupleIteratorT& lhs,
					const TupleIteratorT& rhs) {
      return std::get<Idx>(lhs) != std::get<Idx>(rhs);
    }

    using DerefFunc = DerefT (*)(TupleIteratorT&);
    using ArrowFunc = ArrowT (*)(TupleIteratorT&);
    using IncFunc = void (*)(TupleIteratorT&);
    using NeqFunc = bool (*)(const TupleIteratorT&, const TupleIteratorT&);

    constexpr static std::array<DerefFunc, sizeof...(Is)> derefers{
	{get_and_deref<Is>...}};

    constexpr static std::array<ArrowFunc, sizeof...(Is)> arrowers{
	{get_and_arrow<Is>...}};

    constexpr static std::array<IncFunc, sizeof...(Is)> incrementers{
	{get_and_increment<Is>...}};

    constexpr static std::array<NeqFunc, sizeof...(Is)> neq_comparers{
	{get_and_check_not_equal<Is>...}};

    using tuple_underlying_t = std::remove_reference<
	iterator_deref_t<std::tuple_element_t<0, TupleT>>>;
  };

  Chained(TupleType&& t) : tup_(std::move(t)) {}
  TupleType tup_;

 public:
  Chained(Chained&&) noexcept = default;

  template <typename TupleT>
  class Iterator {
   private:
    std::size_t index_;
    typename IteratorData<TupleT>::TupleIteratorT iters_;
    typename IteratorData<TupleT>::TupleIteratorT ends_;

    void check_for_end_and_adjust() {
      while (index_ < sizeof...(Is) &&
	     !(IteratorData<TupleT>::neq_comparers[index_](iters_, ends_))) {
	++index_;
      }
    }

   public:
    using iterator_category = std::input_iterator_tag;
    using value_type = typename IteratorData<TupleT>::tuple_underlying_t;
    using difference_type = std::ptrdiff_t;
    using pointer = value_type*;
    using reference = value_type&;

    Iterator(std::size_t i,
	     typename IteratorData<TupleT>::TupleIteratorT&& iters,
	     typename IteratorData<TupleT>::TupleIteratorT&& ends)
	: index_{i}, iters_(std::move(iters)), ends_(std::move(ends)) {
      check_for_end_and_adjust();
    }

    decltype(auto) operator*() {
      return IteratorData<TupleT>::derefers[index_](iters_);
    }

    decltype(auto) operator->() {
      return IteratorData<TupleT>::arrowers[index_](iters_);
    }

    Iterator& operator++() {
      IteratorData<TupleT>::incrementers[index_](iters_);
      check_for_end_and_adjust();
      return *this;
    }

    Iterator operator++(int) {
      auto ret = *this;
      ++*this;
      return ret;
    }

    // TODO make const and non-const iterators comparable
    bool operator!=(const Iterator& other) const {
      return index_ != other.index_ ||
	     (index_ != sizeof...(Is) &&
	      IteratorData<TupleT>::neq_comparers[index_](iters_,
							  other.iters_));
    }

    bool operator==(const Iterator& other) const { return !(*this != other); }
  };

  Iterator<TupleType> begin() {
    return {0,
	    {fancy_getters::begin(std::get<Is>(tup_))...},
	    {fancy_getters::end(std::get<Is>(tup_))...}};
  }

  Iterator<TupleType> end() {
    return {sizeof...(Is),
	    {fancy_getters::end(std::get<Is>(tup_))...},
	    {fancy_getters::end(std::get<Is>(tup_))...}};
  }

  Iterator<make_tuple_const_t<TupleType>> begin() const {
    return {0,
	    {fancy_getters::begin(std::as_const(std::get<Is>(tup_)))...},
	    {fancy_getters::end(std::as_const(std::get<Is>(tup_)))...}};
  }

  Iterator<make_tuple_const_t<TupleType>> end() const {
    return {sizeof...(Is),
	    {fancy_getters::end(std::as_const(std::get<Is>(tup_)))...},
	    {fancy_getters::end(std::as_const(std::get<Is>(tup_)))...}};
  }
};

#endif	//__ITERATORTOOLS_CHAIN_HPP__
