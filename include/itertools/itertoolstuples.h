#ifndef __ITERTOOLS_TUPLES_HPP__
#define __ITERTOOLS_TUPLES_HPP__

#include "itertoolsbase.h"
#include "itertoolswrapper.h"

namespace itertools {

template <typename... Ts>
std::tuple<iterator_deref<Ts>...> iterator_tuple_deref_helper(
    const std::tuple<Ts...>&);

template <typename... Ts>
std::tuple<IteratorWrapper<Ts>...> iterator_tuple_type_helper(
    const std::tuple<Ts...>&);

template <typename... Ts>
std::tuple<iterator_deref<const std::remove_reference_t<Ts>>...>
const_iterator_tuple_deref_helper(const std::tuple<Ts...>&);

template <typename... Ts>
std::tuple<IteratorWrapper<const std::remove_reference_t<Ts>>...>
const_iterator_tuple_type_helper(const std::tuple<Ts...>&);

// Given a tuple template argument, evaluates to a tuple of iterators
// for the template argument's contained types.
template <typename TupleType>
using iterator_tuple_type =
    decltype(iterator_tuple_type_helper(std::declval<TupleType>()));

template <typename TupleType>
using const_iterator_tuple_type =
    decltype(const_iterator_tuple_type_helper(std::declval<TupleType>()));

// Given a tuple template argument, evaluates to a tuple of
// what the iterators for the template argument's contained types
// dereference to
template <typename TupleType>
using iterator_deref_tuple =
    decltype(iterator_tuple_deref_helper(std::declval<TupleType>()));

template <typename TupleType>
using const_iterator_deref_tuple =
    decltype(const_iterator_tuple_deref_helper(std::declval<TupleType>()));

// function absorbing all arguments passed to it. used when
// applying a function to a parameter pack but not passing the evaluated
// results anywhere
template <typename... Ts>
void absorb(Ts&&...) {}

}  // namespace itertools

#endif	// __ITERTOOLS_TUPLES_HPP__
