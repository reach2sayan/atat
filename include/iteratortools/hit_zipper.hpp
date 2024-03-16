#ifndef __HYSYS_ITERATORTOOLS_ZIPPER_HPP__
#define __HYSYS_ITERATORTOOLS_ZIPPER_HPP__

#include "hysysiteratortools.hpp"

namespace HysysIteratorTools
	{
    template <typename TupleType, std::size_t... Is>
    class Zipped;

    template <typename TupleType, std::size_t... Is>
    Zipped<TupleType, Is...> zip_impl(TupleType&&, std::index_sequence<Is...>);

    template <typename... Containers>
    decltype(auto) zip(Containers&&... containers);
	}

template <typename TupleType, std::size_t... Is>
class HysysIteratorTools::Zipped
    {
    private:
        TupleType containers_;
        
        friend Zipped HysysIteratorTools::zip_impl<TupleType, Is...>(TupleType&&, std::index_sequence<Is...>);
        Zipped(TupleType&& containers) : containers_(std::move(containers)) {}
    
    public:
        
        Zipped(Zipped&&) noexcept = default;
        
        // template templates here because I need to defer evaluation in the const
        // iteration case for types that don't have non-const begin() and end(). If I
        // passed in the actual types of the tuples of iterators and the type for
        // deref they'd need to be known in the function declarations below.
        
        template <typename TupleTypeT, template <typename> class IteratorTuple, template <typename> class TupleDeref>
        class Iterator 
            {
            public:
                IteratorTuple<TupleTypeT> iters_;
            
            public:
                using iterator_category = std::input_iterator_tag;
                using value_type = TupleDeref<TupleTypeT>;
                using difference_type = std::ptrdiff_t;
                using pointer = value_type*;
                using reference = value_type;

                Iterator(IteratorTuple<TupleTypeT>&& iters) : iters_(std::move(iters)) {}

                Iterator& operator++() 
                    {
                    blackhole(++std::get<Is>(iters_)...);
                    return *this;
                    }
                
                Iterator operator++(int) 
                    {
                    auto ret = *this;
                    ++*this;
                    return ret;
                    }
                
                template <typename T, template <typename> class IT, template <typename> class TD>
                bool operator!=(const Iterator<T, IT, TD>& other) const 
                    {
                    if constexpr (sizeof...(Is) == 0) 
                        {
                        return false;
                        }
                    else 
                        {
                        return (... && (std::get<Is>(iters_) != std::get<Is>(other.iters_)));
                        }
                    }
                
                template <typename T, template <typename> class IT, template <typename> class TD>
                bool operator==(const Iterator<T, IT, TD>& other) const 
                    {
                    return !(*this != other);
                    }
                
                TupleDeref<TupleTypeT> operator*() 
                    {
                    return { (*std::get<Is>(iters_))... };
                    }
                
                auto operator->() -> ArrowProxy<decltype(**this)> 
                    {
                    return { **this };
                    }
            };
        
        Iterator<TupleType, iterator_tuple_t, iterator_deref_tuple> begin() 
            {
            return { {fancy_getters::begin(std::get<Is>(containers_))...} };
            }
        
        Iterator<TupleType, iterator_tuple_t, iterator_deref_tuple> end() 
            {
            return { {fancy_getters::end(std::get<Is>(containers_))...} };
            }
        
        Iterator<make_const_t<TupleType>, const_iterator_tuple_t, const_iterator_deref_tuple> begin() const 
            {
            return { {fancy_getters::begin(std::as_const(std::get<Is>(containers_)))...} };
            }
        
        Iterator<make_const_t<TupleType>, const_iterator_tuple_t, const_iterator_deref_tuple> end() const 
            {
            return { {fancy_getters::end(std::as_const(std::get<Is>(containers_)))...} };
            }
    };

template <typename TupleType, std::size_t... Is>
HysysIteratorTools::Zipped<TupleType, Is...> HysysIteratorTools::zip_impl(TupleType&& containers, std::index_sequence<Is...>) 
    {
    return { std::move(containers) };
    }

template <typename... Containers>
decltype(auto) HysysIteratorTools::zip(Containers&&... containers) 
    {
    return zip_impl(std::tuple<Containers...>{std::forward<Containers>(containers)...}, std::index_sequence_for<Containers...>{});
    }
#endif //__HYSYS_ITERATORTOOLS_ZIPPER_HPP__