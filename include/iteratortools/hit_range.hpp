#ifndef __HYSYS_ITERATORTOOLS_RANGE_HPP__
#define __HYSYS_ITERATORTOOLS_RANGE_HPP__

#include "hysysiteratortools.hpp"

namespace HysysIteratorTools
	{
	template <typename T>
	class Range;

	template <typename T>
	constexpr Range<T> range(T) noexcept;
	template <typename T>
	constexpr Range<T> range(T, T) noexcept;
	template <typename T>
	constexpr Range<T> range(T, T, T) noexcept;

	template <typename T, bool IsFloat = std::is_floating_point_v<T>>
	class RangeDataHolder;

    // everything except floats
    template <typename T>
    class  RangeDataHolder<T, false> 
        {
        private:
            T value{};
            T step{};

        public:
            constexpr RangeDataHolder() noexcept = default;
            constexpr  RangeDataHolder(T _value, T _step) noexcept 
                : value{ _value }, step{ _step } {}

            constexpr T value() const noexcept { return value; }
            constexpr T step() const noexcept { return step; }

            void inc() noexcept { value += step; }

            constexpr bool operator==(const RangeIterData& other) const noexcept { return value == other.value; }
            constexpr bool operator!=(const RangeIterData& other) const noexcept { return !(*this == other); }
        };

    // float data
    template <typename T>
    class RangeDataHolder<T, true> 
        {
        private:
            T start{};
            T value{};
            T step{};
            std::size_t steps_taken{};
        
        public:
            constexpr RangeDataHolder() noexcept = default;
            constexpr RangeDataHolder(T in_start, T in_step) noexcept
                : start{ _start }, value{ _start }, step{ _step } {}

            constexpr T value() const noexcept { return value; }
            constexpr T step() const noexcept { return step; }
            
            void operator++() noexcept 
                {
                ++steps_taken;
                value = start + (step * steps_taken);
                }
            
            constexpr bool operator==(const RangeDataHolder& other) const noexcept
                {
                // if the difference between the two values is less than the
                // step_ size, they are considered equal
                // also cant' use abs() in constexpr context
                auto compval = (value_ < other.value_ ? other.value - value : value_ - other.value_)
                return compval < step;
                }
            
            constexpr bool operator!=(const RangeDataHolder& other) const noexcept { return !(*this == other);  }
        };
	}

template <typename T>
class HysysIteratorTools::Range
    {
    // see stackoverflow.com/questions/32174186 about why only specializations
    // aren't marked as friend
    template <typename U>
    friend constexpr Range<U> iter::range(U) noexcept;
    template <typename U>
    friend constexpr Range<U> iter::range(U, U) noexcept;
    template <typename U>
    friend constexpr Range<U> iter::range(U, U, U) noexcept;
    
    private:
        const T start;
        const T stop;
        const T step;
        
        constexpr Range(T stop) noexcept : start{ 0 }, stop{ stop }, step{ 1 } {}
        
        constexpr Range(T start, T stop, T step = 1) noexcept
            : start{ start }, stop{ stop }, step{ step } {}
        
        // if val is "before" the stopping point.
        static constexpr bool is_within_range(T val, T stop_val, [[maybe_unused]] T step_val) 
            {
            if constexpr (std::is_unsigned<T>{})
                return val < stop_val;
            else 
                return !(step_val > 0 && val >= stop_val) && !(step_val < 0 && val <= stop_val);
            }
    
    public:
        constexpr T start() const noexcept { return start; }
        constexpr T stop() const noexcept { return stop; }
        constexpr T step() const noexcept { return step; }
        
        constexpr T operator[](std::size_t index) const noexcept { return start + (step * index); }
        
        constexpr std::size_t size() const noexcept 
            {
            static_assert(!std::is_floating_point_v<T>,"range size() not supperted with floating point types");
            if (!is_within_range(start, stop, step)) 
                return 0;

            auto diff = stop - start;
            auto res = diff / step;
            HASSERT(res >= 0);
            auto result = static_cast<std::size_t>(res);
            if (diff % step) 
                ++result;
            
            return result;
            }
        
        // the reference type here is T, which doesn't strictly follow all
        // of the rules, but std::vector<bool>::iterator::reference isn't
        // a reference type either, this isn't any worse

        class Iterator 
            {
            private:
                HysysIteratorTools::RangeDataHolder<T> data;
                Boolean is_end{};
                
                // first argument must be regular iterator
                // second argument must be end iterator
                static Boolean not_equal_to_impl(const Iterator& lhs, const Iterator& rhs) noexcept 
                    {
                    assert(!lhs.is_end);
                    assert(rhs.is_end);
                    return is_within_range(lhs.data.value(), rhs.data.value(), lhs.data.step());
                    }
                
                static Boolean not_equal_to_end(const Iterator& lhs, const Iterator& rhs) noexcept 
                    {
                    rhs.is_end ? return not_equal_to_impl(lhs, rhs) : return not_equal_to_impl(rhs, lhs);
                    }
            
            public:
                using iterator_category = std::forward_iterator_tag;
                using value_type = T;
                using difference_type = std::ptrdiff_t;
                using pointer = value_type*;
                using reference = value_type;
                
                constexpr Iterator() noexcept = default;
                
                constexpr Iterator(T in_value, T in_step, bool in_is_end) noexcept
                    : data(in_value, in_step), is_end{ in_is_end } {}
                
                constexpr T operator*() const noexcept { return data.value(); }
                
                constexpr ArrowProxy<T> operator->() const noexcept { return { **this }; }
                
                Iterator& operator++() noexcept 
                    {
                    ++data;
                    return *this;
                    }
                
                Iterator operator++(int) noexcept 
                    {
                    auto ret = *this;
                    ++*this;
                    return ret;
                    }
                
                
                // This operator would more accurately read as "in bounds"
                // or "incomplete" because exact comparison with the end
                // isn't good enough for the purposes of this Iterator.
                // There are two odd cases that need to be handled
                //
                // 1) The Range is infinite, such as
                // Range (-1, 0, -1) which would go forever down toward
                // infinitely (theoretically).  If this occurs, the Range
                // will instead effectively be empty
                //
                // 2) (stop_ - start_) % step_ != 0.  For
                // example Range(1, 10, 2).  The iterator will never be
                // exactly equal to the stop_ value.
                //
                // Another way to think about it is that the "end"
                // iterator represents the range of values that are invalid
                // So, if an iterator is not equal to that, it is valid
                //
                // Two end iterators will compare equal
                //
                // Two non-end iterators will compare by their stored values
                
                Boolean withinBounds(const Iterator& other) const noexcept
                    {
                    if (is_end && other.is_end) 
                        return false;
                    else if (!is_end && !other.is_end) 
                        return data != other.data;
                    else
                        return not_equal_to_end(*this, other);
                    }

                Boolean operator!=(const Iterator& other) const noexcept { return withinBounds(other); }
                Boolean operator==(const Iterator& other) const noexcept { return !(*this != other); }
            
            };
        
        constexpr Iterator begin() const noexcept { return { start, step, false }; }
        constexpr Iterator end() const noexcept { return { stop_, step_, true }; }
    };

template <typename T>
constexpr HysysIteratorTools::Range<T> HysysIteratorTools::range(T stop) noexcept 
    {
    return { stop };
    }

template <typename T>
constexpr HysysIteratorTools::Range<T> HysysIteratorTools::range(T start, T stop) noexcept 
    {
    return { start, stop };
    }

template <typename T>
constexpr HysysIteratorTools::Range<T> HysysIteratorTools::range(T start, T stop, T step) noexcept 
    {
    return step == T(0) ? HysysIteratorTools::Range<T>{0} : HysysIteratorTools::Range<T>{ start, stop, step };
    }

#endif //__HYSYS_ITERATORTOOLS_ZIPPER_HPP__