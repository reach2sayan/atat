#ifndef __ITERTOOLS_RANGE_HPP__
#define __ITERTOOLS_RANGE_HPP__

#include "itertoolsbase.h"

namespace itertools {
template <typename T>
class Range;

template <typename T>
constexpr Range<T> range(T) noexcept;
template <typename T>
constexpr Range<T> range(T, T) noexcept;
template <typename T>
constexpr Range<T> range(T, T, T) noexcept;

template <typename T, bool IsFloat = std::is_floating_point<T>::value>
class RangeIteratorData;

// everything except floats
template <typename T>
class RangeIteratorData<T, false> {
 private:
  T _value{};
  T _step{};

 public:
  constexpr RangeIteratorData() noexcept = default;
  constexpr RangeIteratorData(T value, T step) noexcept
      : _value{value}, _step{step} {}

  constexpr T value() const noexcept { return _value; }
  constexpr T step() const noexcept { return _step; }

  void operator++() noexcept { _value += _step; }

  constexpr bool operator==(const RangeIteratorData& other) const noexcept {
    return _value == other._value;
  }
  constexpr bool operator!=(const RangeIteratorData& other) const noexcept {
    return !(*this == other);
  }
};

// float data
template <typename T>
class RangeIteratorData<T, true> {
 private:
  T _start{};
  T _value{};
  T _step{};
  std::size_t _steps_taken{};

 public:
  constexpr RangeIteratorData() noexcept = default;
  constexpr RangeIteratorData(T start, T step) noexcept
      : _start{start}, _value{start}, _step{step} {}

  constexpr T value() const noexcept { return _value; }

  constexpr T step() const noexcept { return _step; }

  void operator++() noexcept {
    ++_steps_taken;
    _value = _start + (_step * _steps_taken);
  }

  constexpr bool operator==(const RangeIteratorData& other) const noexcept {
    // if the difference between the two values is less than the
    // step_ size, they are considered equal
    return (_value < other._value ? other._value - _value
				  : _value - other._value) < _step;
  }

  constexpr bool operator!=(const RangeIteratorData& other) const noexcept {
    return !(*this == other);
  }
};

}  // namespace itertools

template <typename T>
class itertools::Range {
  // see stackoverflow.com/questions/32174186 about why only specializations
  // aren't marked as friend
  template <typename U>
  friend constexpr Range<U> itertools::range(U) noexcept;
  template <typename U>
  friend constexpr Range<U> itertools::range(U, U) noexcept;
  template <typename U>
  friend constexpr Range<U> itertools::range(U, U, U) noexcept;

 private:
  const T _start;
  const T _stop;
  const T _step;

  constexpr Range(T stop) noexcept : _start{0}, _stop{stop}, _step{1} {}
  constexpr Range(T start, T stop, T step = 1) noexcept
      : _start{start}, _stop{stop}, _step{step} {}

  // if val is "before" the stopping point.
  static constexpr bool is_within_range(T val, T stop_val,
					[[maybe_unused]] T step_val) {
    if constexpr (std::is_unsigned<T>{}) {
      return val < stop_val;
    } else {
      return !(step_val > 0 && val >= stop_val) &&
	     !(step_val < 0 && val <= stop_val);
    }
  }

 public:
  constexpr T start() const noexcept { return _start; }
  constexpr T stop() const noexcept { return _stop; }
  constexpr T step() const noexcept { return _step; }
  constexpr T operator[](std::size_t index) const noexcept {
    return _start + (_step * index);
  }

  constexpr std::size_t size() const noexcept {
    static_assert(!std::is_floating_point_v<T>,
		  "range size() not supperted with floating point types");
    if (!is_within_range(_start, _stop, _step)) {
      return 0;
    }

    auto diff = _stop - _start;
    auto res = (_stop - _start) / _step;
    assert(res >= 0);
    auto result = static_cast<std::size_t>(res);
    if (diff % step()) {
      ++result;
    }
    return result;
  }

  // the reference type here is T, which doesn't strictly follow all
  // of the rules, but std::vector<bool>::iterator::reference isn't
  // a reference type either, this isn't any worse

  class RangeIterator {
   private:
    itertools::RangeIteratorData<T> data;
    bool _is_end{};

    // first argument must be regular iterator
    // second argument must be end iterator
    static bool not_equal_to_impl(const RangeIterator& lhs,
				  const RangeIterator& rhs) noexcept {
      assert(!lhs._is_end);
      assert(rhs._is_end);
      return is_within_range(lhs.data.value(), rhs.data.value(),
			     lhs.data.step());
    }

    static bool not_equal_to_end(const RangeIterator& lhs,
				 const RangeIterator& rhs) noexcept {
      if (rhs._is_end) {
	return not_equal_to_impl(lhs, rhs);
      }
      return not_equal_to_impl(rhs, lhs);
    }

   public:
    using iterator_category = std::forward_iterator_tag;
    using value_type = T;
    using difference_type = std::ptrdiff_t;
    using pointer = value_type*;
    using reference = value_type;

    constexpr RangeIterator() noexcept = default;

    constexpr RangeIterator(T value, T step, bool is_end) noexcept
	: data(value, step), _is_end{is_end} {}

    constexpr T operator*() const noexcept { return data.value(); }

    constexpr ArrowProxy<T> operator->() const noexcept { return {**this}; }

    RangeIterator& operator++() noexcept {
      data.operator++();
      return *this;
    }

    RangeIterator operator++(int) noexcept {
      auto ret = *this;
      ++(*this);
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
    bool operator!=(const RangeIterator& other) const noexcept {
      if (_is_end && other._is_end) {
	return false;
      }

      if (!_is_end && !other._is_end) {
	return data != other.data;
      }
      return not_equal_to_end(*this, other);
    }

    bool operator==(const RangeIterator& other) const noexcept {
      return !(*this != other);
    }
  };

  constexpr RangeIterator begin() const noexcept {
    return {_start, _step, false};
  }
  constexpr RangeIterator end() const noexcept { return {_stop, _step, true}; }
};

template <typename T>
constexpr itertools::Range<T> itertools::range(T stop) noexcept {
  return {stop};
}

template <typename T>
constexpr itertools::Range<T> itertools::range(T start, T stop) noexcept {
  return {start, stop};
}

template <typename T>
constexpr itertools::Range<T> itertools::range(T start, T stop,
					       T step) noexcept {
  return step == T(0) ? Range<T>{0} : Range<T>{start, stop, step};
}
#endif	//__ITERTOOLS_RANGE_HPP__

