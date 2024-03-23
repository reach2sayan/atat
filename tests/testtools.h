#ifndef __ATAT_TESTTOOLS_HPP__
#define __ATAT_TESTTOOLS_HPP__

#include <initializer_list>
#include <type_traits>

namespace ATATTestTools {

template <typename T>
struct reference_t_matches_deref_t
    : std::is_same<typename T::reference, decltype(*std::declval<T&>())> {};
// Can't be momoved, assigned, copied or default constructed
template <typename T>
class MonolithObject {
 private:
  const T value;

 public:
  constexpr MonolithObject(T v) : value{v} {}
  constexpr T get() const { return value; }

  MonolithObject() = delete;
  MonolithObject(const MonolithObject&) = delete;
  MonolithObject(MonolithObject&&) = delete;
  MonolithObject& operator=(const MonolithObject&) = delete;
  MonolithObject& operator=(MonolithObject&&) = delete;
  ~MonolithObject() = default;
};

// IterableType provides a type which allows forward iterator
// operator++(), operator!=(const IterableType&), operator*()
// move constructible only
// not copy constructible, move assignable, or copy assignable
template <typename T>
class IterableType {
 private:
  T* data;
  std::size_t size;
  bool was_moved_from_ = false;
  mutable bool was_copied_from_ = false;

 public:
  IterableType() = default;
  IterableType(std::initializer_list<T> il)
      : data{new T[il.size()]}, size{il.size()} {
    std::size_t i = 0;
    for (auto&& e : il) {
      data[i] = e;
      ++i;
    }
  }

  IterableType& operator=(IterableType&&) = delete;
  IterableType& operator=(const IterableType&) = delete;

#ifndef DEFINE_BASIC_ITERABLE_COPY_CTOR
  IterableType(const IterableType&) = delete;
#else
  IterableType(const IterableType& other)
      : data{new T[other.size]}, size{other.size} {
    other.was_copied_from_ = true;
    auto o_it = begin(other);
    for (auto it = begin(*this); o_it != end(other); ++it, ++o_it) {
      *it = *o_it;
    }
  }
#endif

  IterableType(IterableType&& other) : data{other.data}, size{other.size} {
    other.data = nullptr;
    other.was_moved_from_ = true;
  }

  bool was_moved_from() const { return this->was_moved_from_; }
  bool was_copied_from() const { return this->was_copied_from_; }
  ~IterableType() { delete[] this->data; }
  template <typename U>
  class Iterator {
   private:
    U* p;

   public:
#ifdef DEFINE_DEFAULT_ITERATOR_CTOR
    Iterator() = default;
#endif
    Iterator(U* b) : p{b} {}
    bool operator!=(const Iterator& other) const { return this->p != other.p; }
    Iterator& operator++() {
      ++this->p;
      return *this;
    }
    U& operator*() { return *this->p; }
  };

  IterableType::Iterator<T> begin() { return {data}; }
  IterableType::Iterator<T> end() { return {data + size}; }

#ifdef DEFINE_BASIC_ITERABLE_CONST_BEGIN_AND_END
  IterableType::Iterator<const T> begin(const IterableType& b) {
    return {b.data};
  }

  IterableType::Iterator<const T> end(const IterableType& b) {
    return {b.data + b.size};
  }
#endif

#ifdef DECLARE_REVERSE_ITERATOR
  Iterator<T> rbegin();
  Iterator<T> rend();
#endif	// ifdef DECLARE_REVERSE_ITERATOR
};

}  // namespace ATATTestTools

#endif	//__ATAT_TESTTOOLS_HPP__
