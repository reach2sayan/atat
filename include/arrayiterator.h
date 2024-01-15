#ifndef __ARRAY_ITERATOR_H__
#define __ARRAY_ITERATOR_H__

#include <algorithm>
#include <iterator>
#include <memory>

template <class T>
class ArrayIterator {
 public:
  using iterator_category = std::random_access_iterator_tag;
  using value_type = T;
  using difference_type = ptrdiff_t;
  using pointer = T*;
  using reference = T&;

  ArrayIterator() : m_ptr(nullptr) {}
  ArrayIterator(T* ptr) : m_ptr(ptr){};
  ArrayIterator(const ArrayIterator<T>& other) : m_ptr(other.m_ptr) {}

  T& operator*() const { return *m_ptr; }
  T* operator->() const { return m_ptr; }

  T& operator[](ptrdiff_t index) const { return *(m_ptr + index); }

  ArrayIterator<T>& operator++();
  ArrayIterator<T> operator++(int) const;
  ArrayIterator<T>& operator--();
  ArrayIterator<T> operator--(int) const;
  ArrayIterator<T>& operator+=(const ptrdiff_t movement);
  ArrayIterator<T>& operator-=(const ptrdiff_t movement);

  bool operator==(const ArrayIterator<T>& other) const {
    return m_ptr == other.m_ptr;
  }
  bool operator!=(const ArrayIterator<T>& other) const {
    return m_ptr != other.m_ptr;
  }
  bool operator>(const ArrayIterator<T>& other) const {
    return m_ptr > other.m_ptr;
  }
  bool operator<(const ArrayIterator<T>& other) const {
    return m_ptr < other.m_ptr;
  }
  bool operator>=(const ArrayIterator<T>& other) const {
    return m_ptr >= other.m_ptr;
  }
  bool operator<=(const ArrayIterator<T>& other) const {
    return m_ptr <= other.m_ptr;
  }

  ArrayIterator<T> operator+(const ptrdiff_t movement) const;
  ptrdiff_t operator+(const ArrayIterator<T> other) const {
    return m_ptr + other.m_ptr;
  }

  ArrayIterator<T> operator-(const ptrdiff_t movement) const;
  ptrdiff_t operator-(const ArrayIterator<T> other) const {
    return m_ptr - other.m_ptr;
  }
  operator ArrayIterator<const T>() const;  // implicit conversion

 private:
  T* m_ptr;
};

template <class T>
inline ArrayIterator<T>::operator ArrayIterator<const T>() const {
  ArrayIterator<const T> temp(m_ptr);
  return temp;
}

template <class T>
inline ArrayIterator<T> ArrayIterator<T>::operator--(int) const {
  ArrayIterator<T> temp = *this;
  --m_ptr;
  return temp;
}

template <class T>
inline ArrayIterator<T>& ArrayIterator<T>::operator+=(
    const ptrdiff_t movement) {
  m_ptr += movement;
  return *this;
}

template <class T>
inline ArrayIterator<T>& ArrayIterator<T>::operator-=(
    const ptrdiff_t movement) {
  m_ptr -= movement;
  return *this;
}

template <class T>
inline ArrayIterator<T> ArrayIterator<T>::operator++(int) const {
  ArrayIterator<T> temp = *this;
  ++m_ptr;
  return temp;
}

template <class T>
inline ArrayIterator<T>& ArrayIterator<T>::operator--() {
  --m_ptr;
  return *this;
}

template <class T>
inline ArrayIterator<T>& ArrayIterator<T>::operator++() {
  ++m_ptr;
  return *this;
}

template <class T>
inline ArrayIterator<T> ArrayIterator<T>::operator-(ptrdiff_t movement) const {
  ArrayIterator<T> temp = *this;
  temp -= movement;
  return temp;
}

template <class T>
inline ArrayIterator<T> ArrayIterator<T>::operator+(ptrdiff_t movement) const {
  ArrayIterator<T> temp = *this;
  temp += movement;
  return temp;
}
#endif
