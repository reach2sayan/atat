#ifndef __LINKEDLIST_STL_H__
#define __LINKEDLIST_STL_H__

#include <list>

template <class T>
class LinkedListIterator;

template <class T>
class LinkedList {
 public:
  using iterator_category = std::forward_iterator_tag;
  using difference_type = std::ptrdiff_t;
  using value_type = T;
  using pointer = T*;
  using reference = T&;
  using iterator = LinkedListIterator<T>;
  using const_iterator = LinkedListIterator<const T>;

 private:
  std::list<T> mylist;

 public:
  explicit LinkedList() = default;
  explicit LinkedList(std::size_t n) : mylist{n} {}
  LinkedList(std::size_t n, const T& val) : mylist{n, val} {}

  template <class InputIterator>
  LinkedList(InputIterator first, InputIterator last) : mylist{first, last} {}
  LinkedList(const LinkedList& x) : mylist{x} {}
  LinkedList(LinkedList&& x) : mylist{x} {}

  LinkedList(std::initializer_list<T> il) : mylist{il} {}

  ~LinkedList() = default;

  LinkedList& operator=(const LinkedList& x) { mylist = x; }
  LinkedList& operator=(LinkedList&& x) { mylist = x; }
  LinkedList& operator=(std::initializer_list<T> il) { mylist = il; }

  LinkedListIterator<T> begin() noexcept;
  const LinkedListIterator<T> begin() const noexcept;

  T& front() { return mylist.front(); }
  const T& front() const { return mylist.front(); }
  T& back() { return mylist.back(); }
  const T& back() const { return mylist.back(); }

  int size() const noexcept { return mylist.size(); }

  void push_front(const T& val) { mylist.push_front(val); }
  void push_front(T&& val) { mylist.push_front(val); }

  void push_back(const T& val) { mylist.push_back(val); }
  void push_back(T&& val) { mylist.push_back(val); }

  void clear() noexcept { mylist.clear(); }

  LinkedListIterator<T> erase(LinkedListIterator<const T> position) {
    mylist.erase(position);
  }
  LinkedListIterator<T> erase(LinkedListIterator<const T> first,
			      LinkedListIterator<const T> last) {
    mylist.erase(first, last);
  }

  LinkedListIterator<T> insert(LinkedListIterator<const T> position,
			       const T& val);
  LinkedListIterator<T> insert(LinkedListIterator<const T> position, int n,
			       const T& val);

  template <class InputIterator>
  LinkedListIterator<T> insert(LinkedListIterator<const T> position,
			       InputIterator first, InputIterator last) {
    mylist.insert(position, first, last);
  }

  LinkedListIterator<T> insert(LinkedListIterator<const T> position, T&& val) {
    mylist.insert(position, val);
  }

  LinkedListIterator<T> insert(LinkedListIterator<const T> position,
			       std::initializer_list<T> il) {
    mylist.insert(position, il);
  }

  // legacy ATAT functions
  T* push_front(T* newobj) {
    mylist.push_front(*newobj);
    T* retval = mylist.front();
    return retval;
  }

  T* detach(LinkedListIterator<T> deleteat) {
    if (mylist.empty()) return nullptr;
    T* obj = *deleteat;
    mylist.erase(deleteat);
    return obj;
  }

  T* push_back(T* newobj) {
    mylist.push_back(*newobj);
    T* retval = &(mylist.back());
    return retval;
  }

  T* add(T* newobj, LinkedListIterator<T> insertat) {
    LinkedListIterator<T> retobj = mylist.insert(insertat, *newobj);
    if (retobj == mylist.end())
      return nullptr;
    else
      return &(*retobj);
  }

  int getSize() const { return mylist.size(); }
  int get_size() const { return mylist.size(); }
  void delete_all() { mylist.clear(); }
};

template <class T>
class LinkedListIterator {
 public:
  using iterator_category = typename LinkedList<T>::iterator_category;
  using difference_type = typename LinkedList<T>::difference_type;
  using value_type = typename LinkedList<T>::value_type;
  using pointer = typename LinkedList<T>::pointer;
  using reference = typename LinkedList<T>::reference;
  using iterator = typename LinkedList<T>::iterator;
  using const_iterator = typename LinkedList<T>::const_iterator;

 private:
  T* m_ptr;

 public:
  LinkedListIterator() : m_ptr(nullptr) {}
  LinkedListIterator(T* ptr) : m_ptr(ptr) {}

  T& operator*() const { return *m_ptr; }
  T* operator->() const { return m_ptr; }

  LinkedListIterator<T>& operator++() {
    m_ptr++;
    return *this;
  }
  LinkedListIterator<T> operator++(int) {
    LinkedListIterator<T> tmp = *this;
    ++(*this);
    return tmp;
  }
  friend bool operator==(const LinkedListIterator<T> a,
			 const LinkedListIterator<T> b) {
    return a.m_ptr == b.m_ptr;
  }
  friend bool operator!=(const LinkedListIterator<T> a,
			 const LinkedListIterator<T> b) {
    return a.m_ptr != b.m_ptr;
  }
};

#endif
