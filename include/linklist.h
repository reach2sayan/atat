#ifndef __LINKEDLIST_H__
#define __LINKEDLIST_H__

#include "misc.h"
#include <iostream>

template<class T>
class LinkedList;

template<class T>
class LinkedListIterator;

template<class T>
class LinkedListNode {
 protected:
  friend class LinkedList<T>;
  friend class LinkedListIterator<T>;

  LinkedListNode<T> *next;
  T *obj;
 public:
  LinkedListNode(T *newobj=NULL,LinkedListNode<T> *newnext=NULL) {
    obj=newobj;
    next=newnext;
  }
};

template<class T>
class LinkedList {
 protected:
  friend class LinkedListIterator<T>;
  
  LinkedListNode<T> *head;
  LinkedListNode<T> *tail;
public:
  LinkedList(void) {
    head=new LinkedListNode<T>(NULL,NULL);
    tail=head;
  }
  LinkedList(const LinkedList<T> &ll) {
    head=new LinkedListNode<T>(NULL,NULL);
    tail=head;
    LinkedListIterator<T> i(ll);
    for ( ;i; i++) {
      (*this) << new T(*i);
    }
  }
  ~LinkedList() {
    LinkedListIterator<T> i(*this);
    while (i)
      delete detach(i);
    delete head;
  }
  T *push_back(T *newobj) {
    tail->next=new LinkedListNode<T>(newobj,tail->next);
    tail=tail->next;
    return newobj;
  }
  LinkedList<T>& operator<< (T *newobj) {
    push_back(newobj);
    return *this;
  }
  T *add(T *newobj, LinkedListIterator<T> &insertat) {
    if (insertat.previous) {
      insertat.previous->next=new LinkedListNode<T>(newobj,insertat.previous->next);
      if(insertat.previous==tail)
	tail=insertat.previous->next;
      return newobj;
    }
    else
      return NULL;
  }
  T *push_front(T *newobj) {
    LinkedListIterator<T> i(*this);
    return add(newobj,i);
  }
  T *detach(LinkedListIterator<T> &deleteat) {
    LinkedListNode<T> *previous=deleteat.previous;
    if (previous && previous->next) {
      T *detachedobj=previous->next->obj;
      LinkedListNode<T> *detachednode=previous->next;
      if( deleteat.get_current() == tail )
	tail = previous;
      previous->next=previous->next->next;
      delete detachednode;
      return detachedobj;
    }
    else
      return NULL;
  }
  int getSize(void) const {
    int nb=0;
    LinkedListIterator<T> i(*this);
    while (i) {
      nb++;
      i++;
    }
    return nb;
  }
  int get_size(void) const {
    return getSize();
  }
  void delete_all(void) {
    LinkedListIterator<T> i(*this);
    while (i) delete detach(i);
  }
  void detach_all(void) {
    LinkedListIterator<T> i(*this);
    while (i) detach(i);
  }
};

template <class T>
class LinkedListIterator {
 protected:
  friend class LinkedList<T>;

  LinkedListNode<T> *previous;

  LinkedListNode<T> *get_current(void) {return (previous ? previous->next : NULL);}

 public:
  LinkedListIterator(void) {
    previous=NULL;
  }
  LinkedListIterator(const LinkedList<T> &linkedlist) {
    init(linkedlist);
  }
  void init(const LinkedList<T> &linkedlist) {
    previous=linkedlist.head;
  }
  LinkedListIterator(const LinkedListIterator<T> &iter) {
    previous=iter.previous;
  }
  void operator=(const LinkedListIterator<T> &iter) {
    previous=iter.previous;
  }
  LinkedListIterator<T>& operator++(int) {
    if (get_current()){
      previous=get_current();
    }
    return *this;
  }
  LinkedListIterator<T> operator+(int i) {
    LinkedListIterator<T> iter=*this;
    while (i--) iter++;
    return iter;
  }
  T& operator *(void) {
    return *(get_current()->obj);
  }
  operator T * () {
    if (get_current())
      return get_current()->obj;
    else
      return NULL;
  }
  T* operator -> () {
    if (get_current())
      return get_current()->obj;
    else
      return NULL;
  }
  T* setObject(T *newvalue) {
    if (get_current()) {
      get_current()->obj=newvalue;
      return newvalue;
    } else
      return NULL;
  }
};

// utilitary functions

template<class T>
void add_copy(LinkedList<T> *dest, const LinkedList<T> &src) {
  LinkedListIterator<T> i(src);
  while (i) {
    (*dest) << new T(*i);
    i++;
  }
}

template<class T>
void transfer_list(LinkedList<T> *dest, LinkedList<T> *src) {
  LinkedListIterator<T> i(*src);
  while (i) {
    (*dest) << src->detach(i);
  }
}

template<class T>
inline void add_at_beginning(LinkedList<T> *ll,T *obj) {
  LinkedListIterator<T> beginning(*ll);
  ll->add(obj,beginning);
}

template<class T>
ostream& operator << (ostream &s, const LinkedList<T> &a) {
  s << a.get_size() << endl;
  LinkedListIterator<T> i(a);
  for (; i; i++) {
    s << *i << endl;
  }
  return s;
}

template<class T>
istream& operator >> (istream &s, LinkedList<T> &a) {
  int size=0;
  s >> size;
  for (int i=0; i<size; i++) {
    T *p=new T;
    s >> (*p);
    a << p;
  }
  return s;
}

#endif
