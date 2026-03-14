#include "kmeci.h"

namespace {

template <typename T>
void accumulate_nested(T *dest, const T &src) {
  *dest += src;
}

template <typename T>
void accumulate_nested(Array<T> *dest, const Array<T> &src) {
  if (dest->get_size() != src.get_size()) {
    dest->resize(src.get_size());
  }
  for (int i = 0; i < src.get_size(); i++) {
    accumulate_nested(&((*dest)(i)), src(i));
  }
}

}  // namespace

void MultiKSpaceECI::get_k_space_eci(
    Array<Array<Array<Array<Array<Complex>>>>> *p_ft_eci,
    const Array<Real> &x) {
  LinkedListIterator<KSpaceECI> iter(*this);
  if (!iter) {
    p_ft_eci->resize(0);
    return;
  }

  iter->get_k_space_eci(p_ft_eci, x);
  iter++;

  for (; iter; iter++) {
    Array<Array<Array<Array<Array<Complex>>>>> current;
    iter->get_k_space_eci(&current, x);
    accumulate_nested(p_ft_eci, current);
  }
}
