#ifndef EVENT_TYPE_H
#define EVENT_TYPE_H

#include <iostream>

#include "types_def.h"

namespace event_model
{
template <typename T, int N>
struct Event
{
  Vector<T, N> c;
  T ts;
  int polarity;

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  Event(void) {}
  Event(const Ref<const Vector<T, N> >& c, const T& ts, const int polarity)
      : c(c), ts(ts), polarity(polarity)
  {
  }

  friend std::ostream&
  operator<<(std::ostream& os, const Event<T, N>& ev)
  {
    return os << ev.ts << ' ' << ev.c.transpose() << ' ' << ev.polarity;
  }
  friend std::istream&
  operator>>(std::istream& is, Event<T, N>& ev)
  {
    is >> ev.ts;
    for (int i = 0; i < ev.c.size(); ++i)
    {
      is >> ev.c(i);
    }
    is >> ev.polarity;
    if (ev.polarity == 0)
    {
      ev.polarity = -1;
    }
    return is;
  }
};

template <typename T, int N>
using Events = StdVector<Event<T, N> >;
}  // namespace event_model

#endif  // EVENT_TYPE_H
