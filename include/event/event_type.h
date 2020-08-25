#ifndef EVENT_TYPE_H
#define EVENT_TYPE_H

#include <Eigen/Core>
#include <iostream>
#include <vector>

#include "types_def.h"

namespace event
{

template<typename T, int N=Eigen::Dynamic>
struct Event
{

  vec<T, N> c;
  T ts;
  int polarity;

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  Event(void)
  {}
  Event(
    const Eigen::Ref<const vec<T, N> >& c, const T& ts, const int polarity):
    c(c), ts(ts), polarity(polarity)
  {}

  friend std::ostream&
  operator<<(
    std::ostream& os, const Event<T, N>& ev)
  {
    return os<<ev.ts<<' '<<ev.c.transpose()<<' '<<ev.polarity;
  }
  friend std::istream&
  operator>>(
    std::istream& is, Event<T, N>& ev)
  {
    is>>ev.ts;
    for(int i=0; i<ev.c.size(); ++i)
    {
      is>>ev.c(i);
    }
    is>>ev.polarity;
    if(ev.polarity!=1)
    {
      ev.polarity=-1;
    }
    return is;
  }

};

template<typename T, int N=Eigen::Dynamic>
using Events=typename std::vector<Event<T, N> >;
template<typename T, int N=Eigen::Dynamic>
using EventsIter=typename std::vector<Event<T, N> >::iterator;
template<typename T, int N=Eigen::Dynamic>
using EventsConstIter=typename std::vector<Event<T, N> >::const_iterator;

}

#endif // EVENT_TYPE_H
