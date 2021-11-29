#ifndef EVENT_EMIN_EVENT_CONVERSION_H
#define EVENT_EMIN_EVENT_CONVERSION_H

#include <cassert>
#include <cmath>

#include "EventEMin/event/type.h"
#include "EventEMin/types_def.h"

namespace EventEMin
{
template <typename T, int N>
void
events2cv(const Events<T, N>& evs, const int width, const int height,
          CvMatrix& img)
{
  assert(N >= 2);
  const int nEvents = static_cast<int>(evs.size());

  img.create(height, width, CV_TYPE(T, 1));
  img.setTo(0.0);
  for (int k = 0; k < nEvents; ++k)
  {
    const int x = std::round(evs[k].c(0));
    const int y = std::round(evs[k].c(1));
    if (0 <= x && x < width && 0 <= y && y < height)
    {
      img.at<T>(y, x) += static_cast<T>(evs[k].polarity);
    }
  }
}

template <typename T, int N>
void
events2cv(const Events<T, N>& evs, const Ref<const Vector<T> >& w,
          const int width, const int height, CvMatrix& img)
{
  assert(N >= 2);
  const int nEvents = static_cast<int>(evs.size());
  assert(nEvents == w.size());

  img.create(height, width, CV_TYPE(T, 1));
  img.setTo(0.0);
  for (int k = 0; k < nEvents; ++k)
  {
    const int x = std::round(evs[k].c(0));
    const int y = std::round(evs[k].c(1));
    if (0 <= x && x < width && 0 <= y && y < height)
    {
      img.at<T>(y, x) += w(k) * static_cast<T>(evs[k].polarity);
    }
  }
}

template <typename Derived, int N>
void
event2eigen(const Event<typename Derived::Scalar, N>& ev,
            const DenseBase<Derived>& c, typename Derived::Scalar& ts)
{
  const_cast<DenseBase<Derived>&>(c) << ev.c;
  ts = ev.ts;
}

template <typename Derived, int N>
void
event2eigen(const Event<typename Derived::Scalar, N>& ev,
            const DenseBase<Derived>& c, typename Derived::Scalar& ts,
            int& polarity)
{
  event2eigen(ev, c, ts);
  polarity = ev.polarity;
}

template <typename T, int N>
void
events2eigen(const Events<T, N>& evs, Matrix<T>& c, Vector<T>& ts)
{
  const int nEvents = static_cast<int>(evs.size());
  assert(nEvents > 0);

  c.resize(N, nEvents);
  ts.resize(nEvents);
  for (int k = 0; k < nEvents; ++k)
  {
    event2eigen(evs[k], c.col(k), ts(k));
  }
}

template <typename T, int N>
void
events2eigen(const Events<T, N>& evs, Matrix<T>& c, Vector<T>& ts,
             Vector<int>& polarity)
{
  const int nEvents = static_cast<int>(evs.size());
  assert(nEvents > 0);

  c.resize(N, nEvents);
  ts.resize(nEvents);
  polarity.resize(nEvents);
  for (int k = 0; k < nEvents; ++k)
  {
    event2eigen(evs[k], c.col(k), ts(k), polarity(k));
  }
}
}  // namespace EventEMin

#endif  // EVENT_EMIN_EVENT_CONVERSION_H
