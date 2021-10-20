#ifndef EVENT_UNDISTORT_H
#define EVENT_UNDISTORT_H

#include <cassert>
#include <cmath>
#include <fstream>
#include <opencv2/opencv.hpp>

#include "event/type.h"
#include "types_def.h"

namespace event_model
{
template <typename T, int N>
void
undistort(const int xMin, const int xMax, const int yMin, const int yMax,
          const CvMatrix& map, const Events<T, N>& evs, Events<T, N>& uevs)
{
  const int nEvents = static_cast<int>(evs.size());
  assert(N >= 2);

  uevs.clear();
  uevs.reserve(evs.size());

  Event<T, N> ev;
  for (int k = 0; k < nEvents; ++k)
  {
    CvVector<T, 2> mapPoint = map.at<CvVector<T, 2> >(std::round(evs[k].c(1)),
                                                      std::round(evs[k].c(0)));
    const int x = std::round(mapPoint[0]);
    const int y = std::round(mapPoint[1]);
    if (xMin <= x && x < xMax && yMin <= y && y < yMax)
    {
      ev = evs[k];
      ev.c(0) = mapPoint[0];
      ev.c(1) = mapPoint[1];
      uevs.push_back(ev);
    }
  }
}
}  // namespace event_model

#endif  // EVENT_UNDISTORT_H
