#ifndef EVENT_EMIN_EVENT_SHOW_H
#define EVENT_EMIN_EVENT_SHOW_H

#include <cassert>
#include <opencv2/highgui/highgui.hpp>
#include <string>

#include "EventEMin/event/conversion.h"
#include "EventEMin/event/type.h"
#include "EventEMin/image/conversion.h"

namespace EventEMin
{
template <typename T, int N>
void
showGray(const Events<T, N>& evs, const int width, const int height,
         const std::string& imgName = "Image of Events")
{
  assert(N >= 2);
  CvMatrix img;
  events2cv<T>(evs, width, height, img);
  cv2gray(img);
  cv::imshow(imgName, img);
}

template <typename T, int N, typename Derived>
void
showGray(const Events<T, N>& evs, const DenseBase<Derived>& w, const int width,
         const int height, const std::string& imgName = "Image of Events")
{
  assert(N >= 2);
  assert(static_cast<int>(evs.size()) == w.size());

  CvMatrix img;
  events2cv<T>(evs, w, width, height, img);
  cv2gray(img);
  cv::imshow(imgName, img);
}
}  // namespace EventEMin

#endif  // EVENT_EMIN_EVENT_SHOW_H
