#ifndef EVENT_SHOW_H
#define EVENT_SHOW_H

#include <cassert>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/opencv.hpp>
#include <string>

#include "event/event_type.h"
#include "event/event_conversion.h"
#include "image/image_conversion.h"

namespace event
{

template<typename T>
void
showGray(
	const Events<T, 2>& evs,
	const int width, const int height,
	const std::string& imgName="Image of Events")
{
	cv::Mat img;
	events2cv<T>(evs, width, height, img);
	cv2gray(img);
	cv::imshow(imgName, img);
}

template<typename T>
void
showGray(
  const Events<T, 2>& evs,
  const vec<T>& w,
  const int width, const int height,
  const std::string& imgName="Image of Events")
{
  assert(static_cast<int>(evs.size())==w.size());

  cv::Mat img;
  events2cv<T>(evs, w, width, height, img);
  cv2gray(img);
  cv::imshow(imgName, img);
}

}

#endif // EVENT_SHOW_H
