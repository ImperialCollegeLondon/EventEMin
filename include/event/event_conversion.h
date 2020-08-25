#ifndef EVENT_CONVERSION_H
#define EVENT_CONVERSION_H

#include <cassert>
#include <cmath>
#include <Eigen/Core>
#include <opencv2/opencv.hpp>

#include "event/event_type.h"
#include "types_def.h"

namespace event
{

template<typename T>
void
events2cv(
  const Events<T, 2>& evs,
  const int width, const int height,
  cv::Mat& img)
{
  img.create(height, width, CV_MAKETYPE(cv::DataDepth<T>::value, 1));
  img.setTo(0.0);
  for(int k=0; k<static_cast<int>(evs.size()); ++k)
  {
    const int x=static_cast<int>(std::round(evs.at(k).c(0)));
    const int y=static_cast<int>(std::round(evs.at(k).c(1)));
    if(0<=x && x<width && 0<=y && y<height)
    {
      img.at<T>(y, x)+=static_cast<T>(evs.at(k).polarity);
    }
  }
}

template<typename T>
void
events2cv(
  const Events<T, 2>& evs,
  const vec<T>& w,
  const int width, const int height,
  cv::Mat& img)
{
  assert(static_cast<int>(evs.size())==w.size());

  img.create(height, width, CV_MAKETYPE(cv::DataDepth<T>::value, 1));
  img.setTo(0.0);
  for(int k=0; k<static_cast<int>(evs.size()); ++k)
  {
    const int x=static_cast<int>(std::round(evs.at(k).c(0)));
    const int y=static_cast<int>(std::round(evs.at(k).c(1)));
    if(0<=x && x<width && 0<=y && y<height)
    {
      img.at<T>(y, x)+=w(k)*static_cast<T>(evs.at(k).polarity);
    }
  }
}

template<typename Derived, int N=Eigen::Dynamic>
void
event2eigen(
	const Event<typename Derived::Scalar, N>& ev,
	Eigen::DenseBase<Derived> const & c,
	typename Derived::Scalar& ts)
{
	const_cast<Eigen::DenseBase<Derived>&>(c)<<ev.c;
	ts=ev.ts;
}

template<typename Derived, int N=Eigen::Dynamic>
void
event2eigen(
	const Event<typename Derived::Scalar, N>& ev,
	Eigen::DenseBase<Derived> const & c,
	typename Derived::Scalar& ts,
	int& polarity)
{
	event2eigen(ev, c, ts);
	polarity=ev.polarity;
}

template<typename T, int N=Eigen::Dynamic>
void
events2eigen(
  const Events<T, N>& evs,
  mtx<T>& c,
  vec<T>& ts)
{
  assert(0<evs.size());
  c.resize(evs.at(0).c.size(), evs.size());
  ts.resize(evs.size());
  for(int k=0; k<static_cast<int>(evs.size()); ++k)
  {
    event2eigen(evs.at(k), c.col(k), ts(k));
  }
}

template<typename T, int N=Eigen::Dynamic>
void
events2eigen(
  const Events<T, N>& evs,
  mtx<T>& c,
  vec<T>& ts,
  vec<int>& polarity)
{
  assert(0<evs.size());
  c.resize(evs.at(0).c.size(), evs.size());
  ts.resize(evs.size());
  polarity.resize(evs.size());
  for(int k=0; k<static_cast<int>(evs.size()); ++k)
  {
    event2eigen(evs.at(k), c.col(k), ts(k), polarity(k));
  }
}

}

#endif // EVENT_CONVERSION_H
