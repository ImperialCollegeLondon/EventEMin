#ifndef EVENT_UNDISTORT_H
#define EVENT_UNDISTORT_H

#include <cmath>
#include <fstream>
#include <opencv2/opencv.hpp>

#include "event/event_type.h"

namespace event
{

template<typename T>
void
undistort(
	const int xMin, const int xMax,
	const int yMin, const int yMax,
	const cv::Mat& map,
	const Events<T, 2>& evs,
	Events<T, 2>& uevs)
{
	uevs.clear();
	uevs.reserve(evs.size());

	Event<T, 2> ev;
	for(int k=0; k<static_cast<int>(evs.size()); ++k)
	{
		cv::Vec<T, 2> mapPoint=map.at<cv::Vec<T, 2> >(
			static_cast<int>(std::round(evs.at(k).c(1))),
			static_cast<int>(std::round(evs.at(k).c(0))));
		if(xMin<=static_cast<int>(mapPoint[0]) && static_cast<int>(mapPoint[0])<=xMax &&
      yMin<=static_cast<int>(mapPoint[1]) && static_cast<int>(mapPoint[1])<=yMax)
    {
    	ev.c(0)=mapPoint[0];
    	ev.c(1)=mapPoint[1];
    	ev.ts=evs.at(k).ts;
    	ev.polarity=evs.at(k).polarity;
    	uevs.push_back(ev);
    }
	}
}

template<typename T>
void
undistort(
  const int xMin, const int xMax,
  const int yMin, const int yMax,
  const cv::Mat& map,
  const int nEvents,
  std::ifstream& fin,
  Events<T, 2>& uevs)
{
  uevs.clear();
  uevs.reserve(nEvents);

  Event<T, 2> ev;
  while(static_cast<int>(uevs.size())<nEvents && fin>>ev)
  {
    cv::Vec<T, 2> mapPoint=map.at<cv::Vec<T, 2> >(
      static_cast<int>(std::round(ev.c(1))),
      static_cast<int>(std::round(ev.c(0))));
    if(xMin<=static_cast<int>(mapPoint[0]) && static_cast<int>(mapPoint[0])<=xMax &&
      yMin<=static_cast<int>(mapPoint[1]) && static_cast<int>(mapPoint[1])<=yMax)
    {
      ev.c(0)=mapPoint[0];
      ev.c(1)=mapPoint[1];
      uevs.push_back(ev);
    }
  }
}

template<typename T>
void
undistort(
  const int xMin, const int xMax,
  const int yMin, const int yMax,
  const cv::Mat& map,
  const int nEvents,
  std::ifstream& fin,
  mtx<T>& c,
  vec<T>& ts,
  vec<int>& polarity)
{
  c.resize(2, nEvents);
  ts.resize(nEvents);
  polarity.resize(nEvents);

  int n=0;
  Event<T, 2> ev;
  while(n<nEvents && fin>>ev)
  {
    cv::Vec<T, 2> mapPoint=map.at<cv::Vec<T, 2> >(
      static_cast<int>(std::round(ev.c(1))),
      static_cast<int>(std::round(ev.c(0))));
    if(xMin<=static_cast<int>(mapPoint[0]) && static_cast<int>(mapPoint[0])<=xMax &&
      yMin<=static_cast<int>(mapPoint[1]) && static_cast<int>(mapPoint[1])<=yMax)
    {
      c.col(n)<<mapPoint[0], mapPoint[1];
      ts(n)=ev.ts;
      polarity(n)=ev.polarity;
      ++n;
    }
  }
}

}

#endif // EVENT_UNDISTORT_H
