#ifndef IMAGE_CONVERSION_H
#define IMAGE_CONVERSION_H

#include <cassert>
#include <Eigen/Core>
#include <opencv2/opencv.hpp>

#include "types_def.h"

template<typename T, typename D>
void
eigen2cv(
	const D& c,
	const int width, const int height,
	cv::Mat& img)
{
  assert(c.rows()==2);

	img.create(height, width, CV_MAKETYPE(cv::DataDepth<T>::value, 1));
  img.setTo(0.0);
  const T one=T(1.0);
  for(int k=0; k<c.cols(); ++k)
  {
  	const int x=static_cast<int>(std::round(c(0, k)));
  	const int y=static_cast<int>(std::round(c(1, k)));
		if(0<=x && x<width && 0<=y && y<height)
    {
      img.at<T>(y, x)+=one;
    }  
  }
}

template<typename T, typename D1, typename D2>
void
eigen2cv(
  const D1& c,
  const D2& w,
  const int width, const int height,
  cv::Mat& img)
{
  assert(c.rows()==2);
  assert(c.cols()==w.size());

  img.create(height, width, CV_MAKETYPE(cv::DataDepth<T>::value, 1));
  img.setTo(0.0);
  for(int k=0; k<c.cols(); ++k)
  {
    const int x=static_cast<int>(std::round(c(0, k)));
    const int y=static_cast<int>(std::round(c(1, k)));
    if(0<=x && x<width && 0<=y && y<height)
    {
      img.at<T>(y, x)+=static_cast<T>(w(k));
    }  
  }
}

void
cv2gray(
  const cv::Mat& img,
  cv::Mat& imgGray)
{
  double min, max;
  cv::minMaxLoc(img, &min, &max);
  const double a=255.0/(max-min+1.0e-32), b=-a*min;
  img.convertTo(imgGray, CV_8UC1, a, b);
}

void
cv2gray(
  cv::Mat& img)
{
	cv2gray(img, img);
}

#endif // IMAGE_CONVERSION_H
