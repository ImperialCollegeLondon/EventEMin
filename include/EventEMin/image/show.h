#ifndef EVENT_EMIN_IMAGE_SHOW_H
#define EVENT_EMIN_IMAGE_SHOW_H

#include <cassert>
#include <opencv2/highgui/highgui.hpp>
#include <string>

#include "EventEMin/image/conversion.h"
#include "EventEMin/types_def.h"

namespace EventEMin
{
template <typename Derived>
void
showColour(const DenseBase<Derived>& c, const StdVector<cv::Vec3b>& colour,
           const int width, const int height, CvMatrix& img,
           const std::string& imgName = "Image of Points")
{
  assert(2 <= c.rows());
  assert(c.cols() == static_cast<int>(colour.size()));

  eigen2cv(c, colour, width, height, img);
  cv::imshow(imgName, img);
}

template <typename Derived>
void
showColour(const DenseBase<Derived>& c, const StdVector<cv::Vec3b>& colour,
           const int width, const int height,
           const std::string& imgName = "Image of Points")
{
  CvMatrix img;
  showColour(c, colour, width, height, img, imgName);
}

template <typename Derived>
void
showGray(const DenseBase<Derived>& c, const int width, const int height,
         const std::string& imgName = "Image of Points")
{
  assert(c.rows() >= 2);

  CvMatrix img;
  eigen2cv(c, width, height, img);
  cv2gray(img);
  cv::imshow(imgName, img);
}

template <typename DerivedC, typename DerivedW>
void
showGray(const DenseBase<DerivedC>& c, const DenseBase<DerivedW>& w,
         const int width, const int height,
         const std::string& imgName = "Image of Points")
{
  assert(c.rows() >= 2);
  assert(c.cols() == w.size());

  CvMatrix img;
  eigen2cv(c, w, width, height, img);
  cv2gray(img);
  cv::imshow(imgName, img);
}

template <typename Derived>
void
showGray(const DenseBase<Derived>& c,
         const std::string& imgName = "Image of Points")
{
  typedef typename Derived::Scalar T;

  const CvMatrix img(c.cols(), c.rows(), CV_TYPE(T, 1), const_cast<T*>(&c(0)));
  CvMatrix imgGray;
  cv2gray(img, imgGray);
  cv::imshow(imgName, imgGray);
}

void
showOpticalFlowDirection(const int width, const int height,
                         const std::string& imgName = "Optical Flow Direction")
{
  const int w = (width >> 1) - 1, h = (height >> 1) - 1;

  CvMatrix flow[2];
  flow[0].create(height, width, CV_32F);
  flow[1].create(height, width, CV_32F);
  for (int j = 0; j < height; ++j)
  {
    for (int i = 0; i < width; ++i)
    {
      flow[0].at<float>(j, i) = static_cast<float>(w - i);
      flow[1].at<float>(j, i) = static_cast<float>(h - j);
    }
  }

  CvMatrix bgr;
  flow2colour(flow[0], flow[1], bgr, 1.0e6);
  cv::imshow(imgName, bgr);
}
}  // namespace EventEMin

#endif  // EVENT_EMIN_IMAGE_SHOW_H
