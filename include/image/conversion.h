#ifndef IMAGE_CONVERSION_H
#define IMAGE_CONVERSION_H

#include <cassert>

#include "types_def.h"

namespace event_model
{
template <typename Derived>
void
eigen2cv(const DenseBase<Derived>& c, const int width, const int height,
         CvMatrix& img, const double inc = -1.0)
{
  assert(c.rows() >= 2);

  typedef typename Derived::Scalar T;

  img.create(height, width, CV_TYPE(T, 1));
  img.setTo(0.0);
  for (int k = 0; k < c.cols(); ++k)
  {
    const int x = std::round(c(0, k));
    const int y = std::round(c(1, k));
    if (0 <= x && x < width && 0 <= y && y < height)
    {
      img.at<T>(y, x) += static_cast<T>(inc);
    }
  }
}

template <typename DerivedC, typename DerivedW>
void
eigen2cv(const DenseBase<DerivedC>& c, const DenseBase<DerivedW>& w,
         const int width, const int height, CvMatrix& img)
{
  assert(c.rows() >= 2);
  assert(c.cols() == w.size());

  typedef typename DerivedC::Scalar T;

  img.create(height, width, CV_TYPE(T, 1));
  img.setTo(0.0);
  for (int k = 0; k < c.cols(); ++k)
  {
    const int x = std::round(c(0, k));
    const int y = std::round(c(1, k));
    if (0 <= x && x < width && 0 <= y && y < height)
    {
      img.at<T>(y, x) += static_cast<T>(w(k));
    }
  }
}

template <typename Derived>
void
eigen2cv(const DenseBase<Derived>& c, const StdVector<cv::Vec3b>& colour,
         const int width, const int height, CvMatrix& img)
{
  assert(2 <= c.rows());
  assert(c.cols() == static_cast<int>(colour.size()));

  img.create(height, width, CV_8UC3);
  img.setTo(255);
  for (int k = 0; k < c.cols(); ++k)
  {
    const int x = std::round(c(0, k));
    const int y = std::round(c(1, k));
    if (0 <= x && x < width && 0 <= y && y < height)
    {
      // img.at<cv::Vec3b>(y, x) = colour[k];
      cv::circle(img, cv::Point(x, y), 1, colour[k], 1);
    }
  }
}

void
cv2gray(const CvMatrix& img, CvMatrix& imgGray)
{
  double min, max;
  cv::minMaxLoc(img, &min, &max);
  const double a = 255.0 / (max - min + 1.0e-32), b = -a * min;
  img.convertTo(imgGray, CV_8UC1, a, b);
}

void
cv2gray(CvMatrix& img)
{
  cv2gray(img, img);
}

void
flow2colour(const CvMatrix& flowX, const CvMatrix& flowY, CvMatrix& bgr,
            const double& thresh = 5.0, const double& norm = 1.0)
{
  CvMatrix magnitude, angle, magnNorm;
  cv::cartToPolar(flowX, flowY, magnitude, angle, true);
  cv::threshold(magnitude, magnitude, thresh, thresh, cv::THRESH_TRUNC);
  cv::normalize(magnitude, magnNorm, 0.0, norm, cv::NORM_MINMAX);
  angle *= 180.0 / (360.0 * 255.0);
  // build hsv image
  CvMatrix _hsv[3], hsv, hsv8;
  _hsv[0] = angle;
  _hsv[1] = CvMatrix::ones(angle.size(), CV_32F);
  _hsv[2] = magnNorm;
  cv::merge(_hsv, 3, hsv);
  hsv.convertTo(hsv8, CV_8U, 255.0);
  cv::cvtColor(hsv8, bgr, cv::COLOR_HSV2BGR);
}

void
flow2colour(const CvMatrix& flow, CvMatrix& bgr, const double& thresh = 5.0,
            const double& norm = 1.0)
{
  CvMatrix flowParts[2];
  cv::split(flow, flowParts);
  flow2colour(flowParts[0], flowParts[1], bgr, thresh, norm);
}
}  // namespace event_model

#endif  // IMAGE_CONVERSION_H
