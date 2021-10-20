#ifndef IMAGE_UNDISTORT_H
#define IMAGE_UNDISTORT_H

#include "types_def.h"

namespace event_model
{
template <typename T>
void
initUndistort(const int width, const int height, const CvMatrix& camParams,
              const CvMatrix& distCoeffs, CvMatrix& map)
{
  const int size = width * height;

  CvMatrix points(size, 1, CV_TYPE(T, 2));
  for (int y = 0; y < height; ++y)
  {
    const int yy = y * width;
    for (int x = 0; x < width; ++x)
    {
      points.at<CvVector<T, 2> >(yy + x) = CvVector<T, 2>(x, y);
    }
  }

  CvMatrix mapPoints(size, 1, CV_TYPE(T, 2));
  cv::undistortPoints(points, mapPoints, camParams, distCoeffs, cv::noArray(),
                      camParams);

  map.create(height, width, CV_TYPE(T, 2));
  for (int y = 0; y < height; ++y)
  {
    const int yy = y * width;
    for (int x = 0; x < width; ++x)
    {
      map.at<CvVector<T, 2> >(y, x) = mapPoints.at<CvVector<T, 2> >(yy + x);
    }
  }
}
}  // namespace event_model

#endif  // IMAGE_UNDISTORT_H
