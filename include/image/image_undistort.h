#ifndef IMAGE_UNDISTORT_H
#define IMAGE_UNDISTORT_H

#include <Eigen/Core>
#include <fstream>
#include <opencv2/core/eigen.hpp>
#include <opencv2/opencv.hpp>
#include <string>

template<typename T>
void
initUndistort(
  const int width, const int height,
  const cv::Mat& camParams, const cv::Mat& distCoeffs,
  cv::Mat& map)
{
  const int size=width*height;

  cv::Mat points(size, 1, CV_MAKETYPE(cv::DataDepth<T>::value, 2));
  for(int y=0; y<height; ++y)
  {
    const int yy=y*width;
    for(int x=0; x<width; ++x)
    {
      points.at<cv::Vec<T, 2> >(yy+x)=cv::Vec<T, 2>(x, y);
    }
  }

  cv::Mat mapPoints(size, 1, CV_MAKETYPE(cv::DataDepth<T>::value, 2));
  cv::undistortPoints(points, mapPoints, camParams, distCoeffs, cv::noArray(), camParams);

  map.create(height, width, CV_MAKETYPE(cv::DataDepth<T>::value, 2));
  for(int y=0; y<height; ++y)
  {
    const int yy=y*width;
    for(int x=0; x<width; ++x)
    {
      map.at<cv::Vec<T, 2> >(y, x)=mapPoints.at<cv::Vec<T, 2> >(yy+x);
    }
  }
}

template<typename T>
bool
loadCamParams(
  const std::string& fname,
  int& width, int& height,
  cv::Mat& camParams, cv::Mat& distCoeffs)
{
  std::ifstream fin(fname.c_str());
  if(!fin.is_open())
  {
    return false;
  }

  camParams.create(3, 3, CV_MAKETYPE(cv::DataDepth<T>::value, 1));
  distCoeffs.create(5, 1, CV_MAKETYPE(cv::DataDepth<T>::value, 1));
  fin>>camParams.at<T>(0, 0)>>camParams.at<T>(1, 1)>>camParams.at<T>(0, 2)>>camParams.at<T>(1, 2);
  camParams.at<T>(0, 1)=0.0, camParams.at<T>(1, 0)=0.0;
  camParams.at<T>(2, 0)=0.0, camParams.at<T>(2, 1)=0.0, camParams.at<T>(2, 2)=1.0;
  fin>>distCoeffs.at<T>(0)>>distCoeffs.at<T>(1)>>distCoeffs.at<T>(2)>>distCoeffs.at<T>(3)>>distCoeffs.at<T>(4);
  fin>>width>>height;

  return true;
}

template<typename T>
bool
loadCamParams(
  const std::string& fname,
  int& width, int& height,
  cv::Mat& camParams)
{
  cv::Mat distCoeffs;
  return loadCamParams<T>(fname, width, height, camParams, distCoeffs);
}

template<typename T>
bool
loadCamParams(
  const std::string& fname,
  int& width, int& height,
  mtx<T, 3, 3>& camParams)
{
  cv::Mat camParamsCV;
  const bool status=loadCamParams<T>(fname, width, height, camParamsCV);
  cv::cv2eigen(camParamsCV, camParams);
  return status;
}

#endif // IMAGE_UNDISTORT_H
