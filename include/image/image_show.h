#ifndef IMAGE_SHOW_H
#define IMAGE_SHOW_H

#include <cassert>
#include <Eigen/Core>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/opencv.hpp>
#include <string>

#include "image/image_conversion.h"
#include "types_def.h"

template<typename T, typename D>
void
showGray(
	const D& c,
	const int width, const int height,
	const std::string& imgName="Image of Points")
{
	assert(c.rows()==2);

	cv::Mat img;
	eigen2cv<T>(c, width, height, img);
	cv2gray(img);
	cv::imshow(imgName, img);
}

template<typename T, typename D1, typename D2>
void
showGray(
	const D1& c,
	const D2& w,
	const int width, const int height,
	const std::string& imgName="Image of Points")
{
	assert(c.rows()==2);
	assert(c.cols()==w.size());

	cv::Mat img;
	eigen2cv<T>(c, w, width, height, img);
	cv2gray(img);
	cv::imshow(imgName, img);
}

template<typename D>
void
showGray(
 	const D& c,
 	const std::string& imgName="Image of Points")
{
	typedef typename D::Scalar T;
 	const cv::Mat img(c.cols(), c.rows(),
		CV_MAKETYPE(cv::DataDepth<T>::value, 1), const_cast<T*>(c.data()));
 	cv::Mat imgGray;
	cv2gray(img, imgGray);
	cv::imshow(imgName, imgGray);
}

#endif // IMAGE_SHOW_H
