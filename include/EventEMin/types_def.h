#ifndef EVENT_EMIN_TYPES_DEF_H
#define EVENT_EMIN_TYPES_DEF_H

#include <Eigen/Core>
#include <opencv2/core/eigen.hpp>
#include <opencv2/opencv.hpp>
#include <unsupported/Eigen/AutoDiff>
#include <unsupported/Eigen/CXX11/Tensor>
#include <vector>

namespace EventEMin
{
// std definitions
template <typename T>
using StdVector = typename std::vector<T, Eigen::aligned_allocator<T> >;

// OpenCV definitions
typedef cv::Mat CvMatrix;
template <typename T, int N>
using CvVector = cv::Vec<T, N>;

#define CV_TYPE(T, N) CV_MAKETYPE(cv::DataDepth<T>::value, N)

// Eigen definitions
constexpr int Dynamic = Eigen::Dynamic;

typedef Eigen::Index Index;

template <typename T, int Rows>
using Array = typename Eigen::array<T, Rows>;
template <typename T, int Rows = Dynamic, int Cols = Dynamic>
using Matrix = typename Eigen::Matrix<T, Rows, Cols>;
template <typename T, int Cols = Dynamic>
using RowVector = typename Eigen::Matrix<T, 1, Cols>;
template <typename T, int rank>
using Tensor = typename Eigen::Tensor<T, rank>;
template <typename T, int Rows = Dynamic>
using Vector = typename Eigen::Matrix<T, Rows, 1>;

template <typename T>
using DenseBase = typename Eigen::DenseBase<T>;
template <typename T>
using MatrixBase = typename Eigen::MatrixBase<T>;

template <typename T>
using Map = typename Eigen::Map<T>;

template <typename T>
using Ref = Eigen::Ref<T>;
template <typename T>
using TensorRef = typename Eigen::TensorRef<T>;
}  // namespace EventEMin

#endif  // EVENT_EMIN_TYPES_DEF_H
