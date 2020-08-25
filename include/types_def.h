#ifndef TYPES_DEF_H
#define TYPES_DEF_H

#include <Eigen/Core>
#include <unsupported/Eigen/CXX11/Tensor>

template<typename T, int rows=Eigen::Dynamic, int cols=Eigen::Dynamic>
using mtx=typename Eigen::Matrix<T, rows, cols>;
template<typename T, int rows=Eigen::Dynamic>
using vec=typename Eigen::Matrix<T, rows, 1>;

template<typename T, int rows>
using array=typename Eigen::array<T, rows>;
template<typename T, int rank>
using tensor=typename Eigen::Tensor<T, rank>;

#endif // TYPES_DEF_H
