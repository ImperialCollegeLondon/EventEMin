#ifndef CONVOLUTION_H
#define CONVOLUTION_H

#include <algorithm>
#include <cassert>
#include <cmath>
#include <utility>

#include "types_def.h"

namespace event_model
{
template <typename T>
struct OffsetOp
{
  T
  operator()(const T& val) const
  {
    return (val - 1) >> 1;
  }

  T
  operator()(const T& val1, const T& val2) const
  {
    return val1 + val2 - 1;
  }
};

template <typename T>
struct AddOffsetOp
{
  T
  operator()(const T& val1, const T& val2) const
  {
    return val1 + val2;
  }
};

template <typename T>
struct SubOffsetOp
{
  T
  operator()(const T& val1, const T& val2) const
  {
    return val1 - val2 - 1;
  }
};

template <typename Container, typename Functor>
Container
transform(const Container& c, Functor&& op)
{
  Container container;
  std::transform(c.begin(), c.end(), container.begin(), op);
  return container;
}

template <typename Container, typename Functor>
Container
transform(const Container& c1, const Container& c2, Functor&& op)
{
  Container container;
  std::transform(c1.begin(), c1.end(), c2.begin(), container.begin(), op);
  return container;
}

template <typename T, int N, typename U = T>
class Convolution
{
 private:
  const Array<Index, N> dim_;
  const Array<Index, N> kdim_;
  const Array<Index, N> offset_;

 protected:
  const U lambda_;

  Tensor<T, N> val_;
  Tensor<U, N> ts_;

 public:
  Convolution(const Array<Index, N>& dim, const Array<Index, N>& kdim,
              const U& lambda = U(2.0 * M_PI), const U& tsRef = U(0.0))
      : dim_(dim),
        kdim_(kdim),
        offset_(std::move(transform(kdim_, OffsetOp<Index>()))),
        lambda_(lambda),
        val_(std::move(transform(dim_, kdim_, OffsetOp<Index>()))),
        ts_(dim_)
  {
    reset(tsRef);
  }

  U
  lambda(void) const
  {
    return lambda_;
  }
  const TensorRef<const Tensor<T, N>>
  val(void) const
  {
    return val_;
  }
  T
  val(const Array<Index, N>& ind) const
  {
    return val_(std::move(transform(ind, offset_, AddOffsetOp<Index>())));
  }
  const TensorRef<const Tensor<U, N>>
  ts(void) const
  {
    return ts_;
  }
  U
  ts(const Array<Index, N>& ind) const
  {
    return ts_(ind);
  }

  void
  conv(const Array<Index, N>& ind, const U& ts, const T& val,
       const Tensor<U, N>& kernel)
  {
    // found no better way of doing this using Tensor
    const U expts = std::exp(-lambda_ * (ts - ts_(ind)));
    const Array<Index, N> dim(
        std::move(transform(ind, kdim_, AddOffsetOp<Index>())));
    Array<Index, N> indTemp;
    iterate(expts, val, kernel, 0, dim, ind, indTemp);
    ts_(ind) = ts;
  }

  // found no efficient way of implementing update() using Tensor

  void
  reset(const U& ti = U(0.0))
  {
    val_.setZero();
    ts_.setConstant(ti);
  }

 protected:
  void
  iterate(const U& expts, const T& val, const Tensor<U, N>& kernel, const int d,
          const Array<Index, N>& dim, const Array<Index, N>& iind,
          Array<Index, N>& ind)
  {
    if (d >= N)
    {
      const Array<Index, N> kind(
          std::move(transform(dim, ind, SubOffsetOp<Index>())));
      val_(ind) *= expts;
      val_(ind) += val * kernel(kind);
      return;
    }

    for (ind[d] = iind[d]; ind[d] < dim[d]; ++ind[d])
    {
      iterate(expts, val, kernel, d + 1, dim, iind, ind);
    }
  }

 private:
};

template <typename T, typename U = T>
class Convolution2
{
 private:
  const int width_, height_;
  const int kwidth_, kheight_;
  const int xoffset_, yoffset_;

 protected:
  const U lambda_;

  Matrix<T> val_;
  Matrix<U> ts_;

 public:
  Convolution2(const int width, const int height, const int kwidth,
               const int kheight, const U& lambda = U(2.0 * M_PI),
               const U& tsRef = U(0.0))
      : width_(width),
        height_(height),
        kwidth_(kwidth),
        kheight_(kheight),
        xoffset_((kwidth_ - 1) >> 1),
        yoffset_((kheight_ - 1) >> 1),
        lambda_(lambda),
        val_(width_ + kwidth_ - 1, height_ + kheight_ - 1),
        ts_(width_, height_)
  {
    reset(tsRef);
  }

  int
  width(void) const
  {
    return width_;
  }
  int
  height(void) const
  {
    return height_;
  }
  U
  lambda(void) const
  {
    return lambda_;
  }
  const Ref<const Matrix<T>>
  val(void) const
  {
    return val_.block(xoffset_, yoffset_, width(), height());
  }
  T
  val(const int x, const int y) const
  {
    assert(-xoffset_ <= x && x < width() + xoffset_);
    assert(-yoffset_ <= y && y < height() + yoffset_);
    return val_(x + xoffset_, y + yoffset_);
  }
  const Ref<const Matrix<U>>
  ts(void) const
  {
    return ts_;
  }
  U
  ts(const int x, const int y) const
  {
    assert(0 <= x && x < width());
    assert(0 <= y && y < height());
    return ts_(x, y);
  }

  void
  conv(const int x, const int y, const U& ts, const U& lambda, const T& val,
       const Matrix<U>& kernel)
  {
    assert(0 <= x && x < width());
    assert(0 <= y && y < height());
    assert(kernel.rows() <= kwidth_);
    assert(kernel.cols() <= kheight_);

    val_.template block(x, y, kernel.rows(), kernel.cols()) *=
        std::exp(-lambda * (ts - ts_(x, y)));
    val_.template block(x, y, kernel.rows(), kernel.cols()) +=
        val * kernel.reverse();
    ts_(x, y) = ts;
  }

  void
  conv(const int x, const int y, const U& ts, const T& val,
       const Matrix<U>& kernel)
  {
    conv(x, y, ts, lambda(), val, kernel);
  }

  void
  update(const U& ts)
  {
    const Matrix<U> tsimg((-lambda() * (ts - ts_.array())).exp());

    // core image block
    val_.block(xoffset_, yoffset_, width(), height()).array() *= tsimg.array();

    /* reflective time update */

    // top left corner image block
    val_.topLeftCorner(xoffset_, yoffset_).array() *=
        tsimg.block(1, 1, xoffset_, yoffset_).reverse().array();
    // top right corner image block
    val_.topRightCorner(xoffset_, yoffset_).array() *=
        tsimg.block(1, height() - 1 - yoffset_, xoffset_, yoffset_)
            .reverse()
            .array();
    // bottom left corner image block
    val_.bottomLeftCorner(xoffset_, yoffset_).array() *=
        tsimg.block(width() - 1 - xoffset_, 1, xoffset_, yoffset_)
            .reverse()
            .array();
    // bottom right corner image block
    val_.bottomRightCorner(xoffset_, yoffset_).array() *=
        tsimg
            .block(width() - 1 - xoffset_, height() - 1 - yoffset_, xoffset_,
                   yoffset_)
            .reverse()
            .array();

    // top image rows
    val_.topRows(xoffset_).middleCols(yoffset_, height()).array() *=
        tsimg.middleRows(xoffset_, xoffset_).colwise().reverse().array();
    // bottom image rows
    val_.bottomRows(xoffset_).middleCols(yoffset_, height()).array() *=
        tsimg.middleRows(width() - 1 - xoffset_, xoffset_)
            .colwise()
            .reverse()
            .array();
    // left image columns
    val_.leftCols(yoffset_).middleRows(xoffset_, width()).array() *=
        tsimg.middleCols(yoffset_, yoffset_).rowwise().reverse().array();
    // right image columns
    val_.rightCols(yoffset_).middleRows(xoffset_, width()).array() *=
        tsimg.middleCols(height() - 1 - yoffset_, yoffset_)
            .rowwise()
            .reverse()
            .array();

    resetTime(ts);
  }

  void
  reset(const U& ti = U(0.0))
  {
    val_.setZero();
    resetTime(ti);
  }
  void
  resetTime(const U& ti = U(0.0))
  {
    ts_.setConstant(ti);
  }
};
}  // namespace event_model

#endif  // CONVOLUTION_H
