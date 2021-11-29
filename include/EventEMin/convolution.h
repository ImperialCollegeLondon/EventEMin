#ifndef EVENT_EMIN_CONVOLUTION_H
#define EVENT_EMIN_CONVOLUTION_H

#include <algorithm>
#include <cassert>
#include <cmath>
#include <utility>

#include "EventEMin/types_def.h"

namespace EventEMin
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
 public:
  typedef Tensor<U, N> Kernel;

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
       const Kernel& kernel)
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

template <typename T, typename U>
class Convolution<T, 2, U>
{
 public:
  typedef Matrix<U> Kernel;

 private:
  const Array<Index, 2> dim_;
  const Array<Index, 2> kdim_;
  const Array<Index, 2> offset_;

 protected:
  const U lambda_;

  Matrix<T> val_;
  Matrix<U> ts_;

 public:
  Convolution(const Array<Index, 2>& dim, const Array<Index, 2>& kdim,
              const U& lambda = U(2.0 * M_PI), const U& tsRef = U(0.0))
      : dim_(dim),
        kdim_(kdim),
        offset_(std::move(transform(kdim_, OffsetOp<Index>()))),
        lambda_(lambda),
        val_(dim_[0] + kdim_[0] - 1, dim_[1] + kdim_[1] - 1),
        ts_(dim_[0], dim_[1])
  {
    reset(tsRef);
  }

  int
  width(void) const
  {
    return dim_[0];
  }
  int
  height(void) const
  {
    return dim_[1];
  }
  U
  lambda(void) const
  {
    return lambda_;
  }
  const Ref<const Matrix<T>>
  val(void) const
  {
    return val_.block(offset_[0], offset_[1], width(), height());
  }
  T
  val(const Array<Index, 2>& ind) const
  {
    assert(-offset_[0] <= ind[0] && ind[0] < width() + offset_[0]);
    assert(-offset_[1] <= ind[1] && ind[1] < height() + offset_[1]);
    return val_(ind[0] + offset_[0], ind[1] + offset_[1]);
  }
  const Ref<const Matrix<U>>
  ts(void) const
  {
    return ts_;
  }
  U
  ts(const Array<Index, 2>& ind) const
  {
    assert(0 <= ind[0] && ind[0] < width());
    assert(0 <= ind[1] && ind[1] < height());
    return ts_(ind[0], ind[1]);
  }

  void
  conv(const Array<Index, 2>& ind, const U& ts, const U& lambda, const T& val,
       const Kernel& kernel)
  {
    assert(0 <= ind[0] && ind[0] < width());
    assert(0 <= ind[1] && ind[1] < height());
    assert(kernel.rows() <= kdim_[0]);
    assert(kernel.cols() <= kdim_[1]);

    val_.template block(ind[0], ind[1], kernel.rows(), kernel.cols()) *=
        std::exp(-lambda * (ts - ts_(ind[0], ind[1])));
    val_.template block(ind[0], ind[1], kernel.rows(), kernel.cols()) +=
        val * kernel.reverse();
    ts_(ind[0], ind[1]) = ts;
  }

  void
  conv(const Array<Index, 2>& ind, const U& ts, const T& val,
       const Kernel& kernel)
  {
    conv(ind, ts, lambda(), val, kernel);
  }

  void
  update(const U& ts)
  {
    const Matrix<U> tsimg((-lambda() * (ts - ts_.array())).exp());

    // core image block
    val_.block(offset_[0], offset_[1], width(), height()).array() *=
        tsimg.array();

    /* reflective time update */

    // top left corner image block
    val_.topLeftCorner(offset_[0], offset_[1]).array() *=
        tsimg.block(1, 1, offset_[0], offset_[1]).reverse().array();
    // top right corner image block
    val_.topRightCorner(offset_[0], offset_[1]).array() *=
        tsimg.block(1, height() - 1 - offset_[1], offset_[0], offset_[1])
            .reverse()
            .array();
    // bottom left corner image block
    val_.bottomLeftCorner(offset_[0], offset_[1]).array() *=
        tsimg.block(width() - 1 - offset_[0], 1, offset_[0], offset_[1])
            .reverse()
            .array();
    // bottom right corner image block
    val_.bottomRightCorner(offset_[0], offset_[1]).array() *=
        tsimg
            .block(width() - 1 - offset_[0], height() - 1 - offset_[1],
                   offset_[0], offset_[1])
            .reverse()
            .array();

    // top image rows
    val_.topRows(offset_[0]).middleCols(offset_[1], height()).array() *=
        tsimg.middleRows(offset_[0], offset_[0]).colwise().reverse().array();
    // bottom image rows
    val_.bottomRows(offset_[0]).middleCols(offset_[1], height()).array() *=
        tsimg.middleRows(width() - 1 - offset_[0], offset_[0])
            .colwise()
            .reverse()
            .array();
    // left image columns
    val_.leftCols(offset_[1]).middleRows(offset_[0], width()).array() *=
        tsimg.middleCols(offset_[1], offset_[1]).rowwise().reverse().array();
    // right image columns
    val_.rightCols(offset_[1]).middleRows(offset_[0], width()).array() *=
        tsimg.middleCols(height() - 1 - offset_[1], offset_[1])
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
}  // namespace EventEMin

#endif  // EVENT_EMIN_CONVOLUTION_H
