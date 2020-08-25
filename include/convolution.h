#ifndef CONVOLUTION_H
#define CONVOLUTION_H

#include <algorithm>
#include <cassert>
#include <cmath>
#include <Eigen/Core>
#include <unsupported/Eigen/CXX11/Tensor>
#include <utility>

namespace convolution
{

template<typename T>
struct OffsetOp
{
  T
  operator()(
    const T& val) const
  {
    return (val-1)>>1;
  }

  T
  operator()(
    const T& val1, const T& val2) const
  {
    return val1+val2-1;
  }
};

template<typename T>
struct AddOffsetOp
{
  T
  operator()(
    const T& val1, const T& val2) const
  {
    return val1+val2;
  }
};

template<typename T>
struct SubOffsetOp
{
  T
  operator()(
    const T& val1, const T& val2) const
  {
    return val1-val2-1;
  }
};

template<typename Container, typename Functor>
Container
transform(
  const Container& c, Functor&& op)
{
  Container container;
  std::transform(c.begin(), c.end(), container.begin(), op);
  return container;
}

template<typename Container, typename Functor>
Container
transform(
  const Container& c1, const Container& c2, Functor&& op)
{
  Container container;
  std::transform(c1.begin(), c1.end(), c2.begin(), container.begin(), op);
  return container;
}

template<typename T, int N, typename U=T>
class Convolution
{

public:

  typedef Eigen::Tensor<T, N> tensorT;
  typedef Eigen::Tensor<U, N> tensorU;
  typedef Eigen::array<Eigen::Index, N> arrayN;

private:

  const arrayN dim_;
  const arrayN kdim_;
  const arrayN offset_;

protected:

  const U lambda_;

  tensorT tensor_;
  tensorU ttensor_;

public:

  Convolution(
    const arrayN& dim,
    const arrayN& kdim,
    const U& lambda=U(2.0*M_PI),
    const U& tref=U(0.0)):
    dim_(dim),
    kdim_(kdim),
    offset_(std::move(transform(kdim_, OffsetOp<Eigen::Index>()))),
    lambda_(lambda),
    tensor_(std::move(transform(dim_, kdim_, OffsetOp<Eigen::Index>()))), ttensor_(dim_)
  {
    setZero();
    setZeroTime(tref);
  }

  U
  lambda(void) const
  {
    return lambda_;
  }
  const Eigen::TensorRef<const tensorT>
  tensor(void) const
  {
    return tensor_;
  }
  T
  tensor(
    const arrayN& ind) const
  {
    return tensor_(std::move(transform(ind, offset_, AddOffsetOp<Eigen::Index>())));
  }
  const Eigen::TensorRef<const tensorU>
  ttensor(void) const
  {
    return ttensor_;
  }
  U
  ttensor(
    const arrayN& ind) const
  {
    return ttensor_(ind);
  }

  void
  conv(
    const arrayN& ind,
    const U& ts,
    const T& pval,
    const tensorU& kernel)
  {
    // found no better way of doing this using Eigen::Tensor
    const U expts=std::exp(-lambda_*(ts-ttensor(ind)));
    const arrayN dim(std::move(transform(ind, kdim_, AddOffsetOp<Eigen::Index>())));
    arrayN indTemp;
    iterate(expts, pval, kernel, 0, dim, ind, indTemp);
    ttensor_(ind)=ts;
  }

  // found no efficient way of implementing update() using Eigen::Tensor

  void
  setZero(void)
  {
    tensor_.setZero();
  }
  void
  setZeroTime(
    const U& ti=U(0.0))
  {
    ttensor_.setConstant(ti);
  }

protected:

  void
  iterate(
    const U& expts,
    const T& pval,
    const tensorU& kernel,
    const int d,
    const arrayN& dim, const arrayN& iind,
    arrayN& ind)
  {
    if(d>=N)
    {
      const arrayN kind(std::move(transform(dim, ind, SubOffsetOp<Eigen::Index>())));
      tensor_(ind)*=expts;
      tensor_(ind)+=pval*kernel(kind);
      return;
    }

    for(ind.at(d)=iind.at(d); ind.at(d)<dim.at(d); ++ind.at(d))
    {
      iterate(expts, pval, kernel, d+1, dim, iind, ind);
    }
  }

private:

};

template<typename T, typename U=T>
class Convolution2
{

public:
  
  typedef Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> mtxT;
  typedef Eigen::Matrix<U, Eigen::Dynamic, Eigen::Dynamic> mtxU;

private:

  const int width_, height_;
  const int kwidth_, kheight_;
  const int xoffset_, yoffset_;

protected:

  const U lambda_;

  mtxT mtx_;
  mtxU tmtx_;

public:

  Convolution2(
    const int width, const int height,
    const int kwidth, const int kheight,
    const U& lambda=U(2.0*M_PI),
    const U& tref=U(0.0)):
    width_(width), height_(height),
    kwidth_(kwidth), kheight_(kheight),
    xoffset_((kwidth_-1)>>1), yoffset_((kheight_-1)>>1),
    lambda_(lambda),
    mtx_(width_+kwidth_-1, height_+kheight_-1), tmtx_(width_, height_)
  {
    setZero();
    setZeroTime(tref);
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
  const Eigen::Ref<const mtxT>
  mtx(void) const
  {
    return mtx_.block(xoffset_, yoffset_, width(), height());
  }
  T
  mtx(
    const int x, const int y) const
  {
    assert(-xoffset_<=x && x<width()+xoffset_ && -yoffset_<=y && y<height()+yoffset_);
    return mtx_(x+xoffset_, y+yoffset_);
  }
  const Eigen::Ref<const mtxU>
  tmtx(void) const
  {
    return tmtx_;
  }
  U
  tmtx(
    const int x, const int y) const
  {
    assert(0<=x && x<width() && 0<=y && y<height());
    return tmtx_(x, y);
  }

  void
  conv(
    const int x, const int y,
    const U& ts, const U& lambda,
    const T& pval,
    const mtxU& kernel)
  {
    assert(0<=x && x<width() && 0<=y && y<height());
    assert(kernel.rows()<=kwidth_ && kernel.cols()<=kheight_);

    mtx_.template block(x, y, kernel.rows(), kernel.cols())*=std::exp(-lambda*(ts-tmtx(x, y)));
    mtx_.template block(x, y, kernel.rows(), kernel.cols())+=pval*kernel.reverse();
    tmtx_(x, y)=ts;
  }

  void
  conv(
    const int x, const int y,
    const U& ts,
    const T& pval,
    const mtxU& kernel)
  {
    conv(x, y, ts, lambda(), pval, kernel);
  }

  void
  update(
    const U& ts)
  {
    const mtxU tsimg((-lambda()*(ts-tmtx().array())).exp());

    /* core image block */
    mtx_.block(xoffset_, yoffset_, width(), height()).array()*=tsimg.array();

    /* reflective time update */

    /* top left corner image block */
    mtx_.topLeftCorner(xoffset_, yoffset_).array()*=tsimg.block(1, 1, xoffset_, yoffset_).reverse().array();
    /* top right corner image block */
    mtx_.topRightCorner(xoffset_, yoffset_).array()*=tsimg.block(1, height()-1-yoffset_, xoffset_, yoffset_).reverse().array();
    /* bottom left corner image block */
    mtx_.bottomLeftCorner(xoffset_, yoffset_).array()*=tsimg.block(width()-1-xoffset_, 1, xoffset_, yoffset_).reverse().array();
    /* bottom right corner image block */
    mtx_.bottomRightCorner(xoffset_, yoffset_).array()*=tsimg.block(width()-1-xoffset_, height()-1-yoffset_, xoffset_, yoffset_).reverse().array();
    
    /* top image rows */
    mtx_.topRows(xoffset_).middleCols(yoffset_, height()).array()*=tsimg.middleRows(xoffset_, xoffset_).colwise().reverse().array();
    /* bottom image rows */
    mtx_.bottomRows(xoffset_).middleCols(yoffset_, height()).array()*=tsimg.middleRows(width()-1-xoffset_, xoffset_).colwise().reverse().array();
    /* left image columns */
    mtx_.leftCols(yoffset_).middleRows(xoffset_, width()).array()*=tsimg.middleCols(yoffset_, yoffset_).rowwise().reverse().array();
    /* right image columns */
    mtx_.rightCols(yoffset_).middleRows(xoffset_, width()).array()*=tsimg.middleCols(height()-1-yoffset_, yoffset_).rowwise().reverse().array();
    
    setZeroTime(ts);
  }

  void
  setZero(void)
  {
    mtx_.setZero();
  }
  void
  setZeroTime(
    const U& ti=U(0.0))
  {
    tmtx_.setConstant(ti);
  }

protected:

private:

};

}

#endif // CONVOLUTION_H
