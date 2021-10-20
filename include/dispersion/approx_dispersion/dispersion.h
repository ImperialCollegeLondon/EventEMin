#ifndef APPROX_DISPERSION_IMPL_H
#define APPROX_DISPERSION_IMPL_H

#include "convolution.h"
#include "dispersion/dispersion.h"
#include "event/transform.h"

namespace event_model
{
template <typename Derived>
class ApproxDispersion : public DispersionBase<ApproxDispersion<Derived> >
{
 public:
  typedef typename DispersionTraits<ApproxDispersion<Derived> >::Model Model;
  typedef typename DispersionTraits<ApproxDispersion<Derived> >::T T;

  enum
  {
    NVars = Model::NVars,
    NDims = Model::NDims
  };

 private:
  const T dimScaleMax_;
  Array<Index, NDims> cMin_, cMax_, dim_;

 protected:
  Array<Index, NDims> kdim_;
  const T halfOffset_, offset_;
  const T lambda_;
  Tensor<T, NDims> kernel_;

 public:
  ApproxDispersion(const T& dimScaleMax, const T& offset = T(0.0),
                   const T& lambda = T(0.0), const int ksize = 3)
      : DispersionBase<ApproxDispersion<Derived> >(),
        dimScaleMax_(dimScaleMax),
        halfOffset_(T(0.5) * offset),
        offset_(offset),
        lambda_(lambda)
  {
    assert(T(0.0) <= offset_);
    assert(0 < ksize && ksize % 2 == 1);
    assert(T(0.0) <= lambda_);

    for (int d = 0; d < NDims; ++d)
    {
      cMin_[d] = 0;
    }

    kdim_.fill(ksize);
    kernel_.resize(kdim_);
    kernel::gaussKernel<T, NDims>(ksize, kernel_);
    const Tensor<T, 0> gaussKernelSum = kernel_.sum();
    kernel_ /= kernel_.constant(gaussKernelSum(0));
  }

  T
  dimScale(const int d) const
  {
    assert(0 <= d && d < NDims);
    return static_cast<T>(cMax_[d]);
  }
  int
  dim(const int d) const
  {
    assert(0 <= d && d < NDims);
    return dim_[d];
  }
  const Array<Index, NDims>&
  dim(void) const
  {
    return dim_;
  }
  T
  halfOffset(void) const
  {
    return halfOffset_;
  }
  T
  offset(void) const
  {
    return offset_;
  }
  const Array<Index, NDims>&
  kdim(void) const
  {
    return kdim_;
  }
  T
  lambda(void) const
  {
    return lambda_;
  }
  const TensorRef<const Tensor<T, NDims> >
  kernel(void) const
  {
    return kernel_;
  }

  template <typename U>
  U
  compute(const Matrix<U>& c) const
  {
    Convolution<U, NDims, T> conv(dim(), kdim(), lambda(), this->tsRef());
    add(c, conv);
    return this->underlying().score(c, conv);
  }

  void
  computeDimScale(void)
  {
    const Vector<T, NDims> cLimDiff(this->c().rowwise().maxCoeff() -
                                    this->c().rowwise().minCoeff());
    const T limDiffMax = cLimDiff.maxCoeff();
    for (int d = 0; d < NDims; ++d)
    {
      dim_[d] = std::round(dimScaleMax_ * cLimDiff(d) / limDiffMax);
      cMax_[d] = dim_[d] - 1;
    }
  }

 protected:
  template <typename U>
  void
  add(const Matrix<U>& c, Convolution<U, NDims, T>& conv) const
  {
    const computeBoundaries<Index, NDims, U> computeBoundaries;
    Vector<U, NDims> val;
    Array<Index, NDims> cl, lcMin, lcMax;
    for (int k = 0; k < this->nPoints(); ++k)
    {
      computeBoundaries(c.col(k), cMin_, cMax_, lcMin, lcMax);
      addIterate<U>(this->ts(k), c.col(k), lcMin, lcMax, val, cl, conv);
    }
    // no conv update
  }

  template <typename U>
  U
  interpolate(const Matrix<U>& c, const Convolution<U, NDims, T>& conv) const
  {
    const computeBoundaries<Index, NDims, U> computeBoundaries;
    Vector<U, NDims> val;
    Array<Index, NDims> cl, lcMin, lcMax;
    U f = U(0.0);
    for (int k = 0; k < this->nPoints(); ++k)
    {
      computeBoundaries(c.col(k), cMin_, cMax_, lcMin, lcMax);
      interpolateIterate<U>(this->ts(k), c.col(k), lcMin, lcMax, conv, val, cl,
                            f);
    }
    return f;
  }

 private:
  template <typename U>
  void
  addIterate(const T& ts, const Ref<const Vector<U, NDims> >& c,
             const Array<Index, NDims>& lcMin, const Array<Index, NDims>& lcMax,
             Vector<U, NDims>& val, Array<Index, NDims>& cl,
             Convolution<U, NDims, T>& conv) const
  {
    addIterate<U>(ts, 0, c, lcMin, lcMax, val, cl, conv);
  }
  template <typename U>
  void
  addIterate(const T& ts, const int d, const Ref<const Vector<U, NDims> >& c,
             const Array<Index, NDims>& lcMin, const Array<Index, NDims>& lcMax,
             Vector<U, NDims>& val, Array<Index, NDims>& cl,
             Convolution<U, NDims, T>& conv) const
  {
    if (d >= NDims)
    {
      conv.conv(cl, ts, val.prod(), kernel());
      return;
    }

    for (cl[d] = lcMin[d]; cl[d] <= lcMax[d]; ++cl[d])
    {
      val(d) = T(1.0) - ((c(d) > cl[d]) ? c(d) - cl[d] : cl[d] - c(d));
      addIterate<U>(ts, d + 1, c, lcMin, lcMax, val, cl, conv);
    }
  }

  template <typename U>
  void
  interpolateIterate(const T& ts, const Ref<const Vector<U, NDims> >& c,
                     const Array<Index, NDims>& lcMin,
                     const Array<Index, NDims>& lcMax,
                     const Convolution<U, NDims, T>& conv,
                     Vector<U, NDims>& val, Array<Index, NDims>& cl, U& f) const
  {
    interpolateIterate<U>(ts, 0, c, lcMin, lcMax, conv, val, cl, f);
  }
  template <typename U>
  void
  interpolateIterate(const T& ts, const int d,
                     const Ref<const Vector<U, NDims> >& c,
                     const Array<Index, NDims>& lcMin,
                     const Array<Index, NDims>& lcMax,
                     const Convolution<U, NDims, T>& conv,
                     Vector<U, NDims>& val, Array<Index, NDims>& cl, U& f) const
  {
    if (d >= NDims)
    {
      f += val.prod() * std::exp(-lambda() * (this->tsEnd() - ts)) *
           conv.val(cl);
      return;
    }

    for (cl[d] = lcMin[d]; cl[d] <= lcMax[d]; ++cl[d])
    {
      val(d) = T(1.0) - ((c(d) > cl[d]) ? c(d) - cl[d] : cl[d] - c(d));
      interpolateIterate<U>(ts, d + 1, c, lcMin, lcMax, conv, val, cl, f);
    }
  }

  Derived&
  underlying(void)
  {
    return static_cast<Derived&>(*this);
  }
  const Derived&
  underlying(void) const
  {
    return static_cast<const Derived&>(*this);
  }
};
}

#endif  // APPROX_DISPERSION_IMPL_H
