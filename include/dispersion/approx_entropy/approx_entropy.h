#ifndef APPROX_ENTROPY_H
#define APPROX_ENTROPY_H

#include "convolution.h"
#include "dispersion/dispersion.h"

namespace event_model
{

namespace dispersion
{

namespace approx_entropy
{

// partial specialisation only works for structs not functions
template<typename I, int N, typename T>
struct computeBoundaries
{

  void
  operator()(
    const Eigen::Ref<const vec<T, N> >& c,
    const array<I, N>& cMin, const array<I, N>& cMax,
    array<I, N>& lMin, array<I, N>& lMax)
  {
    for(I d=0; d<N; ++d)
    {
      lMin.at(d)=std::max(static_cast<I>(std::floor(c(d).value())), cMin.at(d));
      lMax.at(d)=std::min(static_cast<I>(std::ceil(c(d).value())), cMax.at(d));
    }
  }

};
template<typename I, int N>
struct computeBoundaries<I, N, double>
{

  void
  operator()(
    const Eigen::Ref<const vec<double, N> >& c,
    const array<I, N>& cMin, const array<I, N>& cMax,
    array<I, N>& lMin, array<I, N>& lMax)
  {
    for(I d=0; d<N; ++d)
    {
      lMin.at(d)=std::max(static_cast<I>(std::floor(c(d))), cMin.at(d));
      lMax.at(d)=std::min(static_cast<I>(std::ceil(c(d))), cMax.at(d));
    }
  }

};
template<typename I, int N>
struct computeBoundaries<I, N, float>
{

  void
  operator()(
    const Eigen::Ref<const vec<float, N> >& c,
    const array<I, N>& cMin, const array<I, N>& cMax,
    array<I, N>& lMin, array<I, N>& lMax)
  {
    for(I d=0; d<N; ++d)
    {
      lMin.at(d)=std::max(static_cast<I>(std::floor(c(d))), cMin.at(d));
      lMax.at(d)=std::min(static_cast<I>(std::ceil(c(d))), cMax.at(d));
    }
  }

};

}

template<typename Derived>
class ApproxEntropy : public Dispersion<ApproxEntropy<Derived> >
{

public:

  typedef typename DispersionTraits<ApproxEntropy<Derived> >::Model Model;
  typedef typename DispersionTraits<ApproxEntropy<Derived> >::T T;

private:

  array<Eigen::Index, Model::ndims()> cMin_, cMax_;

protected:

  array<Eigen::Index, Model::ndims()> dim_, kdim_;
  const T halfOffset_, offset_;
  const T lambda_;
  tensor<T, Model::ndims()> kernel_;

public:

  ApproxEntropy(
    const array<Eigen::Index, Model::ndims()>& dim,
    const T& offset=T(0.0),
    const T& lambda=T(0.0),
    const int ksize=3):
    Dispersion<ApproxEntropy<Derived> >(),
    dim_(dim),
    halfOffset_(T(0.5)*offset), offset_(offset),
    lambda_(lambda)
  {
    assert(T(0.0)<=offset_);
    assert(0<ksize && ksize%2==1 && 0.0<=lambda_);

    for(int d=0; d<this->ndims(); ++d)
    {
      cMin_.at(d)=0;
      cMax_.at(d)=dim_.at(d)-1;
    }
    kdim_.fill(ksize);
    kernel_.resize(kdim_);
    kernel::gaussKernel<T, this->ndims()>(ksize, kernel_);
    const tensor<T, 0> kernelSum=kernel_.sum();
    kernel_=kernel_/kernel_.constant(kernelSum(0));
  }

  int
  dim(
    const int d) const
  {
    assert(0<=d && d<this->ndims());
    return dim_.at(d);
  }
  const array<Eigen::Index, Model::ndims()>&
  dim(void) const
  {
    return dim_;
  }
  T
  dimScale(
    const int d) const
  {
    assert(0<=d && d<this->ndims());
    return static_cast<T>(dim(d)-1);
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
  const array<Eigen::Index, Model::ndims()>&
  kdim(void) const
  {
    return kdim_;
  }
  T
  lambda(void) const
  {
    return lambda_;
  }
  const Eigen::TensorRef<const tensor<T, Model::ndims()> >
  kernel(void) const
  {
    return kernel_;
  }

  template<typename U>
  U
  compute(
    const mtx<U>& c) const
  {
    convolution::Convolution<U, this->ndims(), T> conv(dim(),
      kdim(), lambda(), this->tsref());

    add(c, conv);
    return this->underlying().fFinal(c, conv);
  }

protected:

  template<typename U>
  void
  add(
    const mtx<U>& c,
    convolution::Convolution<U, Model::ndims(), T>& conv) const
  {
    vec<U> val(this->ndims());
    array<Eigen::Index, this->ndims()> cl, lcMin, lcMax;
    for(int k=0; k<this->npoints(); ++k)
    {
      approx_entropy::computeBoundaries<Eigen::Index, this->ndims(), U>()(c.col(k), cMin_, cMax_, lcMin, lcMax);
      addIterate<U>(this->ts(k), c.col(k), lcMin, lcMax, val, cl, conv);
    }
    // no conv update
  }

  template<typename U>
  U
  interpolate(
    const mtx<U>& c,
    const convolution::Convolution<U, Model::ndims(), T>& conv) const
  {
    vec<U> val(this->ndims());
    array<Eigen::Index, this->ndims()> cl, lcMin, lcMax;
    U f=U(0.0);
    for(int k=0; k<this->npoints(); ++k)
    {
      approx_entropy::computeBoundaries<Eigen::Index, this->ndims(), U>()(c.col(k), cMin_, cMax_, lcMin, lcMax);
      interpolateIterate<U>(this->ts(k), c.col(k), lcMin, lcMax, conv, val, cl, f);
    }
    return f;
  }

private:

  template<typename U>
  void
  addIterate(
    const T& ts,
    const Eigen::Ref<const vec<U, Model::ndims()> >& c,
    const array<Eigen::Index, Model::ndims()>& lcMin, const array<Eigen::Index, Model::ndims()>& lcMax,
    vec<U>& val,
    array<Eigen::Index, Model::ndims()>& cl,
    convolution::Convolution<U, Model::ndims(), T>& conv) const
  {
    addIterate<U>(ts, T(1.0), 0, c, lcMin, lcMax, val, cl, conv);
  }
  template<typename U>
  void
  addIterate(
    const T& ts, const T& s,
    const Eigen::Ref<const vec<U, Model::ndims()> >& c,
    const array<Eigen::Index, Model::ndims()>& lcMin, const array<Eigen::Index, Model::ndims()>& lcMax,
    vec<U>& val,
    array<Eigen::Index, Model::ndims()>& cl,
    convolution::Convolution<U, Model::ndims(), T>& conv) const
  {
    addIterate<U>(ts, s, 0, c, lcMin, lcMax, val, cl, conv);
  }
  template<typename U>
  void
  addIterate(
    const T& ts, const T& s, const int d,
    const Eigen::Ref<const vec<U, Model::ndims()> >& c,
    const array<Eigen::Index, Model::ndims()>& lcMin, const array<Eigen::Index, Model::ndims()>& lcMax,
    vec<U>& val,
    array<Eigen::Index, Model::ndims()>& cl,
    convolution::Convolution<U, Model::ndims(), T>& conv) const
  {
    if(d>=this->ndims())
    {
      conv.conv(cl, ts, s*val.prod(), kernel());
      return;
    }

    for(cl.at(d)=lcMin.at(d); cl.at(d)<=lcMax.at(d); ++cl.at(d))
    {
      val(d)=1.0-((c(d)>cl.at(d)) ? c(d)-cl.at(d) : cl.at(d)-c(d));
      addIterate<U>(ts, s, d+1, c, lcMin, lcMax, val, cl, conv);
    }
  }

  template<typename U>
  void
  interpolateIterate(
    const T& ts,
    const Eigen::Ref<const vec<U, Model::ndims()> >& c,
    const array<Eigen::Index, Model::ndims()>& lcMin, const array<Eigen::Index, Model::ndims()>& lcMax,
    const convolution::Convolution<U, Model::ndims(), T>& conv,
    vec<U>& val,
    array<Eigen::Index, Model::ndims()>& cl,
    U& f) const
  {
    interpolateIterate<U>(ts, T(1.0), 0, c, lcMin, lcMax, conv, val, cl, f);
  }
  template<typename U>
  void
  interpolateIterate(
    const T& ts, const T& s,
    const Eigen::Ref<const vec<U, Model::ndims()> >& c,
    const array<Eigen::Index, Model::ndims()>& lcMin, const array<Eigen::Index, Model::ndims()>& lcMax,
    const convolution::Convolution<U, Model::ndims(), T>& conv,
    vec<U>& val,
    array<Eigen::Index, Model::ndims()>& cl,
    U& f) const
  {
    interpolateIterate<U>(ts, s, 0, c, lcMin, lcMax, conv, val, cl, f);
  }
  template<typename U>
  void
  interpolateIterate(
    const T& ts, const T& s, const int d,
    const Eigen::Ref<const vec<U, Model::ndims()> >& c,
    const array<Eigen::Index, Model::ndims()>& lcMin, const array<Eigen::Index, Model::ndims()>& lcMax,
    const convolution::Convolution<U, Model::ndims(), T>& conv,
    vec<U>& val,
    array<Eigen::Index, Model::ndims()>& cl,
    U& f) const
  {
    if(d>=this->ndims())
    {
      f+=s*val.prod()*std::exp(-lambda()*(this->tsend()-ts))*conv.tensor(cl);
      return;
    }

    for(cl.at(d)=lcMin.at(d); cl.at(d)<=lcMax.at(d); ++cl.at(d))
    {
      val(d)=1.0-((c(d)>cl.at(d)) ? c(d)-cl.at(d) : cl.at(d)-c(d));
      interpolateIterate<U>(ts, s, d+1, c, lcMin, lcMax, conv, val, cl, f);
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

}

#endif // APPROX_ENTROPY_H
