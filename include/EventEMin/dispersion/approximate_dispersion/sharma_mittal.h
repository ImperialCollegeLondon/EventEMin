#ifndef EVENT_EMIN_APPROXIMATE_SHARMA_MITTAL_H
#define EVENT_EMIN_APPROXIMATE_SHARMA_MITTAL_H

#include "EventEMin/dispersion/approximate_dispersion/dispersion.h"

namespace EventEMin
{
namespace batch
{
namespace approximate
{
template <typename M>
class SharmaMittal;
}

template <typename M>
struct DispersionTraits<approximate::Dispersion<approximate::SharmaMittal<M> > >
{
  typedef M Model;
  typedef typename Model::T T;
};

namespace approximate
{
template <typename M>
class SharmaMittal : public Dispersion<SharmaMittal<M> >
{
 public:
  typedef M Model;
  typedef typename Model::T T;

  enum
  {
    NVars = Model::NVars,
    NDims = Model::NDims
  };

 private:
  const T alpha_, beta_, gamma_;
  const T scale_;

 public:
  SharmaMittal(const T& dimScaleMax, const T& offset = T(5.0),
               const T& lambda = T(0.0), const int ksize = 3,
               const T& alpha = T(2.0), const T& beta = T(0.5))
      : Dispersion<SharmaMittal<Model> >(dimScaleMax, offset, lambda, ksize),
        alpha_(alpha),
        beta_(beta),
        gamma_((T(1.0) - beta_) / (T(1.0) - alpha_)),
        scale_(T(1.0) / (T(1.0) - beta_))
  {
    assert(T(0.0) < alpha_ && alpha_ != T(1.0));
    assert(beta_ != T(1.0));
    if constexpr (NDims == 2)
    {
      this->kernel_ = this->kernel_.array().pow(alpha_);
    }
    else
    {
      this->kernel_ = this->kernel_.pow(alpha_);
    }
  }

  T
  alpha(void) const
  {
    return alpha_;
  }
  T
  beta(void) const
  {
    return beta_;
  }
  T
  gamma(void) const
  {
    return gamma_;
  }
  T
  scale(void) const
  {
    return scale_;
  }

  template <typename U>
  U
  score(const Matrix<U>& c, const Convolution<U, NDims, T>& conv) const
  {
    return scale() * pow(this->interpolate(c, conv) / this->nPoints(), gamma());
  }
};
}  // namespace approximate
}  // namespace batch

template <typename M>
using ApproximateSharmaMittal = batch::approximate::SharmaMittal<M>;
}  // namespace EventEMin

#endif  // EVENT_EMIN_APPROX_SHARMA_MITTAL_H
